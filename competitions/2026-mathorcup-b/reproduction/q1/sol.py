from __future__ import annotations
import json
from pathlib import Path


from bisect import bisect_left, bisect_right
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
from math import isfinite, sqrt
from numbers import Number
from typing import Any, Iterable, Mapping, Sequence

A_DEFAULT_YEARS = (2021, 2022, 2023, 2024, 2025)
A_ALLOWED_OFFSETS = tuple(range(-15, 16))
A_DIRECTIONS = {"high", "low", "two_sided"}


def _a_get(mapping: Mapping, key: int, default=None):
    return mapping.get(key, mapping.get(str(key), default))


def _a_finite(value: Any) -> bool:

    if isinstance(value, (bool, complex)):return False
    if isinstance(value, (int, Fraction)):return True
    return isinstance(value, Number) and isfinite(value)


def _a_fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def a_midrank_fraction(sorted_values: Sequence[Number], value: Number) -> Fraction:

    if not sorted_values:
        raise ValueError("midrank reference must be nonempty")
    left, right = bisect_left(sorted_values, value), bisect_right(sorted_values, value)
    return Fraction(left + right, 2 * len(sorted_values))


def _a_direction(q: Fraction, direction: str) -> Fraction:
    if direction == "high":
        return q
    if direction == "low":
        return 1 - q
    if direction == "two_sided":
        return abs(2 * q - 1)
    raise ValueError(f"unknown direction: {direction}")


def _a_mean(values: Sequence[Fraction]) -> Fraction:
    return sum(values, Fraction(0)) / len(values)


def _a_prepare(feature: Mapping[str, Any], data: Mapping[str, Any], years: tuple[int, ...],
             min_endpoints: int) -> dict[str, Any]:
    factor_id = str(feature["id"])
    direction = feature["direction"]
    if direction not in A_DIRECTIONS:
        raise ValueError(f"{factor_id}: direction must be one of {sorted(A_DIRECTIONS)}")
    windows = tuple(int(w) for w in feature["window_hours"])
    if not windows or len(set(windows)) != len(windows) or any(w <= 0 for w in windows):
        raise ValueError(f"{factor_id}: window_hours must be distinct positive integers")
    raw_target = data.get("target", {})
    target = {w: _a_get(raw_target, w) for w in windows}
    per_year = data.get("years", {})
    endpoints, n_valid = {}, {}
    for year in years:
        raw_endpoints = _a_get(per_year, year, {})
        accepted = {}
        for offset in A_ALLOWED_OFFSETS:
            row = _a_get(raw_endpoints, offset, {})
            vals = {w: _a_get(row, w) for w in windows}
            if all(_a_finite(v) for v in vals.values()):
                accepted[offset] = vals
        n_valid[year] = len(accepted)
        if len(accepted) >= min_endpoints:
            endpoints[year] = accepted
    return {"id": factor_id, "direction": direction, "windows": windows,
            "target": target, "target_valid": all(_a_finite(v) for v in target.values()),
            "endpoints": endpoints, "n_valid": n_valid}


def _a_score_event(values: Mapping[int, Number], windows: tuple[int, ...], direction: str,
                 sorted_by_year: Mapping[int, Mapping[int, list]],
                 reference_years: tuple[int, ...]) -> dict[str, Any]:
    q = {w: _a_mean([a_midrank_fraction(sorted_by_year[y][w], values[w])
                   for y in reference_years]) for w in windows}
    a = {w: _a_direction(q[w], direction) for w in windows}
    score_t = max(a.values())
    return {"q": q, "a": a, "T": score_t,
            "winning_windows": [w for w in windows if a[w] == score_t]}


def _a_event_json(event: Mapping[str, Any]) -> dict[str, Any]:
    winners = event["winning_windows"]
    return {"q_by_window": {str(w): float(q) for w, q in event["q"].items()},
            "a_by_window": {str(w): float(a) for w, a in event["a"].items()},
            "T": float(event["T"]), "T_exact": _a_fraction_text(event["T"]),
            "winning_windows": winners, "winning_window_tie": len(winners) > 1}


def _a_score_factor(feature: Mapping[str, Any], prepared: Mapping[str, Any],
                  years: tuple[int, ...], min_endpoints: int, min_years: int
                  ) -> tuple[dict[str, Any], Fraction | None]:
    windows, endpoints = prepared["windows"], prepared["endpoints"]
    valid_years = tuple(sorted(endpoints))
    factor = {"id": prepared["id"], "direction": prepared["direction"],
              "window_hours": list(windows), "eligible_for_top5": bool(feature.get("eligible_for_top5", False)),
              "score_S": None, "score_S_exact": None, "rank_eligible": False,
              "support": {"requested_years": list(years), "valid_years": list(valid_years),
                          "n_valid_endpoints_by_year": {str(y): prepared["n_valid"][y] for y in years},
                          "min_endpoints_per_year": min_endpoints, "min_years": min_years,
                          "endpoint_rule": "all applicable windows jointly valid"},
              "windows": {str(w): {"target_value": prepared["target"][w]} for w in windows},
              "configurations": [], "calibration_year_scores": {}, "constant_flags": {}}
    if not prepared["target_valid"]:
        factor.update(status="not_computable", reason="target missing at least one predefined window")
        return factor, None
    if len(valid_years) < min_years:
        factor.update(status="not_computable", reason=f"only {len(valid_years)} valid historical years; require {min_years}")
        return factor, None

    sorted_by_year = {y: {w: sorted(row[w] for row in endpoints[y].values())
                          for w in windows} for y in valid_years}
    background_constant = {
        w: len({row[w] for y in valid_years for row in endpoints[y].values()}) == 1
        for w in windows}
    all_values_constant = {
        w: background_constant[w] and prepared["target"][w] == sorted_by_year[valid_years[0]][w][0]
        for w in windows}
    window_q, window_a, winning_counts = defaultdict(list), defaultdict(list), defaultdict(int)
    scores_by_h = defaultdict(list)
    all_configs_calibration_T_constant, all_configs_target_tied = True, True
    for h in valid_years:
        for reference in combinations([y for y in valid_years if y != h], 3):
            target_event = _a_score_event(prepared["target"], windows, prepared["direction"],
                                        sorted_by_year, reference)
            pseudo_events = [(offset, _a_score_event(vals, windows, prepared["direction"],
                                                  sorted_by_year, reference))
                             for offset, vals in sorted(endpoints[h].items())]
            pseudo_t = sorted(event["T"] for _, event in pseudo_events)
            config_s = a_midrank_fraction(pseudo_t, target_event["T"])
            scores_by_h[h].append(config_s)
            constant_t = len(set(pseudo_t)) == 1
            target_tied = constant_t and pseudo_t[0] == target_event["T"]
            all_configs_calibration_T_constant &= constant_t
            all_configs_target_tied &= target_tied
            for w in windows:
                window_q[w].append(target_event["q"][w])
                window_a[w].append(target_event["a"][w])
                winning_counts[w] += int(w in target_event["winning_windows"])
            factor["configurations"].append({
                "calibration_year": h, "reference_years": list(reference),
                "target": _a_event_json(target_event),
                "pseudo_events": [{"offset": offset, **_a_event_json(event)}
                                  for offset, event in pseudo_events],
                "n_calibration_endpoints": len(pseudo_events),
                "score_S": float(config_s), "score_S_exact": _a_fraction_text(config_s),
                "calibration_T_all_tied": constant_t,
                "target_tied_with_every_calibration_T": target_tied})
    by_year = {h: _a_mean(scores_by_h[h]) for h in valid_years}
    score_s = _a_mean(list(by_year.values()))
    factor["calibration_year_scores"] = {
        str(h): {"score_S": float(by_year[h]), "score_S_exact": _a_fraction_text(by_year[h]),
                 "configuration_count": len(scores_by_h[h])} for h in valid_years}
    n_config = len(factor["configurations"])
    for w in windows:
        factor["windows"][str(w)].update({
            "q_mean": float(_a_mean(window_q[w])), "q_min": float(min(window_q[w])),
            "q_max": float(max(window_q[w])), "a_mean": float(_a_mean(window_a[w])),
            "a_min": float(min(window_a[w])), "a_max": float(max(window_a[w])),
            "winning_configuration_count": winning_counts[w],
            "winning_configuration_fraction": winning_counts[w] / n_config,
            "background_constant": background_constant[w],
            "background_and_target_constant": all_values_constant[w]})
    factor.update(score_S=float(score_s), score_S_exact=_a_fraction_text(score_s),
                  configuration_count=n_config,
                  constant_flags={"all_windows_background_constant": all(background_constant.values()),
                                  "all_windows_background_and_target_constant": all(all_values_constant.values()),
                                  "all_configurations_calibration_T_all_tied": all_configs_calibration_T_constant,
                                  "all_configurations_target_and_calibration_T_all_tied": all_configs_target_tied})
    if all_configs_target_tied:
        factor.update(status="computable_uninformative", reason="target and all calibration T values tie in every configuration")
    else:
        factor.update(status="computable", reason=None, rank_eligible=True)
    return factor, score_s


def _a_average_ranks(values: Sequence[Number]) -> list[float]:
    ordered = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    start = 0
    while start < len(values):
        stop = start + 1
        while stop < len(values) and values[ordered[stop]] == values[ordered[start]]:
            stop += 1
        rank = (start + 1 + stop) / 2
        for i in ordered[start:stop]:
            ranks[i] = rank
        start = stop
    return ranks


def a_spearman_midrank(x: Sequence[Number], y: Sequence[Number]) -> float | None:

    if len(x) != len(y) or not x:
        raise ValueError("Spearman arrays must be nonempty and equal length")
    rx, ry = _a_average_ranks(x), _a_average_ranks(y)
    mean_x, mean_y = sum(rx) / len(rx), sum(ry) / len(ry)
    dx, dy = [r - mean_x for r in rx], [r - mean_y for r in ry]
    ssx, ssy = sum(v * v for v in dx), sum(v * v for v in dy)
    if not ssx or not ssy:
        return None
    return max(-1.0, min(1.0, sum(a * b for a, b in zip(dx, dy)) / sqrt(ssx * ssy)))


def _a_explicitly_same_information(a: Mapping[str, Any], b: Mapping[str, Any]) -> bool:
    family = a.get("physical_family")
    if not family or family != b.get("physical_family"):
        return False
    explicit_pair = (b["id"] in a.get("duplicate_with", []) or
                     a["id"] in b.get("duplicate_with", []))
    return explicit_pair or (a.get("duplicate_information") is True and
                             b.get("duplicate_information") is True)


def _a_redundancy_pairs(features: Mapping[str, Mapping[str, Any]], prepared: Mapping[str, Any],
                      factor_results: Mapping[str, Any], min_endpoints: int,
                      min_years: int, threshold: float) -> list[dict[str, Any]]:
    pairs = []
    for a_id, b_id in combinations(sorted(features), 2):
        a, b = features[a_id], features[b_id]
        if not _a_explicitly_same_information(a, b):
            continue
        wa, wb = int(a.get("representative_hours", 168)), int(b.get("representative_hours", 168))
        result = {"factor_a": a_id, "factor_b": b_id, "physical_family": a.get("physical_family"),
                  "representative_hours_a": wa, "representative_hours_b": wb,
                  "threshold": threshold, "rho": None, "abs_rho": None, "merge_eligible": False,
                  "n_paired_endpoints_by_year": {}, "common_valid_years": []}
        if wa not in prepared[a_id]["windows"] or wb not in prepared[b_id]["windows"]:
            result.update(status="uncertain", reason="predeclared representative window unavailable")
            pairs.append(result)
            continue
        if factor_results[a_id]["score_S"] is None or factor_results[b_id]["score_S"] is None:
            result.update(status="uncertain", reason="one or both factors not computable")
            pairs.append(result)
            continue
        x, y = [], []
        a_years, b_years = prepared[a_id]["endpoints"], prepared[b_id]["endpoints"]
        for year in sorted(set(a_years) & set(b_years)):
            offsets = sorted(set(a_years[year]) & set(b_years[year]))
            result["n_paired_endpoints_by_year"][str(year)] = len(offsets)
            if len(offsets) >= min_endpoints:
                result["common_valid_years"].append(year)
                x.extend(a_years[year][o][wa] for o in offsets)
                y.extend(b_years[year][o][wb] for o in offsets)
        result["paired_endpoints_used"] = len(x)
        if len(result["common_valid_years"]) < min_years:
            result.update(status="uncertain", reason="insufficient jointly valid historical years/endpoints")
        else:
            rho = a_spearman_midrank(x, y)
            result.update(rho=rho, abs_rho=None if rho is None else abs(rho))
            if rho is None:
                result.update(status="uncertain", reason="constant representative-window margin")
            else:
                result.update(status="evaluated", reason=None, merge_eligible=abs(rho) >= threshold)
        pairs.append(result)
    return pairs


def _a_rank_rows(ids: Iterable[str], scores: Mapping[str, Fraction | None],
               factor_results: Mapping[str, Any]) -> list[dict[str, Any]]:
    ordered = sorted(ids, key=lambda factor_id: (-scores[factor_id], factor_id))
    rows, start = [], 0
    while start < len(ordered):
        stop = start + 1
        while stop < len(ordered) and scores[ordered[stop]] == scores[ordered[start]]:
            stop += 1
        for factor_id in ordered[start:stop]:
            rows.append({"id": factor_id, "score_S": float(scores[factor_id]),
                         "score_S_exact": _a_fraction_text(scores[factor_id]),
                         "rank_min": start + 1, "rank_max": stop,
                         "score_tie_size": stop - start,
                         "eligible_for_top5": factor_results[factor_id]["eligible_for_top5"]})
        start = stop
    return rows


def summarize_one(features: Sequence[Mapping[str, Any]], values_by_factor: Mapping[str, Any], *,
                  years: Sequence[int] = A_DEFAULT_YEARS, min_endpoints: int = 25,
                  min_years: int = 4, redundancy_threshold: float = .90,
                  top_n: int = 5) -> dict[str, Any]:


    years = tuple(sorted(int(y) for y in years))
    if len(set(years)) != len(years):
        raise ValueError("historical years must be distinct")
    if min_years < 4 or not 1 <= min_endpoints <= 31 or top_n < 1:
        raise ValueError("require min_years >= 4, 1 <= min_endpoints <= 31 and top_n >= 1")
    if not 0 <= redundancy_threshold <= 1:
        raise ValueError("redundancy_threshold must be in [0,1]")
    specs = {str(feature["id"]): feature for feature in features}
    if len(specs) != len(features):
        raise ValueError("feature ids must be unique")
    prepared = {factor_id: _a_prepare(spec, values_by_factor.get(factor_id, {}), years, min_endpoints)
                for factor_id, spec in specs.items()}
    factors, scores = {}, {}
    for factor_id, spec in specs.items():
        factors[factor_id], scores[factor_id] = _a_score_factor(
            spec, prepared[factor_id], years, min_endpoints, min_years)
    ranking_all = _a_rank_rows([factor_id for factor_id in specs if factors[factor_id]["rank_eligible"]],
                             scores, factors)
    pairs = _a_redundancy_pairs(specs, prepared, factors, min_endpoints, min_years, redundancy_threshold)
    pair_lookup = {frozenset((pair["factor_a"], pair["factor_b"])): pair
                   for pair in pairs if pair["merge_eligible"]}
    retained, suppressed = [], []
    for row in ranking_all:
        factor_id = row["id"]
        if not row["eligible_for_top5"]:
            continue
        representatives = [r for r in retained if frozenset((factor_id, r)) in pair_lookup]
        if representatives:
            rep = representatives[0]
            pair = pair_lookup[frozenset((factor_id, rep))]
            suppressed.append({"id": factor_id, "representative_id": rep,
                               "reason": "predeclared duplicate information and supported representative-window Spearman threshold",
                               "score_tied_with_representative": scores[factor_id] == scores[rep],
                               "rho": pair["rho"], "common_valid_years": pair["common_valid_years"]})
        else:
            retained.append(factor_id)
    candidate_ranking = _a_rank_rows(retained, scores, factors)
    cutoff = scores[retained[top_n - 1]] if len(retained) >= top_n else None
    top_candidates = [row["id"] for row in candidate_ranking
                      if cutoff is None or scores[row["id"]] >= cutoff]
    boundary = [row["id"] for row in candidate_ranking
                if row["rank_min"] <= top_n < row["rank_max"]]
    return {"method": "A_equal_year_cross_calibrated_window_max", "score_interpretation": "descriptive cross-year calibrated historical anomaly position; not a p-value or disaster probability",
            "settings": {"years": list(years), "min_endpoints": min_endpoints,
                         "min_years": min_years, "redundancy_threshold": redundancy_threshold,
                         "top_n": top_n, "tie_rule": "exact midrank on caller-supplied source-precision values; exact rational ranks",
                         "redundancy_rule": "same explicitly duplicated physical information; pooled endpoint Spearman on predeclared representative windows; sufficient paired support in each of >=4 years"},
            "factors": factors, "ranking_all": ranking_all,
            "ranking_exploratory": [r for r in ranking_all if not r["eligible_for_top5"]],
            "redundancy": {"pairs": pairs, "suppressed": suppressed,
                           "representative_tie_rule": "equal score: lexicographic factor id; score tie remains explicitly recorded"},
            "candidate_ranking": candidate_ranking, "top_candidates": top_candidates,
            "top_core_candidates": [r["id"] for r in candidate_ranking if r["rank_max"] <= top_n],
            "boundary_candidates": boundary,
            "top_cutoff_score_S": None if cutoff is None else float(cutoff),
            "top_candidate_count": len(top_candidates), "available_candidate_count": len(retained),
            "fewer_than_requested_candidates": len(retained) < top_n}


def summarize_selection_frequency(named_runs: Mapping[str, Mapping[str, Any]],
                                  features: Sequence[Mapping[str, Any]] | None = None
                                  ) -> dict[str, Any]:


    ids = sorted({str(f["id"]) for f in features} if features is not None else
                 {factor_id for run in named_runs.values() for factor_id in run.get("factors", {})})
    summaries = {}
    for factor_id in ids:
        records, eligible = [], False
        for name, run in named_runs.items():
            factor = run.get("factors", {}).get(factor_id, {})
            eligible |= bool(factor.get("eligible_for_top5", False))
            computable = factor.get("score_S") is not None
            all_rank = next((r for r in run.get("ranking_all", []) if r["id"] == factor_id), None)
            candidate_rank = next((r for r in run.get("candidate_ranking", []) if r["id"] == factor_id), None)
            records.append({"run": str(name), "computable": computable,
                            "omitted_year_was_valid": factor.get("omitted_year_was_valid"),
                            "no_op": bool(factor.get("no_op", False)),
                            "rankable": bool(factor.get("rank_eligible", False)),
                            "selected": factor_id in run.get("top_candidates", []),
                            "score_S": factor.get("score_S"),
                            "rank_min": None if candidate_rank is None else candidate_rank["rank_min"],
                            "rank_max": None if candidate_rank is None else candidate_rank["rank_max"],
                            "pre_redundancy_rank_min": None if all_rank is None else all_rank["rank_min"],
                            "pre_redundancy_rank_max": None if all_rank is None else all_rank["rank_max"],
                            "status": factor.get("status", "not_computable"),
                            "reason": factor.get("reason", "factor absent from this run")})
        n_computable = sum(r["computable"] for r in records)
        n_selected = sum(r["selected"] for r in records)
        perturbed = [r for r in records if not r["no_op"]]
        n_perturbed_computable = sum(r["computable"] for r in perturbed)
        n_perturbed_selected = sum(r["selected"] for r in perturbed)
        available_ranks = [r for r in records if r["rank_min"] is not None]
        all_ranks = [r for r in records if r["pre_redundancy_rank_min"] is not None]
        summaries[factor_id] = {
            "eligible_for_top5": eligible, "selected_count": n_selected,
            "computable_count": n_computable, "not_computable_count": len(records) - n_computable,
            "rankable_count": sum(r["rankable"] for r in records),
            "selection_fraction": n_selected / n_computable if eligible and n_computable else None,
            "no_op_count": len(records) - len(perturbed),
            "perturbed_selected_count": n_perturbed_selected,
            "perturbed_computable_count": n_perturbed_computable,
            "perturbed_not_computable_count": len(perturbed) - n_perturbed_computable,
            "perturbed_fraction": (n_perturbed_selected / n_perturbed_computable
                                   if eligible and n_perturbed_computable else None),
            "rank_min": min((r["rank_min"] for r in available_ranks), default=None),
            "rank_max": max((r["rank_max"] for r in available_ranks), default=None),
            "pre_redundancy_rank_min": min((r["pre_redundancy_rank_min"] for r in all_ranks), default=None),
            "pre_redundancy_rank_max": max((r["pre_redundancy_rank_max"] for r in all_ranks), default=None),
            "run_records": records}
    return summaries


def leave_one_year_out(features: Sequence[Mapping[str, Any]], values_by_factor: Mapping[str, Any], *,
                       years: Sequence[int] = A_DEFAULT_YEARS, min_endpoints: int = 25,
                       min_years: int = 4, redundancy_threshold: float = .90,
                       top_n: int = 5) -> dict[str, Any]:

    years = tuple(int(y) for y in years)
    runs = {str(omitted): summarize_one(
        features, values_by_factor, years=[y for y in years if y != omitted],
        min_endpoints=min_endpoints, min_years=min_years,
        redundancy_threshold=redundancy_threshold, top_n=top_n) for omitted in years}
    baseline_valid_years = {
        str(spec["id"]): set(_a_prepare(spec, values_by_factor.get(str(spec["id"]), {}),
                                        years, min_endpoints)["endpoints"])
        for spec in features}
    for omitted in years:
        run = runs[str(omitted)]
        run["omitted_year"] = omitted
        for factor_id, factor in run["factors"].items():
            was_valid = omitted in baseline_valid_years[factor_id]
            factor["omitted_year_was_valid"] = was_valid
            factor["no_op"] = not was_valid
    return {"runs": runs, "selection_frequency": summarize_selection_frequency(runs, features),
            "interpretation": "sensitivity to omitted historical years; public-fold and actually-perturbed denominators are distinct; deleting a previously invalid year is a factor-specific no-op; unavailable perturbed folds stay NA; not independent-event sampling or disaster prediction coverage"}


from collections import Counter
from datetime import date, timedelta
from typing import Any

import numpy as np
from scipy.stats import kendalltau, rankdata


SEED = 20261024
WINDOWS_H = (24, 72, 168, 360, 720)
CLIM_YEARS = (2021, 2022, 2023, 2024, 2025)
N_PERM = 2000
N_BOOT = 1000
N_PERM_BOOT = 500
MIN_REFERENCE = 100


ORIGINAL_FACTORS = (
    ("precip_sum", "累计降水", 3),
    ("precip_max", "最大小时雨强", 3),
    ("rain_hours", "有雨小时数", 2),
    ("api", "前期降水指数API", 3),
    ("t2_mean", "2 m 气温", 3),
    ("t2_max", "2 m 日最高气温", 2),
    ("pdd", "站点正积温PDD", 3),
    ("warm_rate", "升温速率", 3),
    ("rn_mean", "净辐射Rn", 3),
    ("rsd_mean", "入射短波Rsd", 2),
    ("rld_mean", "入射长波Rld", 2),
    ("vwc01", "0.1 m 体积含水量", 3),
    ("vwc02", "0.2 m 体积含水量", 2),
    ("soilt01", "0.1 m 土温", 2),
    ("g01", "0.1 m 土壤热通量G", 1),
    ("dthdz_night", "夜间位温梯度", 1),
    ("rig_day", "日间梯度理查森数Ri_g中位数", 1),
    ("zeta_med", "稳定度参数ζ中位数", 1),
    ("h_mean", "感热H", 2),
    ("le_mean", "潜热LE", 2),
    ("p_mean", "气压", 1),
    ("rh2_mean", "2 m 相对湿度", 2),
    ("ws10_mean", "10 m 风速", 1),
    ("ustar_mean", "摩擦速度u*", 1),
)
MECHANISM = {key: mechanism for key, _, mechanism in ORIGINAL_FACTORS}
ORIGINAL_NAMES = {key: name for key, name, _ in ORIGINAL_FACTORS}


def _number(x: Any) -> float | int:
    if isinstance(x, (int, np.integer)):
        return int(x)
    try:
        value = float(x)
        return value if np.isfinite(value) else np.nan
    except (TypeError, ValueError):
        return np.nan


def _finite_values(values: np.ndarray) -> np.ndarray:

    values = np.asarray(values)
    if values.dtype != object:
        return np.isfinite(values)
    return np.fromiter((isinstance(v, (int, np.integer)) or
                        (isinstance(v, (float, np.floating)) and np.isfinite(v))
                       for v in values.flat), dtype=bool, count=values.size).reshape(values.shape)


def _window_value(mapping: dict, hours: int) -> float | int:
    labels = (hours, str(hours), f"{hours}h", f"{hours // 24}d")
    return next((_number(mapping[k]) for k in labels if k in mapping), np.nan)


def _features(features: list[dict] | dict) -> list[dict]:
    if isinstance(features, dict):
        if isinstance(features.get("features"), list):
            source = features["features"]
        elif isinstance(features.get("factors"), list):
            source = features["factors"]
        else:
            source = [dict(v, factor_id=k) for k, v in features.items()]
    else:
        source = features
    out = []
    for feature in source:
        key = feature.get("factor_id", feature.get("id", feature.get("key")))
        if key not in MECHANISM:

            continue
        raw_windows = feature.get("windows_h", feature.get("window_hours", feature.get("windows_hours", feature.get("windows", WINDOWS_H))))
        raw_windows = [] if raw_windows is None else raw_windows
        if isinstance(raw_windows, dict):
            raw_windows = list(raw_windows.values())
        windows = []
        for value in raw_windows:
            if isinstance(value, str) and value.endswith("d"):
                value = float(value[:-1]) * 24
            elif isinstance(value, str) and value.endswith("h"):
                value = value[:-1]
            hours = int(value)
            if hours not in WINDOWS_H:
                raise ValueError(f"Unknown original 参考方法 window: {key}/{hours}")
            if hours not in windows:
                windows.append(hours)
        out.append({**feature, "factor_id": key, "windows_h": windows})
    if len({f["factor_id"] for f in out}) != len(out):
        raise ValueError('Repeated 参考方法 factor IDs')
    return out


def bh_qvalues(p: np.ndarray) -> np.ndarray:

    p = np.asarray(p, float)
    out = np.full(p.shape, np.nan)
    finite = np.isfinite(p)
    pv = p[finite]
    if not len(pv):
        return out
    order = np.argsort(pv, kind="stable")
    q = np.minimum.accumulate((pv[order] * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1]
    values = np.empty(len(pv))
    values[order] = np.minimum(1.0, q)
    out[finite] = values
    return out


def entropy_topsis(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, str]:

    matrix = np.asarray(matrix, float)
    n, columns = matrix.shape
    if n < 2:
        return np.full(n, np.nan), np.full(columns, np.nan), "fewer_than_two_rows"
    ranges = np.ptp(matrix, axis=0)
    scaled = np.divide(matrix - matrix.min(axis=0), ranges,
                       out=np.zeros_like(matrix), where=ranges > 0)
    probability = (scaled + 1e-12) / (scaled + 1e-12).sum(axis=0)
    entropy = -(probability * np.log(probability)).sum(axis=0) / np.log(n)
    diversity = np.maximum(0.0, 1.0 - entropy)
    diversity[ranges == 0] = 0.0
    if not diversity.sum() > 1e-14:
        return np.full(n, np.nan), np.full(columns, np.nan), "all_columns_constant"
    weights = diversity / diversity.sum()
    weighted = scaled * weights
    positive = np.linalg.norm(weighted - weighted.max(axis=0), axis=1)
    negative = np.linalg.norm(weighted - weighted.min(axis=0), axis=1)
    denominator = positive + negative
    scores = np.divide(negative, denominator, out=np.full(n, np.nan), where=denominator > 0)
    return scores, weights, "ok"


def tied_ranks(scores: np.ndarray, k: int = 5) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:

    scores = np.asarray(scores, float)
    rank = np.full(scores.shape, np.nan)
    low = rank.copy()
    high = rank.copy()
    inclusion = rank.copy()
    good = np.isfinite(scores)
    if good.any():
        vals = -scores[good]
        rank[good] = rankdata(vals, method="average")
        low[good] = rankdata(vals, method="min")
        high[good] = rankdata(vals, method="max")
        size = high[good] - low[good] + 1
        inclusion[good] = np.clip((k - low[good] + 1) / size, 0, 1)
    return rank, low, high, inclusion


def _rank_all(anomaly: np.ndarray, p: np.ndarray, mechanism: np.ndarray):
    good = np.isfinite(anomaly) & np.isfinite(p) & (p > 0) & np.isfinite(mechanism)
    score = np.full(len(anomaly), np.nan)
    weights = np.full(3, np.nan)
    if good.any():
        score[good], weights, status = entropy_topsis(np.column_stack((anomaly[good], -np.log10(p[good]), mechanism[good])))
    else:
        status = "no_computable_rows"
    return score, weights, status


def _max_integer_ratio(numerator: np.ndarray, denominator: np.ndarray) -> tuple[np.ndarray, np.ndarray]:


    numerator, denominator = np.asarray(numerator, np.int64), np.asarray(denominator, np.int64)
    if numerator.shape != denominator.shape or not numerator.shape[-1]:
        raise ValueError("Nonempty equal-shaped ratio arrays are required")
    best_numerator, best_denominator = numerator[..., 0].copy(), denominator[..., 0].copy()
    for wi in range(1, numerator.shape[-1]):
        better = numerator[..., wi] * best_denominator > best_numerator * denominator[..., wi]
        best_numerator = np.where(better, numerator[..., wi], best_numerator)
        best_denominator = np.where(better, denominator[..., wi], best_denominator)
    return best_numerator, best_denominator


class _FactorCache:


    def __init__(self, key: str, windows: list[int], values: dict):
        self.key, self.windows = key, windows
        self.years = CLIM_YEARS
        seasonal = values.get("seasonal_years", {})
        if not seasonal:
            raise ValueError('seasonal_years cache is required for original 参考方法 pseudo-events')
        anchor = date(2025, 5, 31)
        self.days = [(anchor + timedelta(days=i)).strftime("%m-%d") for i in range(138)]
        self.day_index = {day: i for i, day in enumerate(self.days)}
        season_start = self.day_index["06-15"]
        season_end = self.day_index["09-30"] + 1
        self.candidates = [(yi, di) for yi in range(len(self.years)) for di in range(season_start, season_end)]
        self.candidate_year = np.array([yi for yi, _ in self.candidates], int)
        self.candidate_day = np.array([di for _, di in self.candidates], int)


        self.values = np.full((len(self.years), len(self.days), len(windows)), np.nan, dtype=object)
        for yi, year in enumerate(self.years):
            year_values = seasonal.get(year, seasonal.get(str(year), {}))
            for di, day in enumerate(self.days):
                row = year_values.get(day, {})
                for wi, hours in enumerate(windows):
                    self.values[yi, di, wi] = _window_value(row, hours)
        self.target = np.array([_window_value(values.get("target", {}), h) for h in windows], dtype=object)
        center = self.day_index["08-26"]
        self.target_ref = self.values[:, center - 15:center + 16, :]
        valid = _finite_values(self.target_ref)
        self.target_n = valid.sum(axis=1)
        self.target_less = ((self.target_ref < self.target[None, None, :]) & valid).sum(axis=1)
        self.target_equal = ((self.target_ref == self.target[None, None, :]) & valid).sum(axis=1)
        shape = (len(self.candidates), len(self.years), len(windows))
        self.pseudo_n = np.zeros(shape, np.int32)
        self.pseudo_less = np.zeros(shape, np.int32)
        self.pseudo_equal = np.zeros(shape, np.int32)
        self.pseudo_values = self.values[self.candidate_year, self.candidate_day, :]
        for ci, (yi, di) in enumerate(self.candidates):
            reference = self.values[:, di - 15:di + 16, :]
            finite = _finite_values(reference)
            target = self.pseudo_values[ci]
            self.pseudo_n[ci] = finite.sum(axis=1)
            self.pseudo_less[ci] = ((reference < target[None, None, :]) & finite).sum(axis=1)
            self.pseudo_equal[ci] = ((reference == target[None, None, :]) & finite).sum(axis=1)
        self.memo: dict[tuple[int, ...], dict] = {}
        finite_refs = _finite_values(self.target_ref)
        self.constant = bool(len(windows) and finite_refs.any() and _finite_values(self.target).all() and
                             np.all((self.target_ref == self.target[None, None, :]) | ~finite_refs))

    def distribution(self, weights: tuple[int, ...]) -> dict:
        if weights in self.memo:
            return self.memo[weights]
        weight = np.asarray(weights, int)
        target_n = (self.target_n * weight[:, None]).sum(axis=0)
        target_twice_counts = ((2 * self.target_less + self.target_equal) * weight[:, None]).sum(axis=0)
        target_q = np.divide(target_twice_counts, 2 * target_n, out=np.full(len(self.windows), np.nan), where=target_n >= MIN_REFERENCE)
        target_q[~_finite_values(self.target)] = np.nan


        target_anomaly_numerator = np.abs(target_twice_counts - target_n)
        window_a = np.divide(target_anomaly_numerator, target_n,
                             out=np.full(len(self.windows), np.nan), where=target_n >= MIN_REFERENCE)
        window_a[~_finite_values(self.target)] = np.nan

        target_complete = len(self.windows) and np.isfinite(window_a).all() and not self.constant
        if target_complete:
            target_max_num, target_max_den = _max_integer_ratio(target_anomaly_numerator, target_n)
            target_max_num, target_max_den = int(target_max_num), int(target_max_den)
            anomaly = target_max_num / target_max_den
        else:
            target_max_num = target_max_den = None
            anomaly = np.nan
        selected_weights = np.broadcast_to(weight, (len(self.candidates), len(weight))).copy()
        selected_weights[np.arange(len(self.candidates)), self.candidate_year] = 0
        pseudo_n = np.einsum("cyw,cy->cw", self.pseudo_n, selected_weights, optimize=False)
        pseudo_twice_counts = np.einsum("cyw,cy->cw", 2 * self.pseudo_less + self.pseudo_equal,
                                        selected_weights, optimize=False)
        pseudo_anomaly_numerator = np.abs(pseudo_twice_counts - pseudo_n)
        pseudo_a = np.divide(pseudo_anomaly_numerator, pseudo_n,
                             out=np.full(pseudo_n.shape, np.nan), where=pseudo_n >= MIN_REFERENCE)
        full_support = _finite_values(self.pseudo_values).all(axis=1) & np.isfinite(pseudo_a).all(axis=1)
        present_year = weight[self.candidate_year] > 0
        valid = full_support & present_year
        pseudo_max_num, pseudo_max_den = _max_integer_ratio(pseudo_anomaly_numerator[valid], pseudo_n[valid])
        pseudo_stat = pseudo_max_num / pseudo_max_den
        result = {"A": anomaly, "A_num": target_max_num, "A_den": target_max_den,
                  "q": target_q, "q_twice_count": target_twice_counts,
                  "a": window_a, "a_num": target_anomaly_numerator, "nref": target_n,
                  "pseudo": pseudo_stat, "pseudo_num": pseudo_max_num, "pseudo_den": pseudo_max_den,
                  "candidate_indices": np.flatnonzero(valid),
                  "n_pool": int(valid.sum()), "n_pool_total": int(present_year.sum())}
        self.memo[weights] = result
        return result


def _draw_p(distribution: dict, draws: int, rng: np.random.Generator) -> tuple[float, int, int]:
    if distribution["A_num"] is None or not distribution["n_pool"]:
        return np.nan, 0, 0
    picks = rng.integers(0, distribution["n_pool"], size=draws)

    greater_equal = int(np.count_nonzero(distribution["pseudo_num"][picks] * distribution["A_den"] >=
                                        distribution["A_num"] * distribution["pseudo_den"][picks]))
    return (1 + greater_equal) / (1 + draws), draws, greater_equal


def _comparison(main: np.ndarray, other: np.ndarray, ids: list[str]) -> dict:
    common = np.isfinite(main) & np.isfinite(other)
    main_r, main_lo, _, main_top = tied_ranks(np.where(common, main, np.nan))
    other_r, other_lo, _, other_top = tied_ranks(np.where(common, other, np.nan))
    tau = kendalltau(main_r[common], other_r[common]).statistic if common.sum() >= 2 else np.nan
    return {"n_common_factors": int(common.sum()), "common_factor_ids": [ids[i] for i in np.flatnonzero(common)],
            "kendall_tau_b": tau,
            "top5_fractional_overlap": float(np.nansum(np.minimum(main_top, other_top))) if common.any() else None,
            "top5_possible_overlap": int(np.sum((main_lo <= 5) & (other_lo <= 5))) if common.any() else None,
            "other_top5_possible": [ids[i] for i in np.flatnonzero(other_lo <= 5)],
            "tie_rule": "average ranks; fractional allocation at Top5 boundary; tau_b"}


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def run_ref(features: list[dict] | dict, values_by_factor: dict, seed: int = SEED, *,
              label: str = "B-aligned-corrected", n_perm: int = N_PERM,
              n_boot: int = N_BOOT, n_perm_boot: int = N_PERM_BOOT) -> dict:


    if min(n_perm, n_boot, n_perm_boot) < 1:
        raise ValueError("Draw and bootstrap counts must be positive")
    selected = _features(features)
    ids = [f["factor_id"] for f in selected]
    mechanism = np.array([MECHANISM[key] for key in ids], float)
    size = len(ids)
    rng = np.random.default_rng(seed)
    caches: list[_FactorCache | None] = []
    unavailable = {}
    for feature in selected:
        key = feature["factor_id"]
        if key not in values_by_factor or not feature["windows_h"]:
            caches.append(None)
            unavailable[key] = "no_window_cache_or_no_applicable_windows"
            continue
        try:
            caches.append(_FactorCache(key, feature["windows_h"], values_by_factor[key]))
        except ValueError as error:
            caches.append(None)
            unavailable[key] = str(error)
    full_weights = (1,) * len(CLIM_YEARS)
    main = [cache.distribution(full_weights) if cache else None for cache in caches]
    anomaly = np.array([d["A"] if d else np.nan for d in main])
    p = np.full(size, np.nan)
    draws_run = np.zeros(size, int)
    exceed = np.zeros(size, int)
    for j, distribution in enumerate(main):
        if distribution:
            p[j], draws_run[j], exceed[j] = _draw_p(distribution, n_perm, rng)
    score, entropy_weight, topsis_status = _rank_all(anomaly, p, mechanism)
    ranks, rank_min, rank_max, top_weight = tied_ranks(score)
    bh = bh_qvalues(p)
    borda_score = np.full(size, np.nan)
    good = np.isfinite(score)
    if good.any():
        borda_score[good] = -(rankdata(-anomaly[good], method="average") +
                              rankdata(p[good], method="average") +
                              rankdata(-mechanism[good], method="average"))
    borda_rank, _, _, _ = tied_ranks(borda_score)
    zscore = np.full(size, np.nan)
    for j, cache in enumerate(caches):
        if cache and 168 in cache.windows:
            wi = cache.windows.index(168)
            reference = cache.target_ref[:, :, wi].ravel()
            reference = reference[_finite_values(reference)]
            if len(reference) >= MIN_REFERENCE and _finite_values(cache.target)[wi]:
                center = reference[0]
                differences = [v - center for v in reference]
                target_difference = cache.target[wi] - center
                common_scale = max([abs(v) for v in differences] + [abs(target_difference)])
                if not common_scale > 0:
                    continue


                centered = np.array([v / common_scale for v in differences], float)
                target_centered = target_difference / common_scale
                deviation = centered.std(ddof=1)
                if deviation > 0:
                    zscore[j] = abs(target_centered - centered.mean()) / deviation
    zrank, _, _, _ = tied_ranks(zscore)

    boot_rng = np.random.default_rng(seed + 1)
    year_draws = boot_rng.choice(len(CLIM_YEARS), size=(n_boot, len(CLIM_YEARS)), replace=True)
    boot_weights = np.array([np.bincount(draw, minlength=len(CLIM_YEARS)) for draw in year_draws])
    boot_ranks = np.full((n_boot, size), np.nan)
    boot_scores = boot_ranks.copy()
    boot_p = boot_ranks.copy()
    boot_top = boot_ranks.copy()
    boot_draws = np.zeros((n_boot, size), int)
    boot_pools = np.zeros((n_boot, size), int)
    boot_status = []
    for b, weights in enumerate(boot_weights):
        weights_tuple = tuple(int(w) for w in weights)
        ab = np.full(size, np.nan)
        pb = np.full(size, np.nan)
        for j, cache in enumerate(caches):
            if cache:
                distribution = cache.distribution(weights_tuple)
                ab[j] = distribution["A"]
                boot_pools[b, j] = distribution["n_pool"]
                pb[j], boot_draws[b, j], _ = _draw_p(distribution, n_perm_boot, boot_rng)
        sb, _, status = _rank_all(ab, pb, mechanism)
        boot_scores[b], boot_p[b] = sb, pb
        boot_ranks[b], _, _, boot_top[b] = tied_ranks(sb)
        boot_status.append(status)
    rows, window_rows, bootstrap_rows = [], [], []
    for j, feature in enumerate(selected):
        key = ids[j]
        cache, distribution = caches[j], main[j]
        if cache is None:
            status = unavailable[key]
        elif cache.constant:
            status = "constant_target_and_reference_windows"
        elif not np.isfinite(anomaly[j]):
            status = "target_missing_or_some_window_reference_n_below_100"
        elif not np.isfinite(p[j]):
            status = "no_complete_pseudo_event_reference_pool"
        elif not np.isfinite(score[j]):
            status = f"topsis_{topsis_status}"
        else:
            status = "computed"
        winning = []
        if distribution and np.isfinite(anomaly[j]):
            winning = [h for wi, h in enumerate(cache.windows)
                       if int(distribution["a_num"][wi]) * distribution["A_den"] ==
                       distribution["A_num"] * int(distribution["nref"][wi])]
        main_row = {"factor_id": key, "factor": feature.get("name", ORIGINAL_NAMES[key]), "method": label,
                    "status": status, "windows_h": feature["windows_h"], "A": anomaly[j],
                    "A_exact_numerator": distribution["A_num"] if distribution else None,
                    "A_exact_denominator": distribution["A_den"] if distribution else None,
                    "value_encoding": "rank_ticks" if "rank_scale" in feature else "physical_value",
                    "rank_scale": feature.get("rank_scale"),
                    "best_windows_h": winning, "perm_p": p[j], "bh_q": bh[j],
                    "bh_family_n": int(np.isfinite(p).sum()), "mech_score": int(mechanism[j]),
                    "topsis_score": score[j], "rank": ranks[j], "rank_min": rank_min[j], "rank_max": rank_max[j],
                    "top5_fractional_allocation": top_weight[j], "borda_rank": borda_rank[j],
                    "baseline_7d_absolute_z": zscore[j], "baseline_7d_rank": zrank[j],
                    "pseudo_pool_valid": distribution["n_pool"] if distribution else 0,
                    "pseudo_pool_total": distribution["n_pool_total"] if distribution else 0,
                    "perm_draws_requested": n_perm, "perm_draws_executed": draws_run[j],
                    "perm_greater_equal": exceed[j] if draws_run[j] else None}
        finite = np.isfinite(boot_ranks[:, j])
        count = int(finite.sum())
        bounds = np.percentile(boot_ranks[finite, j], [2.5, 97.5]) if count else [np.nan, np.nan]
        bootstrap_row = {"factor_id": key, "bootstrap_requested": n_boot,
                         "bootstrap_rank_computable": count, "bootstrap_rank_not_computable": n_boot - count,
                         "rank_ci_low": bounds[0], "rank_ci_high": bounds[1],
                         "top5_fractional_count": float(np.nansum(boot_top[:, j])) if count else None,
                         "top5_frequency_computable": float(np.mean(boot_top[finite, j])) if count else None,
                         "top5_definite_count": int(np.sum(boot_top[finite, j] == 1)),
                         "top5_possible_count": int(np.sum(boot_top[finite, j] > 0)),
                         "pseudo_draws_executed": int(boot_draws[:, j].sum()),
                         "pseudo_draws_per_computable_distribution": n_perm_boot,
                         "pseudo_distributions_computable": int(np.sum(boot_draws[:, j] == n_perm_boot))}
        main_row.update(bootstrap_row)
        rows.append(main_row)
        bootstrap_rows.append(bootstrap_row)
        if cache and distribution:
            for wi, hours in enumerate(cache.windows):
                q = distribution["q"][wi]
                window_rows.append({"factor_id": key, "hours": hours, "event_value": cache.target[wi],
                                    "value_encoding": "rank_ticks" if "rank_scale" in feature else "physical_value",
                                    "rank_scale": feature.get("rank_scale"),
                                    "pooled_midrank_quantile": q, "two_sided_anomaly": distribution["a"][wi],
                                    "quantile_twice_midrank_count": distribution["q_twice_count"][wi],
                                    "anomaly_exact_numerator": distribution["a_num"][wi],
                                    "anomaly_exact_denominator": distribution["nref"][wi],
                                    "direction": "high" if q > .5 else "low" if q < .5 else "middle" if q == .5 else None,
                                    "pooled_reference_n": distribution["nref"][wi],
                                    "reference_n_by_year": dict(zip(CLIM_YEARS, cache.target_n[:, wi].tolist())),
                                    "is_best_window": hours in winning})
    comparisons = {"borda": _comparison(score, borda_score, ids),
                   "baseline_7d_absolute_z": _comparison(score, zscore, ids)}
    rows.sort(key=lambda row: (row["rank"] if np.isfinite(row["rank"]) else np.inf, ids.index(row["factor_id"])))
    result = {"factor_rows": rows, "window_rows": window_rows, "bootstrap_rows": bootstrap_rows,
              "comparisons": comparisons,
              "bootstrap_detail": {"factor_ids": ids, "year_draws": np.array(CLIM_YEARS)[year_draws],
                                   "year_multiplicity": boot_weights, "ranks": boot_ranks, "scores": boot_scores,
                                   "p_values": boot_p, "top5_fractional_allocation": boot_top,
                                   "pseudo_draws_executed": boot_draws, "pseudo_pool_sizes": boot_pools},
              "metadata": {"label": label, "source_version": "MODELING_REPORT §4.4–4.5 / final problem1.py",
                           "seed": seed, "climate_years": CLIM_YEARS, "n_original_factors": 24,
                           "n_input_original_factors": size, "n_computable_topsis_factors": int(np.isfinite(score).sum()),
                           "min_pooled_reference": MIN_REFERENCE, "n_perm_requested": n_perm,
                           "n_boot_year_resamples_executed": n_boot, "n_perm_boot_requested": n_perm_boot,
                           "formal_counts": (n_perm, n_boot, n_perm_boot) == (N_PERM, N_BOOT, N_PERM_BOOT),
                           "execution_class": "formal_counts" if (n_perm, n_boot, n_perm_boot) == (N_PERM, N_BOOT, N_PERM_BOOT) else "reduced_count_selfcheck",
                           "entropy_weights_A_neglogp_mechanism": entropy_weight, "topsis_status": topsis_status,
                           "bootstrap_topsis_status_counts": dict(Counter(boot_status)),
                           "n_boot_year_multiplicity_patterns": len({tuple(w) for w in boot_weights}),
                           "pseudo_endpoint_season": "06-15..09-30, inclusive; source year excluded",
                           "window_reference_rule": "same month-day ±15 days; same full applicable W for target and pseudo",
                           "rank_value_comparison": "exact object-integer comparisons when rank ticks are supplied; baseline centers integers before float standardization",
                           "T_and_p_comparison": "exact integer anomaly numerator/reference denominator; multiwindow max and null>=target use integer cross-products; floats only enter display and TOPSIS",
                           "weight_conflict": "Source §4.5 fits entropy weights to the event A/p matrix, contradicting assumption 3. Computational §4.5 is implemented and this is not labeled climate-only weighting.",
                           "p_interpretation": "Monte Carlo pseudo-event comparison under a seasonal sampling scheme; exchangeability is not established; not disaster probability or predictive accuracy",
                           "corrections": ["D0 cleaning and caller-defined time/support instead of original QC/interpolation",
                                           "calendar month-day matching instead of leap-year dayofyear distance",
                                           "full fixed W required; no subset-window nanmax",
                                           "eligible complete pseudo pools, exactly requested draws, unavailable p stays NA",
                                           "constant factors/degenerate TOPSIS stay NA",
                                           "average ranks and explicit Top5 tie fractions instead of ordinal tie breaks",
                                           "BH family is the reported computable factor family, not fictitious 24 valid tests"],
                           "independent_review_status": "pending; execution is not scientific approval"}}
    return _json_safe(result)


Q1_FROZEN_DICTIONARY = json.loads((Path(__file__).resolve().parent / 'data' / 'features.json').read_text(encoding='utf-8'))

def q1_sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def q1_fp(path):
    p = Path(path).resolve()
    try:
        identity=str(p.relative_to(Path.cwd().resolve()))
    except ValueError:
        identity='external_input/' + p.name
    return {'path': identity, 'bytes': p.stat().st_size, 'sha256': q1_sha(p)}

def q1_dump(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')

def q1_csv(path, rows):
    if not rows:
        Path(path).write_text('status\nno_rows\n', encoding='utf-8')
        return
    with Path(path).open('w', encoding='utf-8-sig', newline='') as f:
        fields=list(dict.fromkeys(k for row in rows for k in row))
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
        for row in rows:
            w.writerow({k: json.dumps(v, ensure_ascii=False) if isinstance(v, (dict, list)) else v for k, v in row.items()})

def q1_hour(value):
    t = dt.datetime.fromisoformat(value)
    if t.tzinfo is None or t.utcoffset() != dt.timedelta(hours=8) or t.minute or t.second or t.microsecond:
        raise ValueError('D0标签须为UTC+8整点且唯一')
    return int((t - Q1_ORIGIN).total_seconds() // 3600)

def q1_cache(arrays, features, endpoint, cutoff_hour, coverage, ledger):
    values = {}; source_features = [f for f in features if f.get('computable')]
    first = dt.date(2021, 5, 31); last = dt.date(2021, 10, 15)
    monthdays = [(first + dt.timedelta(days=d)).strftime('%m-%d') for d in range((last-first).days+1)]
    for f in source_features:
        g = f['id']; data = {'target': {}, 'years': {}, 'seasonal_years': {}}
        end = q1_hour(f'2026-08-26T{cutoff_hour:02d}:00:00+08:00')
        for h in f['window_hours']:
            value, meta = q1_window(arrays, f, end, h, endpoint, coverage)
            data['target'][str(h)] = value
            ledger.append(dict(factor_id=g,endpoint=endpoint,cutoff_hour=cutoff_hour,coverage_threshold=coverage,year=2026,monthday='08-26',window_hours=h,value_exact_raw_grid_fraction=str(value) if value is not None else None,physical_value=float(value/f['source_grid_scale']) if value is not None else None,**meta))
        for year in A_DEFAULT_YEARS:
            seasonal = {}
            for md in monthdays:
                end = q1_hour(f'{year}-{md}T{cutoff_hour:02d}:00:00+08:00'); row = {}
                for h in f['window_hours']:
                    value, meta = q1_window(arrays, f, end, h, endpoint, coverage); row[str(h)] = value
                    ledger.append(dict(factor_id=g,endpoint=endpoint,cutoff_hour=cutoff_hour,coverage_threshold=coverage,year=year,monthday=md,window_hours=h,value_exact_raw_grid_fraction=str(value) if value is not None else None,physical_value=float(value/f['source_grid_scale']) if value is not None else None,**meta))
                seasonal[md] = row
            data['seasonal_years'][str(year)] = seasonal
            data['years'][str(year)] = {str(offset): seasonal[(dt.date(year,8,26)+dt.timedelta(days=offset)).strftime('%m-%d')] for offset in range(-15,16)}
        values[g] = data
    return values

def q1_encode_exact_cache(cache, features):

    for feature in features:
        g=feature['id']
        if g not in Q1_COLUMNS:continue
        rows=[]
        for branch in cache.values():
            data=branch[g];rows.append(data['target'])
            rows.extend(r for year in data['seasonal_years'].values() for r in year.values())
        denominator=1
        for row in rows:
            for value in row.values():
                if value is not None:denominator=math.lcm(denominator,value.denominator)
        for row in rows:
            for h,value in row.items():
                if value is not None:
                    encoded=value*denominator
                    assert encoded.denominator==1
                    row[h]=encoded.numerator
        feature['rank_scale']=feature['source_grid_scale']*denominator
        feature['exact_common_denominator']=denominator

def q1_trim_a(result):
    r = copy.deepcopy(result)
    for f in r.get('factors', {}).values():
        f.pop('configurations', None)
    return r

def q1_common_support(features, values, primary):
    ids = [g for g,f in primary['factors'].items() if f.get('score_S') is not None]
    out = copy.deepcopy(values)
    counts = {}
    for year in A_DEFAULT_YEARS:
        ys = str(year); common = []
        for offset in range(-15,16):
            os_ = str(offset)
            if all(all(v is not None for v in out[g]['years'][ys][os_].values()) for g in ids): common.append(os_)
        counts[ys] = len(common)
        for g in out:
            out[g]['years'][ys] = {o:v for o,v in out[g]['years'][ys].items() if o in common}
    return out, {'factor_ids': ids, 'paired_endpoints_by_year': counts, 'note': '共同全部窗口/因素支持；不满足25/4仍NA'}

def q1_rank_rows(scores):
    valid = sorted([(k,v) for k,v in scores.items() if v is not None], key=lambda z:(-z[1],z[0]))
    result = []
    for k,v in valid:
        result.append({'id':k,'rank_min':1+sum(x>v for _,x in valid),'rank_max':sum(x>=v for _,x in valid),'score':v})
    return result


import datetime as dt
from fractions import Fraction
import math

import numpy as np


P2_YEARS = (2021, 2022, 2023, 2024, 2025)
P2_HEIGHTS = (1, 2, 4, 10, 20)
P2_QUANTILES = ((1, 20), (1, 4), (1, 2), (3, 4), (19, 20))
P2_QKEYS = ('q05', 'q25', 'median', 'q75', 'q95')
P2_CHANNELS = (
    ('t2', 'Ta_2m(degree C)', '2 m气温', '°C', 'temperature'),
    ('rh2', 'RH_2m(%)', '2 m相对湿度', '%', 'humidity'),
    ('ws10', 'WS_10m(m/s)', '10 m标量风速', 'm s⁻¹', 'wind'),
    ('pressure', 'Pressure(hPa)', '气压原列', 'hPa', 'pressure'),
    *((f'soilt{str(z).replace(".", "")}', f'soilT_{z}m(degree C)',
       f'{z} m土温', '°C', 'soil_temperature') for z in (0.1, 0.2, 0.4, 0.8, 1.6)),
    *((f'vwc{str(z).replace(".", "")}', f'VWC_{z}m(m3/m3)',
       f'{z} m体积含水量', 'm³ m⁻³', 'soil_water') for z in (0.1, 0.2, 0.4, 0.8, 1.6)),
    *((x.lower(), f'{x}(W/m2)', f'{x}辐射原通道', 'W m⁻²', 'radiation')
      for x in ('Rsd', 'Rsu', 'Rld', 'Rlu')),
    ('h', 'H(W/m2)', 'H有符号原通道', 'W m⁻²', 'flux'),
    ('le', 'LE(W/m2)', 'LE有符号原通道', 'W m⁻²', 'flux'),
    ('ustar', 'USTAR(m/s)', 'USTAR原通道', 'm s⁻¹', 'friction_velocity'),
    ('g01', 'G_0.1m(W/m2)', '0.1 m板热流原通道', 'W m⁻²', 'ground_heat'),
    ('g02', 'G_0.2m(W/m2)', '0.2 m板热流原通道', 'W m⁻²', 'ground_heat'),
)


def _p2_plain(value):

    if isinstance(value, dict):
        return {str(k): _p2_plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [_p2_plain(v) for v in value]
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating, Fraction)):
        number = float(value)
        return number if math.isfinite(number) else None
    return value


def _p2_equal_year_quantile(groups, probabilities=P2_QUANTILES):


    if len(groups) != len(P2_YEARS) or any(not len(g) for g in groups):
        return [None] * len(probabilities)
    denominator = math.lcm(*(len(g) for g in groups))
    rows = sorted((v, denominator // len(g)) for g in groups for v in g)
    total = len(groups) * denominator
    result = []
    for numerator, divisor in probabilities:
        cumulative = 0
        answer = None
        for value, weight in rows:
            cumulative += weight
            if cumulative * divisor >= numerator * total:
                answer = value
                break
        result.append(answer)
    return result


def _p2_background(values, valid, indices, scale):

    quantiles = {key: [] for key in P2_QKEYS}
    exact = {key: [] for key in P2_QKEYS}
    counts = np.zeros((5, indices.shape[-1]), dtype=int)
    for k in range(indices.shape[-1]):
        groups = []
        for yi in range(5):
            ix = indices[yi, :, k]
            inside = (ix >= 0) & (ix < len(values))
            ix = ix[inside]
            ix = ix[valid[ix]]
            group = [int(v) for v in values[ix]]
            groups.append(group)
            counts[yi, k] = len(group)
        qs = _p2_equal_year_quantile(groups)
        for key, q in zip(P2_QKEYS, qs):
            exact[key].append(q)
            quantiles[key].append(float(q / scale) if q is not None else None)
    return {'quantiles': quantiles, 'exact_scaled_quantiles': exact,
            'history_n_by_year': counts.tolist(),
            'history_years_count': (counts > 0).sum(axis=0).tolist()}


def _p2_values(values, valid, indices, scale):
    return [float(values[i] / scale) if 0 <= i < len(values) and valid[i] else None
            for i in indices]


def _p2_summary(values):
    n = len(values)
    pairs = [(k, float(v)) for k, v in enumerate(values) if v is not None]
    result = {'expected_hours': n, 'available_hours': len(pairs),
              'complete_support': len(pairs) == n,
              'available_fraction': len(pairs) / n if n else None,
              'mean': None, 'min': None, 'max': None, 'slope_per_day': None,
              'first_24h_mean': None, 'last_24h_mean': None,
              'first_24h_available_hours': 0, 'last_24h_available_hours': 0,
              'last_minus_first_24h_mean': None}
    if pairs:
        x = np.array([p[0] for p in pairs], dtype=float)
        y = np.array([p[1] for p in pairs], dtype=float)
        result.update(mean=float(y.mean()), min=float(y.min()), max=float(y.max()))
        if len(x) > 1 and np.ptp(x) > 0:
            result['slope_per_day'] = float(np.dot(x-x.mean(), y-y.mean()) /
                                           np.dot(x-x.mean(), x-x.mean()) * 24)
    for part, seq in (('first', values[:24]), ('last', values[-24:])):
        good = [v for v in seq if v is not None]
        result[f'{part}_24h_available_hours'] = len(good)
        result[f'{part}_24h_mean'] = float(np.mean(good)) if good else None
    if result['first_24h_mean'] is not None and result['last_24h_mean'] is not None:
        result['last_minus_first_24h_mean'] = result['last_24h_mean'] - result['first_24h_mean']
    return result


def _p2_qc_note(identity):
    if identity == 'not_provided':
        return '源文件未提供QC列；D0有效原值；未伪造QC0，也不等于已通过通量专项质控。'
    if identity == 'source_not_read':
        return '所需源列未传入；本列为NA，未构造有效值或QC。'
    return '保留调用方D0有效掩码及真实QC身份；未新增插补、截尾、翻号或异常值筛除。'


def _p2_episode_records(values, valid, thresholds, target_indices, relative,
                        support_ends, direction, scale):
    assessed, excess, exceedance = [], [], []
    for k, i in enumerate(target_indices):
        available = 0 <= i < len(values) and bool(valid[i]) and thresholds[k] is not None
        assessed.append(available)
        if not available:
            exceedance.append(None)
            excess.append(None)
        else:
            delta = (int(values[i]) - thresholds[k]) * (1 if direction == 'high' else -1)
            exceedance.append(bool(delta > 0))
            excess.append(float(delta / scale))
    episodes = []
    start = None
    for k in range(len(target_indices) + 1):
        hit = k < len(target_indices) and exceedance[k] is True
        if hit and start is None:
            start = k
        if not hit and start is not None:
            stop = k
            peak = max(range(start, stop), key=lambda j: excess[j])
            episodes.append({
                'start_support_iso': (support_ends[start] - dt.timedelta(hours=1)).isoformat(),
                'end_support_iso': support_ends[stop-1].isoformat(),
                'start_support_end_relative_hours': relative[start],
                'end_support_end_relative_hours': relative[stop-1],
                'duration_hours': stop-start,
                'peak_support_end_iso': support_ends[peak].isoformat(),
                'peak_value': float(values[target_indices[peak]] / scale),
                'peak_threshold': float(thresholds[peak] / scale),
                'peak_excess': excess[peak],
                'support_complete_within_episode': True,
            })
            start = None
    summaries = {}
    for label, hours in (('30d', 720), ('7d', 168), ('72h', 72)):
        hours = min(hours, len(target_indices))
        offset = len(target_indices) - hours
        sub = exceedance[offset:]
        available = sum(0 <= i < len(values) and bool(valid[i]) for i in target_indices[offset:])
        known = sum(v is not None for v in sub)
        hits = sum(v is True for v in sub)
        lengths = []
        window_episodes = []
        run = 0
        for local_k, hit in enumerate(sub + [None]):
            if hit is True:
                run += 1
            elif run:
                lengths.append(run)
                first_k, last_k = offset+local_k-run, offset+local_k-1
                window_episodes.append({
                    'start_support_iso': (support_ends[first_k]-dt.timedelta(hours=1)).isoformat(),
                    'end_support_iso': support_ends[last_k].isoformat(),
                    'duration_hours': run,
                    'clipped_to_display_window': first_k == offset and offset > 0 and exceedance[offset-1] is True,
                })
                run = 0
        good_k = [k for k in range(offset, len(target_indices)) if assessed[k]]
        peak_k = max(good_k, key=lambda k: excess[k]) if good_k else None
        observed_k = [k for k in range(offset, len(target_indices))
                      if 0 <= target_indices[k] < len(values) and valid[target_indices[k]]]
        physical_k = ((max if direction == 'high' else min)(observed_k,
                       key=lambda k: values[target_indices[k]]) if observed_k else None)
        summaries[label] = {
            'expected_hours': hours, 'observed_available_hours': int(available),
            'threshold_available_hours': sum(t is not None for t in thresholds[offset:]),
            'assessed_hours': known, 'unassessed_hours': hours-known,
            'complete_support': known == hours,
            'known_exceedance_hours': hits if known else None,
            'fraction_among_assessed_hours': hits/known if known else None,
            'longest_known_episode_hours': max(lengths, default=0) if known else None,
            'longest_known_episodes': [r for r in window_episodes
                                      if r['duration_hours'] == max(lengths, default=0)],
            'observed_min': float(min(values[target_indices[k]] for k in observed_k)/scale) if observed_k else None,
            'observed_max': float(max(values[target_indices[k]] for k in observed_k)/scale) if observed_k else None,
            'observed_extreme_value': float(values[target_indices[physical_k]]/scale) if physical_k is not None else None,
            'observed_extreme_support_end': support_ends[physical_k].isoformat() if physical_k is not None else None,
            'max_signed_excess': excess[peak_k] if peak_k is not None else None,
            'value_at_max_signed_excess': float(values[target_indices[peak_k]] / scale) if peak_k is not None else None,
            'threshold_at_max_signed_excess': float(thresholds[peak_k] / scale) if peak_k is not None else None,
            'support_end_at_max_signed_excess': support_ends[peak_k].isoformat() if peak_k is not None else None,
        }
    return {'exceedance': exceedance, 'signed_excess': excess, 'episodes': episodes,
            'summaries': summaries,
            'episode_rule': '严格>高阈值或严格<低阈值；每条为假设下完整1h支持；缺观测或缺历史阈值均断段，绝不视为0/非超阈。'}


def build_processes(arrays, endpoint='H-END', cutoff_hour=10):


    if endpoint not in ('H-END', 'H-START'):
        raise ValueError('endpoint must be H-END or H-START')
    if not isinstance(cutoff_hour, int) or not 0 <= cutoff_hour <= 23:
        raise ValueError('cutoff_hour must be an integer between 0 and 23')
    tz = dt.timezone(dt.timedelta(hours=8))
    origin = dt.datetime(2021, 1, 1, tzinfo=tz)
    cutoff = dt.datetime(2026, 8, 26, cutoff_hour, tzinfo=tz)
    end = int((cutoff-origin).total_seconds() // 3600)
    relative = list(range(-719, 1))
    support_ends = [cutoff + dt.timedelta(hours=h) for h in relative]
    shift = 0 if endpoint == 'H-END' else -1
    target_indices = np.array([end+h+shift for h in relative], dtype=int)
    hist_indices = np.array([[[int((cutoff.replace(year=year)-origin).total_seconds() // 3600) +
                                  day*24+h+shift for h in relative]
                                 for day in range(-15, 16)] for year in P2_YEARS], dtype=int)
    grid_length = max((len(v['values']) for v in arrays.values()), default=end+1)
    cache = {}

    def source(column):
        if column in arrays:
            a = arrays[column]
            if len(a['values']) != len(a['valid']):
                raise ValueError(f'values/valid length mismatch: {column}')
            return a['values'], np.asarray(a['valid'], dtype=bool), 10**int(a['p']), a.get('qc_identity', 'not_supplied_to_process_builder')
        return np.zeros(grid_length, dtype=object), np.zeros(grid_length, dtype=bool), 1, 'source_not_read'

    def series(column):
        if column not in cache:
            values, valid, scale, qc = source(column)
            background = _p2_background(values, valid, hist_indices, scale)
            cache[column] = {'values': _p2_values(values, valid, target_indices, scale),
                             'background': background, 'qc_identity': qc}
        return cache[column]

    channels = []
    for key, column, label, unit, group in P2_CHANNELS:
        record = series(column)
        channel = {'id': key, 'column': column, 'label': label, 'unit': unit,
                   'group': group, 'values': record['values'],
                   'qc_identity': record['qc_identity'], 'qc_note': _p2_qc_note(record['qc_identity']),
                   'history_n_by_year': record['background']['history_n_by_year'],
                   'history_years_count': record['background']['history_years_count'],
                   'summaries': {label: _p2_summary(record['values'][-hours:])
                                 for label, hours in (('30d', 720), ('7d', 168), ('72h', 72))}}
        channel.update({f'background_{name}': values for name, values in record['background']['quantiles'].items()})
        if group in ('radiation', 'flux', 'ground_heat'):
            channel['interpretation'] = ('同名通道跨年可比的条件描述；保留有符号原值，不合成净能量、不取绝对值或翻号。'
                                          'G为给定深度板热流，不是订正地表热流；H/LE正向及处理链未由本统计新核实。')
        elif group == 'friction_velocity':
            channel['interpretation'] = 'USTAR同名原通道条件描述；未补低u*筛选，未推Obukhov长度或完整湍流机制。'
        channels.append(channel)

    tcols = [f'Ta_{z}m(degree C)' for z in P2_HEIGHTS]
    ucols = [f'WS_{z}m(m/s)' for z in P2_HEIGHTS]

    def paired(upper, lower):
        high, vh, sh, qh = source(upper)
        low, vl, sl, ql = source(lower)
        scale = max(sh, sl)
        n = min(len(high), len(low))
        values = np.zeros(n, dtype=object)
        valid = vh[:n] & vl[:n]


        valid[min(len(valid), end+shift+1):] = False
        ix = np.flatnonzero(valid)
        values[ix] = (np.asarray(high[ix], dtype=object) * (scale // sh) -
                      np.asarray(low[ix], dtype=object) * (scale // sl))
        return {'values': values, 'valid': valid, 'scale': scale,
                'qc_identity': {'upper': qh, 'lower': ql},
                'target': _p2_values(values, valid, target_indices, scale),
                'background': _p2_background(values, valid, hist_indices, scale)}

    td = paired(tcols[-1], tcols[1])
    ud = paired(ucols[-1], ucols[1])

    def joint_mask(columns):
        records = [source(c) for c in columns]
        n = min(len(r[1]) for r in records)
        return records, np.logical_and.reduce([r[1][:n] for r in records])

    def snapshot(columns, relative_hour):
        records, joint = joint_mask(columns)
        k = relative.index(relative_hour)
        target = [series(c)['values'][k] for c in columns]
        groups_by_height = [[] for _ in columns]
        counts = []
        for yi in range(5):
            ix = hist_indices[yi, :, k]
            ix = ix[(ix >= 0) & (ix < len(joint))]
            ix = ix[joint[ix]]
            counts.append(len(ix))
            for z, (values, _, scale, _) in enumerate(records):
                groups_by_height[z].append([Fraction(int(values[i]), scale) for i in ix])
        q = [_p2_equal_year_quantile(groups) for groups in groups_by_height]
        return {'values': target, 'valid_count': sum(v is not None for v in target),
                'target_complete_profile': all(v is not None for v in target),
                'background_q25': [v[1] for v in q],
                'background_median': [v[2] for v in q],
                'background_q75': [v[3] for v in q],
                'history_n_by_year': counts,
                'historical_joint_mask': 'all five heights valid at the same hour before aggregation',
                'qc_identity': [r[3] for r in records]}

    profiles = [{'support_end_relative_hours': h,
                 'support_end_iso': (cutoff+dt.timedelta(hours=h)).isoformat(),
                 'temperature': snapshot(tcols, h), 'wind': snapshot(ucols, h),
                 'statistic': 'one nominal hourly record per height; not a 24h mean'}
                for h in (-72, -24, 0)]
    time_height = {'support_end_relative_hours': relative,
                   'heights_m': list(P2_HEIGHTS),
                   'temperature_c': [series(c)['values'] for c in tcols],
                   'wind_ms': [series(c)['values'] for c in ucols],
                   'temperature_difference_c': td['target'],
                   'wind_difference_ms': ud['target'],
                   'temperature_difference_summaries': {label: _p2_summary(td['target'][-hours:])
                        for label, hours in (('30d', 720), ('7d', 168), ('72h', 72))},
                   'wind_difference_summaries': {label: _p2_summary(ud['target'][-hours:])
                        for label, hours in (('30d', 720), ('7d', 168), ('72h', 72))},
                   'temperature_pair_qc_identity': td['qc_identity'],
                   'wind_pair_qc_identity': ud['qc_identity'],
                   'temperature_pair_history_n_by_year': td['background']['history_n_by_year'],
                   'wind_pair_history_n_by_year': ud['background']['history_n_by_year']}
    for prefix, pair in (('temperature', td), ('wind', ud)):
        time_height.update({f'pair_{prefix}_background_{name}': value
                            for name, value in pair['background']['quantiles'].items()})


    def period_profile(columns, period_mask):
        records, joint = joint_mask(columns)
        ix = target_indices[period_mask]
        ix = ix[(ix >= 0) & (ix < len(joint))]
        ix = ix[joint[ix]]
        values = [Fraction(sum(int(v[i]) for i in ix), len(ix)*scale) if len(ix) else None
                  for v, _, scale, _ in records]
        groups_by_height = [[] for _ in columns]
        counts = []
        support_counts = []
        for yi in range(5):
            yearly = [[] for _ in columns]
            yearly_support = []
            for oi in range(31):
                hi = hist_indices[yi, oi, period_mask]
                hi = hi[(hi >= 0) & (hi < len(joint))]
                hi = hi[joint[hi]]
                yearly_support.append(len(hi))
                if len(hi):
                    for zi, (v, _, scale, _) in enumerate(records):
                        yearly[zi].append(Fraction(sum(int(v[i]) for i in hi), len(hi)*scale))
            counts.append(len(yearly[0]))
            support_counts.append(yearly_support)
            for zi in range(len(columns)):
                groups_by_height[zi].append(yearly[zi])
        qs = [_p2_equal_year_quantile(g) for g in groups_by_height]
        return {'values': values, 'background_q25': [q[1] for q in qs],
                'background_median': [q[2] for q in qs], 'background_q75': [q[3] for q in qs],
                'valid_counts': [len(ix)]*len(columns),
                'expected_hours': int(period_mask.sum()), 'history_n_by_year': counts,
                'history_available_hours_by_year_trajectory': support_counts,
                'joint_all_heights_before_mean': True}

    def period_pair(pair, period_mask):
        ix = target_indices[period_mask]
        ix = ix[(ix >= 0) & (ix < len(pair['valid']))]
        ix = ix[pair['valid'][ix]]
        groups, support_counts = [], []
        for yi in range(5):
            group, yearly_support = [], []
            for oi in range(31):
                hi = hist_indices[yi, oi, period_mask]
                hi = hi[(hi >= 0) & (hi < len(pair['valid']))]
                hi = hi[pair['valid'][hi]]
                yearly_support.append(len(hi))
                if len(hi):
                    group.append(Fraction(sum(int(pair['values'][i]) for i in hi), len(hi)*pair['scale']))
            groups.append(group)
            support_counts.append(yearly_support)
        qs = _p2_equal_year_quantile(groups)
        return {'value': Fraction(sum(int(pair['values'][i]) for i in ix), len(ix)*pair['scale']) if len(ix) else None,
                'valid_count': len(ix), 'expected_hours': int(period_mask.sum()),
                'background_q25': qs[1], 'background_median': qs[2], 'background_q75': qs[3],
                'history_n_by_year': [len(g) for g in groups],
                'history_available_hours_by_year_trajectory': support_counts,
                'joint_pair_before_mean': True}

    periods = []
    for first in (0, 6, 12, 18):

        mask = np.array([first <= (t-dt.timedelta(hours=1)).hour < first+6 for t in support_ends])
        periods.append({'label': f'{first:02d}–{first+6:02d}', 'hours': list(range(first, first+6)),
                        'temperature': period_profile(tcols, mask),
                        'wind': period_profile(ucols, mask),
                        'temperature_difference': period_pair(td, mask),
                        'wind_difference': period_pair(ud, mask)})

    extreme_series = []
    for key, column, label, direction, quantile in (
        ('t2_low_p05', tcols[1], '2 m气温低于历史P05', 'low', 'q05'),
        ('t2_high_p95', tcols[1], '2 m气温高于历史P95', 'high', 'q95'),
        ('ws10_high_p95', ucols[3], '10 m风速高于历史P95', 'high', 'q95'),
    ):
        values, valid, scale, qc = source(column)
        background = series(column)['background']
        thresholds = background['exact_scaled_quantiles'][quantile]
        record = {'id': key, 'label': label, 'column': column, 'direction': direction,
                  'channel_id': 't2' if column == tcols[1] else 'ws10',
                  'unit': '°C' if column == tcols[1] else 'm s⁻¹',
                  'quantile': quantile, 'values': series(column)['values'],
                  'threshold': background['quantiles'][quantile],
                  'threshold_exact_scaled': [str(t) if t is not None else None for t in thresholds],
                  'scale': str(scale), 'history_n_by_year': background['history_n_by_year'],
                  'history_years_count': background['history_years_count'], 'qc_identity': qc}
        record.update(_p2_episode_records(values, valid, thresholds, target_indices,
                                          relative, support_ends, direction, scale))
        extreme_series.append(record)

    note = ('逐目标日期/钟点，2021–2025每年同历日±15日(31个候选小时)，'
            '按每年有效n_y给每条1/(5 n_y)权重；五年任一年n_y=0则分位NA；'
            '逆经验CDF inf{x:F(x)≥p}，不插值；分位带是历史背景，非置信区间。')
    result = {
        'metadata': {'schema_version': 'q1-process-v2.1', 'primary_cutoff': cutoff.isoformat(),
                     'endpoint': endpoint, 'timezone': 'UTC+8', 'origin': origin.isoformat(),
                     'reference_years': list(P2_YEARS), 'seasonal_days_each_side': 15,
                     'hour_label_assumption': 'H-END: label h supports [h-1h,h); H-START: label h supports [h,h+1h)',
                     'time_coordinate': 'support interval end relative to cutoff in hours',
                     'nominal_support_hours': 1, 'actual_processing_support_confirmed': False,
                     'channel_selection': 'fixed scientific process groups; independent of Top5/rank',
                     'future_support_used': False,
                     'background_quantile_note': note,
                     'no_inference': '不据这些过程认定事件触发、致因、真实日照或完整PBL高度。'},
        'near_surface': {'heights_m': list(P2_HEIGHTS), 'profiles': profiles,
                         'time_height': time_height,
                         'diurnal': {'periods': periods,
                                     'scope_note': ('固定UTC+8 00–06/06–12/12–18/18–24，按支持区间起点小时分组，'
                                                    '不等同真实昼夜。目标为30d内共同有效小时均值；'
                                                    '历史为每年31条同季30d轨迹各自的同固定时段共同有效均值，'
                                                    '再按年等权经验分位；轨迹重叠，不是独立样本/置信区间。')},
                         'scope_note': '1–20 m近地层记录；T20−T2为温差而非位温梯度；WS20−WS2为标量风速差而非矢量切变；跨列先共同有效配对。'},
        'weather': {'time_relative_hours': relative, 'support_end_iso': [t.isoformat() for t in support_ends],
                    'support_start_iso': [(t-dt.timedelta(hours=1)).isoformat() for t in support_ends],
                    'channels': channels, 'seasonal_background_note': note,
                    'display_windows_hours': [720, 168, 72],
                    'summary_note': ('均值/极值/OLS仅描述已有效支持，附可用时长；不补齐缺格。'
                                     '末24h与首24h差/斜率不等于末段单调变化。窗口用于固定过程展示，不新增排名择窗。'),
                    'precipitation': {'status': 'not_computed_unknown_increment_or_accumulation',
                                      'note': '当前Precip增量/累计/复位语义未闭合；未求和/差分/雨强，不能据此声称无雨或无强雨。'}},
        'extremes': {'threshold_rule': note,
                     'predefined_tests': ['T2<P05', 'T2>P95', 'WS10>P95'],
                     'revision_context': 'V2为已看V1结果后的设计修订；本组阈值固定后不按2026命中情况挑阈值，不称事前盲定。',
                     'time_relative_hours': relative, 'series': extreme_series,
                     'support_note': '完整支持/可用支持分别报告；已知超阈小时数仅为可评估小时内计数；缺失/缺阈值断段，不推断缺口内事件。'},
    }
    return _p2_plain(result)


def process_synthetic_selftest():


    groups = [[0], [0], [0], [100], [100]*31]
    assert _p2_equal_year_quantile(groups, ((1, 2),))[0] == 0
    assert _p2_equal_year_quantile([[], [0], [0], [0], [0]]) == [None]*5

    big = 2**60
    assert _p2_equal_year_quantile([[big+i] for i in range(5)], ((1, 2),))[0] == big+2
    tz = dt.timezone(dt.timedelta(hours=8))
    ends = [dt.datetime(2026, 8, 26, h, tzinfo=tz) for h in range(1, 6)]
    vals = np.array([11, 11, 11, 11, 11], dtype=object)
    mask = np.array([True, True, False, True, True])
    rec = _p2_episode_records(vals, mask, [10]*5, np.arange(5), list(range(-4, 1)), ends, 'high', 1)
    assert [e['duration_hours'] for e in rec['episodes']] == [2, 2]
    assert rec['exceedance'][2] is None

    origin = dt.datetime(2021, 1, 1, tzinfo=tz)
    cutoff = dt.datetime(2026, 8, 26, 10, tzinfo=tz)
    end = int((cutoff-origin).total_seconds()//3600)
    n = end+3
    data = {}
    for h in P2_HEIGHTS:
        for prefix, suffix, base in (('Ta', '(degree C)', h), ('WS', '(m/s)', h*2)):
            data[f'{prefix}_{h}m{suffix}'] = {'values': np.full(n, base, dtype=object),
                                           'valid': np.ones(n, dtype=bool), 'p': 0,
                                           'qc_identity': 'synthetic'}

    data['Ta_20m(degree C)']['values'][end] = 200
    data['Ta_2m(degree C)']['valid'][end] = False
    data['Ta_20m(degree C)']['values'][end+1:] = 99999
    result = build_processes(data)
    th = result['near_surface']['time_height']
    assert th['temperature_c'][-1][-1] == 200
    assert th['temperature_difference_c'][-1] is None
    assert th['temperature_difference_c'][-2] == 18
    assert result['near_surface']['diurnal']['periods'][1]['temperature']['values'][-1] == 20
    assert max(v for row in th['temperature_c'] for v in row if v is not None) == 200
    assert result['weather']['support_end_iso'][-1] == cutoff.isoformat()

    reference_before = result['weather']['channels'][0]['background_q95']
    data['Ta_2m(degree C)']['values'][end-719:end+1] = big
    changed = build_processes(data, endpoint='H-START')
    assert changed['weather']['channels'][0]['background_q95'] == reference_before
    assert changed['near_surface']['time_height']['temperature_c'][-1][-1] == 20
    import json
    json.dumps(result, allow_nan=False)
    return {'status': 'pass', 'checks': ['exact_equal_year_inverse_cdf', 'all_five_years_required',
            'integer_threshold_precision', 'missing_hour_breaks_episode', 'pairwise_common_mask',
            'complete_support_before_cutoff', 'H_END_vs_H_START', 'historical_threshold_excludes_2026',
            'strict_json_no_nan'], 'scope': 'synthetic author self-check; not independent scientific acceptance'}


import copy as v2b_copy
import math as v2b_math
import numbers as v2b_numbers
from fractions import Fraction as v2b_Fraction


V2B_NATIVE_V1_IDS = (
    't2_mean', 't2_max', 'pdd', 'warm_rate', 'vwc01', 'vwc02',
    'soilt01', 'p_mean', 'rh2_mean', 'ws10_mean',
)


def v2b_finite(value):
    return (isinstance(value, v2b_numbers.Real) and not isinstance(value, bool)
            and v2b_math.isfinite(value))


def v2b_plain(value):

    if value is None or isinstance(value, (str, bool)):
        return value
    if isinstance(value, dict):
        return {str(key): v2b_plain(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [v2b_plain(item) for item in value]
    if isinstance(value, v2b_numbers.Integral):
        return int(value)
    if isinstance(value, v2b_numbers.Real):
        return float(value) if v2b_math.isfinite(value) else None
    raise TypeError('comparison output is not strict JSON: ' + type(value).__name__)


def v2b_specs(features):
    if isinstance(features, dict):
        features = features.get('factors', features.get('features', []))
    out = {str(feature['id']): feature for feature in features}
    if len(out) != len(features):
        raise ValueError('duplicate V2 comparison feature id')
    return out


def v2b_rows_by_id(method):
    return {str(row.get('factor_id', row.get('id'))): row
            for row in method.get('factor_rows', [])}


def v2b_a_scores(result, ids):
    scores = {}
    for factor_id in ids:
        factor = result.get('factors', {}).get(factor_id, {})
        if factor.get('rank_eligible') and v2b_finite(factor.get('score_S')):
            exact = factor.get('score_S_exact')
            scores[factor_id] = (v2b_Fraction(exact) if exact is not None
                                 else factor['score_S'])
    return scores


def v2b_ranks(scores):
    return {row['id']: row for row in q1_rank_rows(scores)}


def v2b_interval(rows, factor_id):
    row = rows.get(factor_id)
    if row is None or not v2b_finite(row.get('rank_min')):
        return None
    return {'rank_min': int(row['rank_min']), 'rank_max': int(row['rank_max'])}


def v2b_reuse_check(baseline, validation, seed):

    if not isinstance(validation, dict) or validation.get('status') != 'pass':
        return False, 'no passed caller-side input/code/configuration validation'
    digest = validation.get('file_sha')
    if (not isinstance(digest, str) or len(digest) != 64
            or any(char not in '0123456789abcdefABCDEF' for char in digest)):
        return False, 'validated baseline file SHA256 is absent or malformed'
    controls = validation.get('verified_controls')
    if not isinstance(controls, (dict, list, tuple)) or not controls:
        return False, 'verified_controls must contain caller-side evidence, not a boolean'
    if not isinstance(baseline, dict):
        return False, 'baseline analysis result is absent'
    method = baseline.get('methods', {}).get('B_aligned')
    if not isinstance(method, dict):
        return False, 'baseline methods.B_aligned is absent'
    rows = method.get('factor_rows', [])
    ids = [row.get('factor_id') for row in rows]
    if len(ids) != len(V2B_NATIVE_V1_IDS) or set(ids) != set(V2B_NATIVE_V1_IDS):
        return False, 'baseline B is not the fixed native V1 ten-factor family'
    metadata = method.get('metadata', {})
    expected = {'seed': seed, 'n_perm_requested': 2000,
                'n_boot_year_resamples_executed': 1000, 'n_perm_boot_requested': 500}
    if any(metadata.get(key) != value for key, value in expected.items()):
        return False, 'baseline seed or formal draw configuration differs'
    if metadata.get('formal_counts') is not True:
        return False, 'baseline does not record the required formal counts'
    return True, 'reuse authorized by caller-side evidence bound to baseline SHA256'


def v2b_fit(specs, values, ids, seed, label, role):

    if len(ids) < 2:
        return {'factor_rows': [
                    {'factor_id': factor_id, 'status': 'insufficient_comparison_family',
                     'topsis_score': None, 'rank': None, 'rank_min': None,
                     'rank_max': None, 'perm_draws_executed': 0}
                    for factor_id in ids],
                'window_rows': [], 'bootstrap_rows': [], 'comparisons': {},
                'metadata': {
                    'label': label, 'family_role': role, 'input_factor_ids': list(ids),
                    'n_input_original_factors': len(ids),
                    'n_computable_topsis_factors': 0, 'seed': seed,
                    'execution_in_v2': False, 'execution_class': 'not_executed_insufficient_family',
                    'reason': 'TOPSIS requires at least two rows; no support threshold relaxed',
                    'n_perm_requested': 2000, 'n_boot_year_resamples_requested': 1000,
                    'n_perm_boot_requested': 500, 'n_boot_year_resamples_executed': 0,
                    'new_main_pseudo_draws_executed': 0,
                    'new_bootstrap_pseudo_draws_executed': 0,
                    'entropy_weights_A_neglogp_mechanism': [None, None, None],
                    'topsis_status': 'fewer_than_two_input_rows'}}
    result = qf_run_ref([specs[factor_id] for factor_id in ids], values, seed=seed,
                       label=label, n_perm=2000, n_boot=1000, n_perm_boot=500)
    metadata = result.setdefault('metadata', {})
    metadata.update({
        'family_role': role, 'input_factor_ids': list(ids), 'execution_in_v2': True,
        'reuse_status': 'new_actual_fit',
        'new_main_pseudo_draws_executed': sum(
            row.get('perm_draws_executed', 0) for row in result.get('factor_rows', [])),
        'new_bootstrap_pseudo_draws_executed': sum(
            row.get('pseudo_draws_executed', 0) for row in result.get('bootstrap_rows', [])),
        'fitting_scope': 'entropy weights, TOPSIS, BH and each bootstrap fit use exactly this input family; unavailable rows remain NA',
        'proxy_policy': 'only original24 factors with existing original mechanism scores; no proxy mechanism score invented',
        'D0_alignment': 'caller-supplied D0 cache and applicable W; no filling, clipping or recleaning',
        'original_pdd_unit_alignment': 'dictionary degree-C h; constant /24 to degree-C d leaves ranks and B unchanged',
    })
    if result.pop('_final_reused_modex',False):
        metadata.update(execution_in_final=False,reuse_status='frozen_参考方法_reused_after_exact_input_check',
                        new_main_pseudo_draws_executed=0,new_bootstrap_pseudo_draws_executed=0,
                        n_boot_year_resamples_executed_in_final=0,
                        count_fields_scope='source-run counts retained; no new 参考方法 draws in this run')
    else:
        metadata['execution_in_final']=True
    for row in result.get('window_rows', []):
        factor_id = row['factor_id']
        scale = specs[factor_id].get('rank_scale')
        if row.get('event_value') is not None and isinstance(scale, int) and scale > 0:
            row['event_value_physical'] = row['event_value'] / scale
    return result


def v2b_common_table(ids, specs, e1, e2, b_result, c_index, *, scope):
    e1_ranks = v2b_ranks(v2b_a_scores(e1, ids))
    e2_ranks = v2b_ranks(v2b_a_scores(e2, ids))
    b_rows = v2b_rows_by_id(b_result)
    b_ranks = v2b_ranks({factor_id: b_rows[factor_id]['topsis_score']
                        for factor_id in ids if factor_id in b_rows
                        and v2b_finite(b_rows[factor_id].get('topsis_score'))})
    c_ranks = v2b_ranks({factor_id: c_index[factor_id]['A'] for factor_id in ids
                        if factor_id in c_index and v2b_finite(c_index[factor_id].get('A'))})
    e1_native = {row['id']: row for row in e1.get('ranking_all', [])}
    e2_native = {row['id']: row for row in e2.get('ranking_all', [])}
    e1_dedup = {row['id']: row for row in e1.get('candidate_ranking', [])}
    e2_dedup = {row['id']: row for row in e2.get('candidate_ranking', [])}
    return [{
        'id': factor_id, 'label': specs[factor_id].get('label', factor_id),
        'E1_rank': v2b_interval(e1_ranks, factor_id),
        'E2_rank': v2b_interval(e2_ranks, factor_id),
        'B_rank': v2b_interval(b_ranks, factor_id),
        'C_rank': v2b_interval(c_ranks, factor_id),
        'E1_score_S': e1.get('factors', {}).get(factor_id, {}).get('score_S'),
        'E2_score_S': e2['factors'][factor_id].get('score_S'),
        'B_topsis_score': b_rows.get(factor_id, {}).get('topsis_score'),
        'C_fixed7d_A': c_index[factor_id].get('A'),
        'E1_native_rank': v2b_interval(e1_native, factor_id),
        'E2_native_rank': v2b_interval(e2_native, factor_id),
        'E1_native_deduplicated_rank': v2b_interval(e1_dedup, factor_id),
        'E2_native_deduplicated_rank': v2b_interval(e2_dedup, factor_id),
        'C_native_rank': v2b_interval(c_index, factor_id),
        'scope': scope,
        'E1_scope': 'E1 original10 intersect displayed IDs only; newly restored factors are NA',
    } for factor_id in ids]


def build_v2_comparisons(features, values, e1, e2, c_rows, baseline=None,
                         baseline_validation=None, seed=20261024):


    specs = v2b_specs(features)
    absent = set(V2B_NATIVE_V1_IDS) - set(specs)
    if absent:
        raise ValueError('V2 dictionary lacks native V1 factors: ' + ', '.join(sorted(absent)))
    original_ids = [factor[0] for factor in ORIGINAL_FACTORS]
    expanded_ids = [factor_id for factor_id in original_ids if factor_id in specs
                    and str(specs[factor_id].get('status', '')).upper() != 'NA'
                    and specs[factor_id].get('window_hours')]
    reuse, reuse_reason = v2b_reuse_check(baseline, baseline_validation, seed)
    if reuse:
        native = v2b_copy.deepcopy(baseline['methods']['B_aligned'])
        native.setdefault('metadata', {}).update({
            'family_role': 'native_V1_original10_preserved',
            'input_factor_ids': list(V2B_NATIVE_V1_IDS),
            'execution_in_v2': False, 'reuse_status': 'reused_verified_prior_run',
            'reuse_source': {'field': 'methods.B_aligned',
                             'baseline_file_sha256': baseline_validation['file_sha'],
                             'caller_validation': v2b_copy.deepcopy(baseline_validation)},
            'new_main_pseudo_draws_executed': 0,
            'new_bootstrap_pseudo_draws_executed': 0,
            'count_interpretation': 'original executed counts are preserved as source-run facts; no new V2 draws for this native result',
        })
    else:
        native = v2b_fit(specs, values, list(V2B_NATIVE_V1_IDS), seed,
                         'B-native-original10-V2-actual-recalculation',
                         'native_V1_original10_recalculated')
        native['metadata']['reuse_not_taken_reason'] = reuse_reason
    expanded = v2b_fit(specs, values, expanded_ids, seed,
                       'B-expanded-original24-V2-actual-fit',
                       'expanded_semantically_available_original24')
    b_expanded = v2b_rows_by_id(expanded)
    a_ids = set(v2b_a_scores(e2, original_ids))
    b_ids = {factor_id for factor_id, row in b_expanded.items()
             if v2b_finite(row.get('topsis_score'))}
    c_index = {str(row['id']): row for row in c_rows}
    c_ids = {factor_id for factor_id, row in c_index.items()
             if factor_id in original_ids and row.get('rank_eligible',True) and not row.get('uninformative',False) and not row.get('reason')
             and v2b_finite(row.get('q')) and v2b_finite(row.get('A'))}
    common_ids = [factor_id for factor_id in original_ids
                  if factor_id in a_ids & b_ids & c_ids]
    common = v2b_fit(specs, values, common_ids, seed,
                     'B-common-original24-V2-actual-refit',
                     'common_original24_before_uniform_redundancy')
    common_set = set(common_ids)
    applied, not_applied = [], []
    for row in e2.get('redundancy', {}).get('suppressed', []):
        factor_id, representative = row['id'], row['representative_id']
        if factor_id not in common_set:
            continue
        if representative in common_set:
            applied.append(v2b_copy.deepcopy(row))
        else:
            not_applied.append({**v2b_copy.deepcopy(row),
                                'comparison_action': 'retained: E2 representative absent from this common pool'})
    removed = {row['id'] for row in applied}
    dedup_ids = [factor_id for factor_id in common_ids if factor_id not in removed]
    dedup = v2b_fit(specs, values, dedup_ids, seed,
                    'B-common-original24-uniform-deduplicated-V2-actual-refit',
                    'common_original24_after_uniform_E2_redundancy')
    methods = {'B_native_v1': native, 'B_expanded_original24': expanded,
               'B_common_refit': common, 'B_common_deduplicated_refit': dedup,
               'C_factors': v2b_copy.deepcopy(c_rows)}
    if reuse and isinstance(baseline.get('methods', {}).get('B_original_cutoff'), dict):
        original_cutoff = v2b_copy.deepcopy(baseline['methods']['B_original_cutoff'])
        original_cutoff.setdefault('metadata', {}).update({
            'execution_in_v2': False, 'reuse_status': 'prior_V1_C09_reference_only',
            'baseline_file_sha256': baseline_validation['file_sha'],
            'new_main_pseudo_draws_executed': 0,
            'new_bootstrap_pseudo_draws_executed': 0,
            'V2_scope': 'preserved source reference; excluded from every V2 common fit and comparison',
        })
        methods['B_native_v1_original_cutoff_reference'] = original_cutoff
    original10_scores1 = v2b_a_scores(e1, V2B_NATIVE_V1_IDS)
    original10_scores2 = v2b_a_scores(e2, V2B_NATIVE_V1_IDS)
    paired10 = [factor_id for factor_id in V2B_NATIVE_V1_IDS
                if factor_id in original10_scores1 and factor_id in original10_scores2]
    ranks1 = v2b_ranks({factor_id: original10_scores1[factor_id] for factor_id in paired10})
    ranks2 = v2b_ranks({factor_id: original10_scores2[factor_id] for factor_id in paired10})
    result = {
        'methods': methods,
        'comparisons_v2': v2b_common_table(
            dedup_ids, specs, e1, e2, dedup, c_index,
            scope='uniform E2 redundancy on common original24 IDs; B actually refitted; A/C scores retained'),
        'comparisons_common_before_redundancy': v2b_common_table(
            common_ids, specs, e1, e2, common, c_index,
            scope='common original24 IDs before redundancy; B actually refitted; A/C scores retained'),
        'comparison_factor_ids': dedup_ids,
        'common_before_redundancy_ids': common_ids,
        'E1_E2_original10_common': [
            {'id': factor_id, 'label': specs[factor_id].get('label', factor_id),
             'E1_rank': v2b_interval(ranks1, factor_id),
             'E2_rank': v2b_interval(ranks2, factor_id),
             'E1_score_S': float(original10_scores1[factor_id]),
             'E2_score_S': float(original10_scores2[factor_id]),
             'scope': 'same rankable original10 IDs; unchanged exact A scores, ranks reindexed without redundancy'}
            for factor_id in paired10],
        'native_scope': {
            'E1': {'factor_ids': list(e1.get('factors', {})), 'settings': e1.get('settings'),
                   'rankable_ids': list(v2b_a_scores(e1, specs)),
                   'deduplicated_ids': [row['id'] for row in e1.get('candidate_ranking', [])]},
            'E2': {'factor_ids': list(e2.get('factors', {})), 'settings': e2.get('settings'),
                   'rankable_ids': list(v2b_a_scores(e2, specs)),
                   'deduplicated_ids': [row['id'] for row in e2.get('candidate_ranking', [])]},
            'B_native': {'factor_ids': list(V2B_NATIVE_V1_IDS), 'reused': reuse,
                         'reuse_decision_reason': reuse_reason},
            'B_expanded': {'factor_ids': expanded_ids,
                           'original24_semantic_NA_ids': [factor_id for factor_id in original_ids if factor_id not in expanded_ids],
                           'excluded_new_proxy_ids': [factor_id for factor_id in specs if factor_id not in original_ids]},
            'C': {'factor_ids': list(c_index), 'common_legal_original_ids': [factor_id for factor_id in original_ids if factor_id in c_ids],
                  'rule': 'caller-computed fixed 168h D0 C10; annual>=15 endpoints and >=4 years; source gap rule retained'},
            'pool_differences': 'A: factor-specific >=25 joint full-W endpoints/year, >=4 years, annual equal-weight 3+1 calibration; B: pooled >=100 same-date references, Jun15-Sep30 pseudo pool, own full-W support, event-fitted entropy(A,-log10p,original mechanism); C: fixed7d with own valid support. Common IDs do not equal identical annual/pseudo observations or algorithms.',
        },
        'refit_scope': {
            'common_rule': 'rankable non-all-tied E2 intersect computable expanded-B TOPSIS intersect legally computable fixed7d-C, restricted to original24',
            'common_before_redundancy_ids': common_ids,
            'uniform_deduplicated_ids': dedup_ids,
            'redundancy_applied': applied, 'redundancy_not_applied_missing_representative': not_applied,
            'redundancy_source': 'E2 existing suppression: predeclared duplicated physical information AND sufficiently supported historical representative-window abs(Spearman rho)>=0.90; apply only when representative belongs to the same common pool',
            'B_fit_roles': ['B_expanded_original24', 'B_common_refit', 'B_common_deduplicated_refit'],
            'counts_rule': 'each executed new fit requests 2000 main pseudo, 1000 year-block resamples and 500 pseudo per computable bootstrap distribution; actual draws and NA remain in its core rows; <2 input factors is not executed',
            'shared_input_scope': 'same per-factor definitions, full applicable W, C10/H-END D0 cache; observation/support differences remain explicit',
            'A_C_scope': 'A and C own scores retained; exact A fractions used to reindex ranks on listed IDs; no combined score',
            'E1_scope': 'original10 only; E1 newly restored-factor ranks remain NA and E1 available-subset ranks are not a full common-family comparison',
        },
        'comparison_note': 'Native methods preserve their own pools/weights/support. B common-family entropy, TOPSIS and BH are real new fits; A/C rankings restrict unchanged own scores. New proxies remain in A/C native outputs, outside B and common-original24 tables. Rank agreement/overlap is not accuracy or proof of V2 superiority.',
    }
    return v2b_plain(result)


import argparse
import copy
import csv
import datetime as dt
import decimal
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys
import time

Q1_TZ=dt.timezone(dt.timedelta(hours=8))
Q1_ORIGIN=dt.datetime(2021,1,1,tzinfo=Q1_TZ)
Q1_GRID_LENGTH=52584
Q2_E1_IDS=('t2_mean','t2_max','pdd','warm_rate','vwc01','vwc02','soilt01','p_mean','rh2_mean','ws10_mean')
Q2_PROXY_IDS=('t20_t2_diff','ws20_ws2_diff')
Q1_COLUMNS={f['id']:(f['raw_columns'][0]['table'],f['raw_columns'][0]['column'],f['raw_columns'][0].get('qc_column'),0)
            for f in Q1_FROZEN_DICTIONARY['factors'] if f.get('computable')}
Q2_OUTPUTS=('configuration.json','feature_dictionary.json','analysis_results.json','process_v2.json',
            'window_cache.json.gz','window_quality.csv.gz','all_factors.csv','primary_windows.csv',
            'single_sided_hypotheses.csv','experiments.csv','sensitivity.csv','loo_frequency.csv',
            'method_comparison.csv','c_fixed7d.csv','extremes.csv','process_support.csv')


def q1_features():
    d=copy.deepcopy(Q1_FROZEN_DICTIONARY)
    for f in d['factors']:
        f['windows_h']=f['window_hours']
        if f.get('computable'):
            assert f['direction']=='two_sided' and f['eligible_for_top5']
            f['source_grid_scale']=1;f['rank_scale']=1
            f['rank_encoding']='exact_source_grid_Fraction_then_lossless_integer_LCM'
    d['design_status']='看过V1目标结果后的研究修订；不称盲定或独立终评'
    return d


def q2_requested_columns(features):
    wanted={c['column'] for f in features if f.get('computable') for c in f['raw_columns']}
    wanted.update(p[1] for p in P2_CHANNELS)
    wanted.update(f'Ta_{z}m(degree C)' for z in P2_HEIGHTS)
    wanted.update(f'WS_{z}m(m/s)' for z in P2_HEIGHTS)
    return wanted


def q1_read_d0(data_dir,features):
    base=Path(data_dir);summary=json.loads((base/'run_summary_all.json').read_text())
    with (base/'attachment01_02/variable_dictionary.csv').open(encoding='utf-8-sig') as stream:
        metadata=list(csv.DictReader(stream))
    wanted=q2_requested_columns(features);arrays={};counts={}
    lookup={r['variable']:r for r in metadata if r['dataset'] in ('MET','SOIL','RADM','FLUX') and r['variable'] in wanted}
    assert set(lookup)==wanted,'所需D0变量/元数据不完整'
    for table in ('MET','SOIL','RADM','FLUX'):
        file=base/'attachment01_02/station'/(table+'.csv.gz')
        expected=summary['output_sha256']['attachment01_02/station/'+table+'.csv.gz']
        if isinstance(expected,dict):expected=expected['sha256']
        assert q1_sha(file)==expected,'D0已有清理版本不匹配：'+table
        specs={col:r for col,r in lookup.items() if r['dataset']==table}
        for col,r in specs.items():
            no_qc=r['qc_column'] in ('not_provided','',None)
            arrays[col]={'values':np.zeros(Q1_GRID_LENGTH,dtype=object),'valid':np.zeros(Q1_GRID_LENGTH,dtype=bool),'p':0,
                         'unit':r['unit'],'table':table,'qc_column':None if no_qc else r['qc_column'],
                         'qc_identity':'not_provided' if no_qc else 'provided',
                         'qc_note':'源QC列未提供，不能推定QC9或QC0；沿用D0保留值' if no_qc else 'D0原QC0有效掩码，不恢复排除值'}
        seen=set()
        with gzip.open(file,'rt',encoding='utf-8-sig',newline='') as stream:
            for row in csv.DictReader(stream):
                i=q1_hour(row['timestamp_label']);assert i not in seen and 0<=i<Q1_GRID_LENGTH
                seen.add(i)
                for col,r in specs.items():
                    token=row[col]
                    if token=='':continue
                    source=arrays[col];qc=source['qc_column']
                    if qc is not None:assert decimal.Decimal(row[qc])==0,'D0源QC政策不匹配'
                    number=decimal.Decimal(token);assert number.is_finite()
                    places=max(source['p'],-number.as_tuple().exponent)
                    if places>source['p']:
                        source['values']*=10**(places-source['p']);source['p']=places
                    exact=Fraction(number)*10**places;assert exact.denominator==1
                    source['values'][i]=exact.numerator;source['valid'][i]=True
        counts[table]={'labels':len(seen),'channels':{col:{'valid':int(arrays[col]['valid'].sum()),'qc_identity':arrays[col]['qc_identity'],'qc_column':arrays[col]['qc_column'],'decimal_places':arrays[col]['p']} for col in specs}}
    return arrays,counts


def q2_inputs(arrays,feature):
    g=feature['id']
    cached=arrays.setdefault('_q2_feature_cache',{})
    if g in cached:return cached[g]
    names=['Ta_20m(degree C)','Ta_2m(degree C)'] if g=='t20_t2_diff' else ['WS_20m(m/s)','WS_2m(m/s)'] if g=='ws20_ws2_diff' else [feature['raw_columns'][0]['column']]
    p=max(arrays[n]['p'] for n in names)
    raw=arrays[names[0]]['values']*10**(p-arrays[names[0]]['p'])
    mask=arrays[names[0]]['valid'].copy()
    if len(names)==2:
        raw=raw-arrays[names[1]]['values']*10**(p-arrays[names[1]]['p'])
        mask &= arrays[names[1]]['valid']
    cached[g]=(raw,mask,p,names)
    return cached[g]


def q1_window(arrays,feature,end,hours,endpoint,coverage):
    raw_all,valid_all,p,names=q2_inputs(arrays,feature)
    shift=int(endpoint=='H-END');lo,hi=end-hours+shift,end+shift
    meta={'expected':hours,'valid':0,'coverage':0.,'reason':None,'joint_columns':names}
    if lo<0 or hi>Q1_GRID_LENGTH:meta['reason']='outside_D0';return None,meta
    raw,mask=raw_all[lo:hi],valid_all[lo:hi];n=int(mask.sum());meta.update(valid=n,coverage=n/hours)
    if n<math.ceil(coverage*hours):meta['reason']='hour_coverage_below_threshold';return None,meta
    g=feature['id']
    if g in ('pdd','t2_max'):
        if n!=hours:meta['reason']='complete_interval_required';return None,meta
        value=Fraction(int(np.maximum(raw,0).sum())) if g=='pdd' else Fraction(int(raw.max()))
    elif g=='warm_rate':
        days=hours//24;xs=[];ys=[]
        for d in range(days):
            m=mask[d*24:(d+1)*24];r=raw[d*24:(d+1)*24]
            if int(m.sum())>=math.ceil(coverage*24):xs.append(Fraction(2*d+1,2));ys.append(Fraction(int(r[m].sum()),int(m.sum())))
        meta['daily_blocks_valid']=len(xs)
        if len(xs)<max(3,math.ceil(coverage*days)) or not xs or xs[0]!=Fraction(1,2) or xs[-1]!=Fraction(2*days-1,2):
            meta['reason']='slope_block_support_or_endpoints';return None,meta
        xm,ym=sum(xs)/len(xs),sum(ys)/len(ys)
        value=sum((x-xm)*(y-ym) for x,y in zip(xs,ys))/sum((x-xm)**2 for x in xs)
    else:value=Fraction(int(raw[mask].sum()),n)
    assert hi-shift==end and lo-shift==end-hours
    return value,meta


def q2_c_route(arrays,features,values):
    rows=[];gap_exclusions=[]
    for f in features:
        g=f['id'];data=values.get(g,{})
        def value(year,offset,encoded):
            if encoded is None:return None
            if g in ('pdd','t2_max'):return encoded
            _,mask,_,_=q2_inputs(arrays,f);end=q1_hour(f'{year}-08-26T10:00:00+08:00')+24*offset
            run=longest=0
            for valid in mask[end-167:end+1]:run=0 if valid else run+1;longest=max(longest,run)
            if longest>3:gap_exclusions.append({'id':g,'year':year,'offset':offset,'gap':longest});return None
            return encoded
        event=value(2026,0,data.get('target',{}).get('168'));annual={};all_values=[]
        for year,records in data.get('years',{}).items():
            xs=[value(int(year),int(o),r.get('168')) for o,r in records.items()];xs=[x for x in xs if x is not None]
            if len(xs)>=15 and event is not None:annual[year]={'n':len(xs),'cdf':float(a_midrank_fraction(sorted(xs),event))};all_values.extend(xs)
        q=sum(x['cdf'] for x in annual.values())/len(annual) if len(annual)>=4 else None
        uninformative=bool(all_values) and len(set(all_values+[event]))==1
        rows.append({'id':g,'label':f['label'],'unit':f['unit'],'value':event/f['rank_scale'] if event is not None else None,'q':q,'A':abs(2*q-1) if q is not None else None,'annual':annual,'rank_eligible':not uninformative,'uninformative':uninformative,'reason':None if q is not None else 'definition_target_or_15window_4year_support_NA'})
    ranks={r['id']:r for r in q1_rank_rows({r['id']:r['A'] for r in rows if r['rank_eligible']})}
    for row in rows:row.update({k:v for k,v in ranks.get(row['id'],{}).items() if k in ('rank_min','rank_max')})
    return rows,{'mean_gap_rule_hours':3,'gap_exclusions':gap_exclusions,'cutoff_hour':10,'endpoint':'H-END'}


def q2_evidence(features,result,arrays):
    rows=[];index={f['id']:f for f in features};end=q1_hour('2026-08-26T10:00:00+08:00')
    for g,factor in result['factors'].items():
        f=index[g];windows=[];directions=set();winning=[];all_counts={'high':0,'low':0,'center':0};mean_directions=set()
        for h,row in factor['windows'].items():
            q=row.get('q_mean');v=row['target_value']
            counts={'high':0,'low':0,'center':0}
            for configuration in factor['configurations']:
                cq=configuration['target']['q_by_window'][str(h)]
                cd='high' if cq>.5 else 'low' if cq<.5 else 'center'
                counts[cd]+=1;all_counts[cd]+=1
            if q is not None:mean_directions.add('high' if q>.5 else 'low' if q<.5 else 'center')
            windows.append({'hours':int(h),'config_direction_counts':counts,'configuration_direction_inconsistent':sum(n>0 for n in counts.values())>1,'physical_value':v/f['rank_scale'] if v is not None else None,'unit':f['unit'],'q_mean':q,'q_min':row.get('q_min'),'q_max':row.get('q_max'),'a_mean':row.get('a_mean'),'position_direction':'high' if q is not None and q>.5 else 'low' if q is not None and q<.5 else 'center' if q is not None else None,'raw_sign':'positive' if v is not None and v>0 else 'negative' if v is not None and v<0 else 'zero' if v==0 else None})
        for config in factor['configurations']:
            for h in config['target']['winning_windows']:
                q=config['target']['q_by_window'][str(h)];direction='high' if q>.5 else 'low' if q<.5 else 'center';directions.add(direction)
                v=factor['windows'][str(h)]['target_value']
                winning.append({'calibration_year':config['calibration_year'],'reference_years':config['reference_years'],'hours':h,'q':q,'direction':direction,'physical_value':v/f['rank_scale']})

        first,m1=q1_window(arrays,f,end-696,24,'H-END',.9) if g!='warm_rate' else (None,{})
        last,m2=q1_window(arrays,f,end,24,'H-END',.9) if g!='warm_rate' else (None,{})
        rows.append({'id':g,'label':f['label'],'unit':f['unit'],'S':factor['score_S'],'windows':windows,'winning_configurations':winning,'winning_directions':sorted(directions),'direction_inconsistent_across_winning_windows_or_configs':len(directions)>1,'all_window_mean_directions':sorted(mean_directions),'all_config_window_direction_counts':all_counts,'direction_inconsistent_across_all_windows_or_configs':sum(n>0 for n in all_counts.values())>1,'first24_of30d':float(first/f['source_grid_scale']) if first is not None else None,'last24':float(last/f['source_grid_scale']) if last is not None else None,'last_minus_first24':float((last-first)/f['source_grid_scale']) if first is not None and last is not None else None,'first_last_support':[m1,m2]})
    return rows


def q2_original_specs(features):
    specs=[]
    for f in features:
        if f['id'] not in Q2_E1_IDS:continue
        s=copy.deepcopy(f);s['direction']=f['hypothesis_direction'];s['eligible_for_top5']=f['id'] not in ('p_mean','rh2_mean','ws10_mean')
        specs.append(s)
    return specs


def q2_hypotheses(features,old_policy):
    rows=[]
    for f in features:
        if f['id'] not in Q2_E1_IDS or f.get('hypothesis_direction')!='high':continue
        r=old_policy['factors'][f['id']];s=r['score_S']
        rows.append({'id':f['id'],'hypothesis':'原V1偏高假说，不反转','direction':'high','S_one_sided':s,'description':'相对突出线索' if s is not None and s>=.8 else '较弱线索' if s is not None and s>=.5 else '未见偏高支持' if s is not None else '支持不足NA','rule':'描述性分级：S>=.8/.5；不是p检验，不参与V2主榜','windows':r['windows']})
    return rows


def q2_validate_v1(args,cache,features,data,before):
    if args.v1_results is None:return None,None,[]
    directory=Path(args.v1_results);files=[directory/'analysis_results.json',directory/'window_cache.json.gz',directory/'RUN_RECEIPT.json']
    if not all(p.is_file() for p in files):raise ValueError('指定V1结果目录缺少已声明复用文件')
    baseline=json.loads(files[0].read_text());receipt=json.loads(files[2].read_text())
    assert receipt['inputs_unchanged'] and baseline['schema']=='Q1_FORMAL_V1'
    assert json.loads((directory/'configuration.json').read_text())['seed']==args.seed

    required=['run_summary_all.json','attachment01_02/variable_dictionary.csv',
              'attachment01_02/station/MET.csv.gz','attachment01_02/station/SOIL.csv.gz']
    verified_inputs=[]
    for relative in required:
        matches=[x for x in receipt['inputs_before'] if x['path'].endswith('/'+relative)]
        assert len(matches)==1 and matches[0]['sha256']==q1_sha(data/relative), 'V1 D0 source fingerprint differs'
        verified_inputs.append({'input':relative,'sha256':matches[0]['sha256']})
    for name in ('analysis_results.json','window_cache.json.gz','configuration.json'):
        matches=[x for x in receipt['outputs'] if x['path'].endswith('/'+name)]
        assert len(matches)==1 and matches[0]['sha256']==q1_sha(directory/name), 'V1 reusable result fingerprint differs'
    with gzip.open(files[1],'rt',encoding='utf-8') as stream:old=json.load(stream)
    old_scale={f['id']:f['rank_scale'] for f in baseline['features'] if f['id'] in Q2_E1_IDS};new_scale={f['id']:f['rank_scale'] for f in features if f['id'] in Q2_E1_IDS}
    n=0
    for case,records in cache.items():
        ep,ch,cov=case;old_case=old[f'{ep}_C{ch}_{cov}']
        for g in Q2_E1_IDS:
            current=records[g];previous=old_case[g]
            rows=[(current['target'],previous['target'])]
            rows.extend((current['seasonal_years'][y][md],previous['seasonal_years'][y][md]) for y in current['seasonal_years'] for md in current['seasonal_years'][y])
            for a,b in rows:
                assert set(a)==set(b)
                for h,value in a.items():
                    other=b[h];assert (value is None)==(other is None)
                    if value is not None:assert value*old_scale[g]==other*new_scale[g]
                    n+=1
    bindings=[q1_fp(p) for p in files+[directory/'configuration.json']]
    return baseline,{'status':'pass','file_sha':q1_sha(files[0]),'verified_controls':['D0 exact input hashes below','same seed, variables, formulas, endpoints/cutoffs/coverage','all original10 target and seasonal physical Fraction values identical'],'physical_window_comparisons':n,'verified_D0_inputs':verified_inputs,'source_bindings':bindings},files+[directory/'configuration.json']


def q2_write_tables(res,result,process):
    specs=result['features'];primary=result['primary'];rank={r['id']:r for r in primary['ranking_all']};merged={x['id']:x['representative_id'] for x in primary['redundancy']['suppressed']};rows=[];win=[]
    for f in specs:
        g=f['id'];r=primary['factors'].get(g,{});rr=rank.get(g,{})
        rows.append({'id':g,'label':f['label'],'unit':f['unit'],'definition_status':f['status'],'statistic_status':r.get('status','definition_NA'),'rank_eligible':r.get('rank_eligible',False),'merged_into':merged.get(g),'direction':f['direction'],'S':r.get('score_S'),'rank_min':rr.get('rank_min'),'rank_max':rr.get('rank_max'),'top_candidate':g in primary['top_candidates'],'boundary':g in primary['boundary_candidates'],'years':r.get('support',{}).get('valid_years'),'annual_counts':r.get('support',{}).get('n_valid_endpoints_by_year'),'reason':r.get('reason') or (f.get('status_reason') if not f.get('computable') else None)})
        for h,w in r.get('windows',{}).items():win.append({'id':g,'hours':h,'value':w['target_value']/f['rank_scale'] if w['target_value'] is not None else None,'unit':f['unit'],**{k:v for k,v in w.items() if k!='target_value'}})
    q1_csv(res/'all_factors.csv',rows);q1_csv(res/'primary_windows.csv',win)
    q1_csv(res/'single_sided_hypotheses.csv',result['single_sided_hypotheses'])
    experiments=[{'experiment':k,'support_min':r.get('min_endpoints',25),'computable':sum(f['score_S'] is not None for f in r['result']['factors'].values()),'pool_after_redundancy':len(r['result']['candidate_ranking']),'top':r['result']['top_candidates'],'boundary':r['result']['boundary_candidates']} for k,r in result['experiments'].items()]
    q1_csv(res/'experiments.csv',experiments)
    q1_csv(res/'sensitivity.csv',[{'branch':b['name'],'min_endpoints':b.get('min_endpoints',25),'computable':sum(f['score_S'] is not None for f in b['result']['factors'].values()),'pool_after_redundancy':len(b['result']['candidate_ranking']),'top':b['result']['top_candidates'],'boundary':b['result']['boundary_candidates']} for b in result['branches']])
    q1_csv(res/'loo_frequency.csv',[{'id':g,**r} for g,r in result['loo']['selection_frequency'].items()])
    q1_csv(res/'method_comparison.csv',result['comparisons_v2']);q1_csv(res/'c_fixed7d.csv',result['methods']['C_factors'])
    q1_csv(res/'extremes.csv',[{'id':e['id'],'label':e['label'],'unit':e['unit'],'period':period,**summary} for e in process['extremes']['series'] for period,summary in e['summaries'].items()])
    q1_csv(res/'process_support.csv',[{'channel':c['column'] if 'column' in c else c.get('id'),'details':c.get('support',c.get('summaries',{}))} for c in process.get('weather',{}).get('channels',[])])


def q1_execute(args):
    started=dt.datetime.now(dt.timezone.utc).isoformat();t0=time.monotonic();out=Path(args.output_dir);res=out/'results';res.mkdir(parents=True,exist_ok=True)
    data=Path(args.data_dir);dictionary=q1_features();features=dictionary['factors'];active=[f for f in features if f.get('computable')]
    inputs=[Path(__file__), Path(__file__).resolve().parent / 'data/features.json', Path(__file__).resolve().parent / 'data/reference.json',data/'run_summary_all.json',data/'attachment01_02/variable_dictionary.csv']+[data/'attachment01_02/station'/(t+'.csv.gz') for t in ('MET','SOIL','RADM','FLUX')]
    if args.v1_results:inputs += [Path(args.v1_results)/n for n in ('analysis_results.json','window_cache.json.gz','RUN_RECEIPT.json','configuration.json')]
    before=[q1_fp(p) for p in inputs];arrays,counts=q1_read_d0(data,features)
    for f in active:
        _,_,places,_=q2_inputs(arrays,f);f['source_grid_scale']=10**places;f['rank_scale']=10**places
    ledger=[];cache={}
    for ep,ch,cov in [('H-END',10,.9),('H-START',10,.9),('H-END',0,.9),('H-START',0,.9),('H-END',10,.95),('H-END',9,.9)]:
        print('V2_WINDOW_CACHE',ep,ch,cov,flush=True);cache[(ep,ch,cov)]=q1_cache(arrays,active,ep,ch,cov,ledger)
    q1_encode_exact_cache(cache,features);q1_dump(res/'feature_dictionary.json',dictionary)
    with gzip.open(res/'window_quality.csv.gz','wt',encoding='utf-8',newline='') as stream:
        wr=csv.DictWriter(stream,fieldnames=sorted({k for r in ledger for k in r}));wr.writeheader();wr.writerows(ledger)
    del ledger
    print('FIXED_MULTILEVEL_WEATHER_EXTREMES',flush=True)
    process=build_processes({k:v for k,v in arrays.items() if k!='_q2_feature_cache'},endpoint='H-END',cutoff_hour=10);q1_dump(res/'process_v2.json',process)
    primary_values=cache[('H-END',10,.9)];e1spec=[f for f in active if f['id'] in Q2_E1_IDS]
    print('E1_ALL_TWO_SIDED',flush=True);e1=summarize_one(e1spec,primary_values)
    print('E2_EXPANDED_ALL_TWO_SIDED',flush=True);primary=summarize_one(active,primary_values)
    assert all(e1['factors'][g]['score_S_exact']==primary['factors'][g]['score_S_exact'] for g in e1['factors'])
    old_policy=summarize_one(q2_original_specs(features),primary_values)
    experiments={'E0_V1_policy_reproduced':{'min_endpoints':25,'result':q1_trim_a(old_policy)},'E1_original10_two_sided':{'min_endpoints':25,'result':q1_trim_a(e1)},'E2_expanded_main':{'min_endpoints':25,'result':q1_trim_a(primary)}}
    branches=[]
    for ep,ch,cov in [('H-END',10,.9),('H-START',10,.9),('H-END',0,.9),('H-START',0,.9),('H-END',10,.95)]:
        r=primary if (ep,ch,cov)==('H-END',10,.9) else summarize_one(active,cache[(ep,ch,cov)])
        branches.append({'name':f'{ep}_C{ch:02d}_coverage{int(cov*100)}','endpoint':ep,'cutoff_hour':ch,'coverage':cov,'window_mode':'full','min_endpoints':25,'result':q1_trim_a(r)})
    reduced=copy.deepcopy(active)
    for f in reduced:f['window_hours']=[w for w in f['window_hours'] if w!=720]
    branches.append({'name':'drop30d_all_factors','window_mode':'drop30d','min_endpoints':25,'result':q1_trim_a(summarize_one(reduced,primary_values))})
    for cutoff in (20,15):
        for name,spec in [('E1',e1spec),('E2',active)]:
            r=summarize_one(spec,primary_values,min_endpoints=cutoff)
            experiments[f'{name}_uniform_min{cutoff}']={'min_endpoints':cutoff,'result':q1_trim_a(r)}
            if name=='E2':branches.append({'name':f'all_factors_min{cutoff}','min_endpoints':cutoff,'window_mode':'full','result':q1_trim_a(r)})
    for rho in (.85,.95):branches.append({'name':f'rho{rho}','min_endpoints':25,'result':q1_trim_a(summarize_one(active,primary_values,redundancy_threshold=rho))})
    common,common_info=q1_common_support(active,primary_values,primary)
    branches.append({'name':'common_history_support','min_endpoints':25,'common_support':common_info,'result':q1_trim_a(summarize_one(active,common))})
    loo=leave_one_year_out(active,primary_values)
    for k in loo['runs']:loo['runs'][k]=q1_trim_a(loo['runs'][k])
    c_rows,c_quality=q2_c_route(arrays,features,primary_values)
    baseline,validation,_=q2_validate_v1(args,cache,features,data,before)
    print('B_NATIVE_AND_ACTUAL_COMMON_REFITS',flush=True)
    methods=build_v2_comparisons(features,primary_values,e1,primary,c_rows,baseline=baseline,baseline_validation=validation,seed=args.seed)
    methods['methods']['C_factors']=c_rows;methods['methods']['C_quality']=c_quality
    result={'schema':'Q1_FORMAL_V2','research_revision_after_V1':True,'features':features,'primary_branch':{'name':'E2_H-END_C10_min25','endpoint':'H-END','cutoff_hour':10,'coverage':.9,'min_endpoints':25},'primary':primary,'experiments':experiments,'branches':branches,'loo':loo,'single_sided_hypotheses':q2_hypotheses(features,old_policy),'factor_evidence':q2_evidence(features,primary,arrays),'data_counts':counts,'common_support':common_info,'baseline_validation':validation,**methods,'scientific_acceptance':'pending V2 independent total review'}
    q1_dump(res/'configuration.json',{'seed':args.seed,'main_support_min':25,'support_sensitivities_all_factors':[20,15],'windows_hours':[24,72,168,360,720],'primary':'E2 C10 H-END cov90 min25','revision_after_V1':True,'not_independent_2026_evaluation':True,'fixed_process_variables':'23 channels plus five-level T/U, not selected by ranking','fixed_thresholds':'T2 lowP05/highP95 and WS10 highP95; five-year equal-weight inverse empirical CDF','B_counts':[2000,1000,500],'baseline_validation':validation})
    q1_dump(res/'analysis_results.json',result);q1_dump(res/'process_v2.json',process)
    with gzip.open(res/'window_cache.json.gz','wt',encoding='utf-8') as stream:json.dump({f'{ep}_C{ch}_{cov}':v for (ep,ch,cov),v in cache.items()},stream,ensure_ascii=False,allow_nan=False)
    q2_write_tables(res,result,process)
    after=[q1_fp(p) for p in inputs];assert before==after,'实际输入/代码漂移'
    receipt={'record_type':'actual_V2_analysis_process','started_utc':started,'finished_utc':dt.datetime.now(dt.timezone.utc).isoformat(),'elapsed_seconds':time.monotonic()-t0,'pid':os.getpid(),'argv':sys.argv,'environment':{'python':sys.version,'executable':'project .venv Python','executable_sha256':q1_sha(Path(sys.executable).resolve()),'numpy':np.__version__,'scipy':__import__('scipy').__version__,'platform':platform.platform(),'PYTHONDONTWRITEBYTECODE':os.environ.get('PYTHONDONTWRITEBYTECODE')},'inputs_before':before,'inputs_after':after,'inputs_unchanged':True,'outputs':[q1_fp(res/name) for name in Q2_OUTPUTS],'output_scope':'explicit analysis outputs; excludes all later receipts and rendering/writing files','status':'executed_pending_V2_independent_review'}
    q1_dump(res/'RUN_RECEIPT.json',receipt)
    print('V2_COMPLETE',json.dumps({'top':primary['top_candidates'],'boundary':primary['boundary_candidates'],'S':{g:r['score_S'] for g,r in primary['factors'].items()}},ensure_ascii=False),flush=True)


P3_FROZEN_V2 = json.loads((Path(__file__).resolve().parent / 'data' / 'reference.json').read_text(encoding='utf-8'))


from itertools import combinations as p3_combinations
import math as p3_math
import numpy as np

P3_MIN_ENDPOINTS = 25
P3_MIN_YEARS = 4
P3_RETAINED_VARIANCE = 0.85
P3_RANK_FACTOR = 1.0
P3_NEAR_ZERO_FACTOR = 64.0
P3_ADD_RTOL = 1e-10
P3_ADD_ATOL = 1e-12
P3_BRANCHES = (('main', 0), ('k_minus1', -1), ('k_plus1', 1))


def _p3_plain(value):
    if isinstance(value, dict):
        return {str(k): _p3_plain(v) for k, v in value.items()}
    if isinstance(value, (tuple, list, np.ndarray)):
        return [_p3_plain(v) for v in value]
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        number = float(value)
        return number if p3_math.isfinite(number) else None
    return value


def _p3_midrank(reference, value):

    x = np.asarray(reference, dtype=np.float64)
    if not p3_math.isfinite(float(value)) or not len(x) or not np.isfinite(x).all():
        return None
    return float((np.count_nonzero(x < value) + .5 * np.count_nonzero(x == value)) / len(x))


def _p3_ranks(values):

    x = np.asarray(values, dtype=np.float64)
    return [{'rank_min': int(1 + np.count_nonzero(x > v)),
             'rank_max': int(np.count_nonzero(x >= v))} for v in x]


def _p3_fit(reference_years, matrices):
    p = next(iter(matrices.values())).shape[1]
    blocks = [matrices[y] for y in reference_years]
    n = sum(len(x) for x in blocks)
    x = np.vstack(blocks)
    w = np.concatenate([np.full(len(block), 1. / (3 * len(block))) for block in blocks])
    assert len(reference_years) == 3 and np.isclose(w.sum(), 1., rtol=0., atol=1e-14)

    anchor = x[0].copy()
    mu = anchor + np.sum((x - anchor) * w[:, None], axis=0)
    d = x - mu
    variance = np.sum(w[:, None] * d * d, axis=0)
    base = {'reference_years': list(reference_years), 'training_rows': n,
            'training_rows_by_year': {y: len(matrices[y]) for y in reference_years},
            'row_weight_by_year': {y: 1. / (3 * len(matrices[y])) for y in reference_years},
            'weight_sum': float(w.sum()), 'mu': mu, 'variance': variance,
            'moment_convention': 'population weighted moment ddof=0; each reference year total weight 1/3'}
    if not np.isfinite(mu).all() or not np.isfinite(variance).all():
        return {**base, 'status': 'NA', 'reason': 'nonfinite_training_moments'}
    zero = np.flatnonzero(variance <= 0).tolist()
    if zero:
        return {**base, 'status': 'NA', 'reason': 'training_zero_variance_no_variables_deleted',
                'zero_variance_indices': zero}
    sigma = np.sqrt(variance)
    z = d / sigma
    covariance = (z * w[:, None]).T @ z
    covariance = (covariance + covariance.T) * .5
    raw_values, vectors = np.linalg.eigh(covariance)
    order = np.argsort(raw_values)[::-1]
    raw_values, vectors = raw_values[order], vectors[:, order]

    for j in range(p):
        pivot = int(np.argmax(np.abs(vectors[:, j])))
        if vectors[pivot, j] < 0:
            vectors[:, j] *= -1
    largest = float(max(0., raw_values[0]))
    tolerance = P3_RANK_FACTOR * max(n, p) * np.finfo(np.float64).eps * largest
    negative = np.flatnonzero(raw_values < 0).tolist()
    substantial_negative = np.flatnonzero(raw_values < -tolerance).tolist()
    base.update(sigma=sigma, standardized_covariance=covariance,
                raw_eigenvalues=raw_values, eigenvectors=vectors,
                eigenvalue_rank_tolerance=tolerance,
                negative_eigenvalue_indices=negative,
                substantial_negative_eigenvalue_indices=substantial_negative,
                weighted_standardized_center=np.sum(w[:, None] * z, axis=0),
                eigendecomposition_max_abs_error=float(np.max(np.abs(covariance @ vectors - vectors * raw_values))),
                orthogonality_max_abs_error=float(np.max(np.abs(vectors.T @ vectors - np.eye(p)))))
    if substantial_negative or largest <= 0:
        return {**base, 'status': 'NA', 'reason': 'invalid_nonpositive_covariance_spectrum'}
    eigenvalues = np.maximum(raw_values, 0.)
    effective_rank = int(np.count_nonzero(eigenvalues > tolerance))
    if effective_rank < 1:
        return {**base, 'status': 'NA', 'reason': 'effective_rank_zero', 'effective_rank': 0}

    total = float(np.sum(eigenvalues[:effective_rank]))
    explained = np.zeros_like(eigenvalues)
    explained[:effective_rank] = eigenvalues[:effective_rank] / total
    cumulative = np.cumsum(explained)
    main_k = min(effective_rank, int(np.searchsorted(cumulative[:effective_rank], P3_RETAINED_VARIANCE, side='left') + 1))
    base.update(status='ok', reason=None, eigenvalues=eigenvalues,
                effective_rank=effective_rank, main_k=main_k,
                explained_variance_ratio=explained,
                cumulative_explained_variance=cumulative,
                discarded_numerically_near_zero_eigenvalue_indices=list(range(effective_rank, p)),
                retained_variance_at_main_k=float(cumulative[main_k - 1]),
                loading_sign_rule='largest absolute component in each eigenvector is positive; no cross-model loading averaging')
    return base


def _p3_score(model, physical_values, k):
    x = np.asarray(physical_values, dtype=np.float64)
    x = np.atleast_2d(x)
    z = (x - model['mu']) / model['sigma']
    p = z.shape[1]
    pk = model['eigenvectors'][:, :k]
    lam = model['eigenvalues'][:k]
    t = z @ pk
    t2 = np.sum(t * t / lam, axis=1)
    residual = z - t @ pk.T
    q = np.sum(residual * residual, axis=1)
    ct = z * ((t / lam) @ pk.T)
    cq = residual * residual
    ct_error = np.abs(ct.sum(axis=1) - t2)
    cq_error = np.abs(cq.sum(axis=1) - q)
    if not (np.isfinite(z).all() and np.isfinite(t2).all() and np.isfinite(q).all()
            and np.isfinite(ct).all() and np.isfinite(cq).all()):
        raise FloatingPointError('nonfinite_PCA_score_or_contribution')
    if not np.allclose(ct.sum(axis=1), t2, rtol=P3_ADD_RTOL, atol=P3_ADD_ATOL):
        raise ArithmeticError('T2_contributions_do_not_add')
    if not np.allclose(cq.sum(axis=1), q, rtol=P3_ADD_RTOL, atol=P3_ADD_ATOL):
        raise ArithmeticError('Q_contributions_do_not_add')
    near_zero = P3_NEAR_ZERO_FACTOR * np.finfo(np.float64).eps * np.maximum(1., np.sum(z * z, axis=1))
    shares_t, shares_q = [], []
    for i in range(len(x)):
        shares_t.append((ct[i] / t2[i]).tolist() if t2[i] > near_zero[i] else None)
        shares_q.append((cq[i] / q[i]).tolist() if k < p and q[i] > near_zero[i] else None)
    return {'z': z, 'scores': t, 'residual': residual, 'T2': t2, 'Q': q, 'cT': ct, 'cQ': cq,
            'share_T2': shares_t, 'share_Q': shares_q,
            'rank_T2': [_p3_ranks(v) for v in ct], 'rank_Q': [_p3_ranks(v) for v in cq],
            'Q_informative': k < p,
            'Q_reason': None if k < p else 'full_dimension_residual_degenerate',
            'near_zero_share_threshold': near_zero,
            'max_contribution_sum_error_T2': float(ct_error.max(initial=0.)),
            'max_contribution_sum_error_Q': float(cq_error.max(initial=0.))}


def _p3_target_record(score):
    fields = ('z', 'scores', 'residual', 'T2', 'Q', 'cT', 'cQ', 'share_T2', 'share_Q', 'rank_T2', 'rank_Q',
              'near_zero_share_threshold')
    result = {key: score[key][0] for key in fields}
    result.update(Q_informative=score['Q_informative'], Q_reason=score['Q_reason'],
                  max_contribution_sum_error_T2=score['max_contribution_sum_error_T2'],
                  max_contribution_sum_error_Q=score['max_contribution_sum_error_Q'])
    return result


def _p3_equal_h(records, get_value):

    groups = {}
    for rec in records:
        value = get_value(rec)
        if value is not None and p3_math.isfinite(float(value)):
            groups.setdefault(rec['calibration_year'], []).append(float(value))
    by_year = [{'year': y, 'mean': float(np.mean(v)), 'min': min(v), 'max': max(v),
                'valid_configurations': len(v)} for y, v in sorted(groups.items())]
    values = [v for group in groups.values() for v in group]
    return {'mean': float(np.mean([r['mean'] for r in by_year])) if by_year else None,
            'min': min(values) if values else None, 'max': max(values) if values else None,
            'calibration_years': by_year, 'valid_configurations': len(values),
            'valid_calibration_years': len(by_year)}


def _p3_summary(configurations, ids, branch):
    records = [c for c in configurations if c['k_branch'] == branch]
    valid = [c for c in records if c['status'] == 'ok']
    result = {'k_branch': branch, 'expected_configurations': len(records),
              'valid_configurations': len(valid), 'NA_configurations': len(records) - len(valid),
              'metrics': {}, 'contributions': {}}
    for metric, ckey, skey in (('T2', 'cT', 'share_T2'), ('Q', 'cQ', 'share_Q')):
        agg = _p3_equal_h(valid, lambda c: c['positions'][metric])
        by_h = [{**v, 'position_mean': v['mean'], 'position_min': v['min'], 'position_max': v['max']}
                for v in agg['calibration_years']]
        result['metrics'][metric] = {'position_mean': agg['mean'], 'position_min': agg['min'],
                                     'position_max': agg['max'], 'calibration_years': by_h,
                                     'valid_configurations': agg['valid_configurations'],
                                     'valid_calibration_years': agg['valid_calibration_years'],
                                     'NA_configurations': len(records) - agg['valid_configurations'],
                                     'interpretation': 'held-out empirical middle-rank position; not p/probability/accuracy'}
        rows = []
        for j, factor_id in enumerate(ids):
            pos = _p3_equal_h(valid, lambda c: c['positions']['contribution_' + metric][j])
            share = _p3_equal_h(valid, lambda c: c['target'][skey][j] if c['target'][skey] is not None else None)
            zstat = _p3_equal_h(valid, lambda c: c['target']['z'][j])

            rank_records = [c for c in valid if metric != 'Q' or c['target']['Q_informative']]
            ranks = [c['target']['rank_' + metric][j] for c in rank_records]
            negative = sum(c['target'][ckey][j] < 0 for c in rank_records)
            rows.append({'id': factor_id,
                         'normalized_share_mean': share['mean'], 'normalized_share_min': share['min'],
                         'normalized_share_max': share['max'], 'share_valid_configurations': share['valid_configurations'],
                         'share_valid_calibration_years': share['valid_calibration_years'],
                         'position_mean': pos['mean'], 'position_min': pos['min'], 'position_max': pos['max'],
                         'position_calibration_years': pos['calibration_years'],
                         'rank_min': min((v['rank_min'] for v in ranks), default=None),
                         'rank_max': max((v['rank_max'] for v in ranks), default=None),
                         'negative_configurations': negative, 'valid_configurations': len(rank_records),
                         'z_mean': zstat['mean'], 'z_min': zstat['min'], 'z_max': zstat['max'],
                         'z_positive_configurations': sum(c['target']['z'][j] > 0 for c in valid),
                         'z_negative_configurations': sum(c['target']['z'][j] < 0 for c in valid),
                         'share_interpretation': 'signed cT/T2' if metric == 'T2' else 'cQ/Q',
                         'raw_contributions_not_averaged': True})
        result['contributions'][metric] = rows
    result['k_values'] = sorted(set(c['k'] for c in valid))
    return result


def _p3_structure_summary(models, ids):
    valid = [m for m in models if m['status'] == 'ok']
    rows = []
    for i, j in p3_combinations(range(len(ids)), 2):
        xs = [float(m['standardized_covariance'][i, j]) for m in valid]
        rows.append({'id_a': ids[i], 'id_b': ids[j], 'training_fits': len(xs),
                     'correlation_mean': float(np.mean(xs)) if xs else None,
                     'correlation_min': min(xs) if xs else None, 'correlation_max': max(xs) if xs else None,
                     'positive_fits': sum(x > 0 for x in xs), 'negative_fits': sum(x < 0 for x in xs),
                     'interpretation': 'equal-year training standardized covariance; each distinct R counted once; no PC loading average'})
    return rows


def build_pca_v3(feature_ids, matrices_by_year, targets_by_scale, feature_meta=None):


    ids = [str(v) for v in feature_ids]
    if len(ids) < 2 or len(set(ids)) != len(ids):
        raise ValueError('PCA feature IDs must be unique and contain at least two factors')
    p = len(ids)
    result = {'schema': 'Q1_PCA_V3.1', 'feature_ids': ids, 'feature_meta': feature_meta,
              'settings': {'scales_hours': [168, 720], 'min_endpoints_per_year': P3_MIN_ENDPOINTS,
                           'minimum_years': P3_MIN_YEARS, 'reference_years_per_fit': 3,
                           'retained_variance_main': P3_RETAINED_VARIANCE,
                           'retained_variance_role': 'operational convention fixed for this post-V2 exploration, not accuracy or official standard',
                           'moment_convention': 'population weighted variance/covariance ddof=0; row weight 1/(3*n_y)',
                           'rank_rule': 'eigenvalue > max(n_training_rows,p)*float64_eps*lambda_max',
                           'rank_factor': P3_RANK_FACTOR, 'small_negative_rule': 'record and clip [-rank_tol,0) to zero; below -rank_tol is NA',
                           'zero_variance_rule': 'whole reference fit NA; no variable deleted',
                           'Q_degenerate_rule': 'k==p: raw roundoff residual recorded, Q position/contribution position/share/rank are NA',
                           'share_near_zero_rule': 'total <=64*eps*max(1,sum(z^2)): share NA, raw total/contribution retained',
                           'contribution_sum_tolerance': {'rtol': P3_ADD_RTOL, 'atol': P3_ADD_ATOL},
                           'aggregation': 'equal R combinations within calibration year, then equal calibration years; raw totals and contributions never averaged across models',
                           'sensitivity': 'only k-1/k+1 inside [1,effective_rank]; reuse same decomposition',
                           'scope': 'auxiliary linear PCA after seeing V2; no independent disaster evaluation'},
              'scales': {}}
    for hours in (168, 720):
        source = matrices_by_year.get(hours, matrices_by_year.get(str(hours), {}))
        target = targets_by_scale.get(hours, targets_by_scale.get(str(hours)))
        matrices, excluded = {}, []
        for year, rows in sorted(source.items(), key=lambda q: int(q[0])):
            y = str(year)
            if int(y) not in range(2021, 2026):
                raise ValueError('Only 2021--2025 historical years may enter PCA fitting/calibration')
            x = np.asarray(rows, dtype=np.float64)
            if x.ndim != 2 or x.shape[1] != p:
                raise ValueError('Each historical PCA matrix must preserve the complete feature vector/order')
            if not np.isfinite(x).all():
                raise ValueError('PCA common-support matrices must contain no missing/nonfinite values')
            if len(x) < P3_MIN_ENDPOINTS:
                excluded.append({'year': y, 'rows': len(x), 'reason': 'common_endpoints_below_25'})
            else:
                matrices[y] = x
        years = tuple(sorted(matrices, key=int))
        scale = {'hours': hours, 'feature_ids': ids, 'years': list(years),
                 'rows_by_year': {y: len(matrices[y]) for y in years}, 'excluded_years': excluded,
                 'status': 'pending', 'models': [], 'configurations': [], 'summaries': {}, 'structure_summary': []}
        result['scales'][str(hours)] = scale
        if len(years) < P3_MIN_YEARS:
            scale.update(status='NA', reason='fewer_than_four_common_support_years')
            continue
        tx = np.asarray(target, dtype=np.float64) if target is not None else np.array([])
        if tx.shape != (p,) or not np.isfinite(tx).all():
            scale.update(status='NA', reason='target_complete_physical_vector_unavailable')
            continue
        scale.update(status='ok', reason=None, target_physical=tx)
        model_cache = {}
        for reference_years in p3_combinations(years, 3):
            model_id = f'h{hours}_R' + '-'.join(reference_years)
            try:
                with np.errstate(over='raise', invalid='raise', divide='raise'):
                    model = _p3_fit(reference_years, matrices)
            except (FloatingPointError, np.linalg.LinAlgError) as error:
                model = {'status': 'NA', 'reason': type(error).__name__ + ':' + str(error),
                         'reference_years': list(reference_years)}
            model['model_id'] = model_id
            model_cache[reference_years] = model
            scale['models'].append(model)
        for h in years:
            for reference_years in p3_combinations(tuple(y for y in years if y != h), 3):
                assert h not in reference_years
                model = model_cache[reference_years]
                for branch, delta in P3_BRANCHES:
                    rec = {'id': f'h{hours}_h{h}_R' + '-'.join(reference_years) + '_' + branch,
                           'hours': hours, 'reference_years': list(reference_years),
                           'calibration_year': h, 'model_id': model['model_id'], 'k_branch': branch,
                           'calibration_rows': len(matrices[h]), 'status': 'pending'}
                    scale['configurations'].append(rec)
                    if model['status'] != 'ok':
                        rec.update(status='NA', reason='reference_fit_NA:' + model['reason'], k=None)
                        continue
                    k = model['main_k'] + delta
                    rec['k'] = k
                    if not 1 <= k <= model['effective_rank']:
                        rec.update(status='NA', reason='adjacent_k_outside_valid_rank_range')
                        continue
                    try:
                        with np.errstate(over='raise', invalid='raise', divide='raise'):
                            ts = _p3_score(model, tx, k)
                            hs = _p3_score(model, matrices[h], k)
                        q_valid = bool(ts['Q_informative'])
                        rec.update(status='ok', reason=None, target=_p3_target_record(ts),
                                   calibration={key: hs[key] for key in ('z', 'scores', 'residual', 'T2', 'Q', 'cT', 'cQ', 'share_T2', 'share_Q')},
                                   positions={'T2': _p3_midrank(hs['T2'], ts['T2'][0]),
                                              'Q': _p3_midrank(hs['Q'], ts['Q'][0]) if q_valid else None,
                                              'contribution_T2': [_p3_midrank(hs['cT'][:, j], ts['cT'][0, j]) for j in range(p)],
                                              'contribution_Q': [_p3_midrank(hs['cQ'][:, j], ts['cQ'][0, j]) if q_valid else None for j in range(p)]},
                                   checks={'reference_and_calibration_disjoint': h not in reference_years,
                                           'target_not_in_fit': True,
                                           'T2_additivity_target_max_error': ts['max_contribution_sum_error_T2'],
                                           'Q_additivity_target_max_error': ts['max_contribution_sum_error_Q'],
                                           'T2_additivity_calibration_max_error': hs['max_contribution_sum_error_T2'],
                                           'Q_additivity_calibration_max_error': hs['max_contribution_sum_error_Q']})
                    except (FloatingPointError, ArithmeticError) as error:
                        rec.update(status='NA', reason=type(error).__name__ + ':' + str(error))
        scale['summaries'] = {name: _p3_summary(scale['configurations'], ids, name) for name, _ in P3_BRANCHES}
        scale['structure_summary'] = _p3_structure_summary(scale['models'], ids)
        scale['fit_cache_count'] = len(model_cache)
        scale['decomposition_reused_for_calibration_years_and_adjacent_k'] = True
        scale['status'] = 'ok' if scale['summaries']['main']['valid_configurations'] else 'NA'
        if scale['status'] == 'NA':
            scale['reason'] = 'all_reference_fits_or_main_configurations_NA'
    return _p3_plain(result)


def pca_v3_selftest():

    rng = np.random.default_rng(20261026)
    mats = {str(y): rng.normal(size=(25 + (y % 3), 3)) * np.array([2., .5, 4.]) + np.array([10., -2., 20.]) for y in range(2021, 2026)}
    model = _p3_fit(('2021', '2022', '2023'), mats)
    assert model['status'] == 'ok'
    assert all(abs(model['row_weight_by_year'][y] * len(mats[y]) - 1/3) < 1e-14 for y in model['reference_years'])
    assert np.allclose(np.diag(model['standardized_covariance']), 1., atol=1e-13)
    assert model['orthogonality_max_abs_error'] < 1e-12
    assert model['eigendecomposition_max_abs_error'] < 1e-12
    target = np.array([13., -1., 23.])
    for k in range(1, model['effective_rank'] + 1):
        sc = _p3_score(model, target, k)
        assert np.allclose(sc['cT'].sum(axis=1), sc['T2'])
        assert np.allclose(sc['cQ'].sum(axis=1), sc['Q'])
        if k == 3:
            assert sc['share_Q'] == [None] and not sc['Q_informative']
    with_target_changed = build_pca_v3(['a', 'b', 'c'], {168: mats, 720: mats}, {168: target + 100, 720: target + 100})
    original = build_pca_v3(['a', 'b', 'c'], {168: mats, 720: mats}, {168: target, 720: target})
    assert original['scales']['168']['models'] == with_target_changed['scales']['168']['models']

    altered = {y: v.copy() for y, v in mats.items()}
    altered['2025'] += 1000
    assert np.array_equal(_p3_fit(('2021', '2022', '2023'), mats)['mu'], _p3_fit(('2021', '2022', '2023'), altered)['mu'])
    constant = {str(y): np.column_stack([np.arange(25), np.ones(25), np.arange(25)**2]) for y in range(2021, 2026)}
    assert _p3_fit(('2021', '2022', '2023'), constant)['status'] == 'NA'

    manual = {'mu': np.zeros(2), 'sigma': np.ones(2), 'eigenvectors': np.array([[1., 1.], [1., -1.]]) / np.sqrt(2), 'eigenvalues': np.array([1.5, .5])}
    neg = _p3_score(manual, np.array([2., -1.]), 1)
    assert neg['cT'][0, 1] < 0 and np.isclose(neg['cT'].sum(), neg['T2'][0])
    import json
    json.dumps(original, allow_nan=False)
    return {'status': 'pass', 'checks': ['year_equal_population_moments', 'reference_only_standardization',
                                        'eigen_orthogonality', 'signed_contribution_additivity',
                                        'target_does_not_change_fit', 'held_out_year_does_not_change_excluding_R_fit',
                                        'zero_training_variance_is_NA', 'full_dimension_Q_is_uninformative',
                                        'strict_JSON'],
            'scope': 'bounded author self-check; not independent statistical or scientific acceptance'}


def p3_inputs(windows):
    ids=P3_FROZEN_V2['feature_ids'];specs={f['id']:f for f in P3_FROZEN_V2['features']}
    support_rows=[];n_by_year={};n_two_scale={};offsets={};matrices={str(h):{} for h in (168,720)};exact={str(h):{} for h in (168,720)}
    target_valid=all(windows[g]['target'].get(str(h)) is not None for g in ids for h in (168,720))
    targets={str(h):[float(windows[g]['target'][str(h)]) if windows[g]['target'].get(str(h)) is not None else None for g in ids] for h in (168,720)}
    for year in map(str,range(2021,2026)):
        common=[];dual=0
        for offset in range(-15,16):
            o=str(offset);two_missing=[f'{g}/{h}' for g in ids for h in (168,720) if windows[g]['years'][year][o].get(str(h)) is None]
            wide_missing=[f'{g}/{h}' for g in ids for h in specs[g]['window_hours'] if windows[g]['years'][year][o].get(str(h)) is None]
            if not two_missing:dual+=1
            eligible=not wide_missing
            if eligible:common.append(o)
            endpoint=dt.datetime(int(year),8,26,10,tzinfo=Q1_TZ)+dt.timedelta(days=offset)
            support_rows.append({'year':year,'offset_days':offset,'endpoint_id':f'{year}:{offset:+d}','endpoint_UTC8':endpoint.isoformat(),'two_scale_complete':not two_missing,'V2_all_applicable_windows_complete':eligible,'included_in_PCA':eligible,'two_scale_missing_fields':two_missing,'V2_all_window_missing_fields':wide_missing,'exclusion_reason':None if eligible else 'V2_all_window_common_support_conservatively_reused'})
        n_by_year[year]=len(common);n_two_scale[year]=dual;offsets[year]=common
        for h in (168,720):
            exact[str(h)][year]=[[str(windows[g]['years'][year][o][str(h)]) for g in ids] for o in common]
            matrices[str(h)][year]=[[float(windows[g]['years'][year][o][str(h)]) for g in ids] for o in common]
    supported=[y for y,n in n_by_year.items() if n>=25]
    for row in support_rows:
        if row['year'] not in supported:
            row['included_in_PCA']=False
            if row['V2_all_applicable_windows_complete']:row['exclusion_reason']='year_common_endpoints_below25'
    matrices={h:{y:v for y,v in data.items() if y in supported} for h,data in matrices.items()}
    support={'rule':'V2 all applicable windows common support, conservatively reused without adding newly eligible two-scale rows','n_by_year':n_by_year,'two_scale_available_n_by_year':n_two_scale,'included_offsets_by_year':{y:offsets[y] for y in supported},'included_endpoint_ids_by_year':{y:[f'{y}:{int(o):+d}' for o in offsets[y]] for y in supported},'valid_years':supported,'minimum_endpoints_per_year':25,'minimum_years':4,'included_rows_total':sum(n_by_year[y] for y in supported),'same_endpoints_for_both_scales':True,'extra_two_scale_available_not_used':[r for r in support_rows if r['two_scale_complete'] and not r['included_in_PCA']],'target_complete':target_valid,'note':'Main PCA freezes approved V2 all-window shared support. 7/30-day-only availability is disclosed, not a sensitivity fit; no target-dependent selection.'}
    inputs={'feature_ids':ids,'matrices_by_year':matrices,'targets_by_scale':targets,'physical_values_exact':exact,'target_physical_values_exact':{str(h):[str(windows[g]['target'][str(h)]) if windows[g]['target'].get(str(h)) is not None else None for g in ids] for h in (168,720)},'support':support}
    return inputs,support_rows

def p3_tables(res,result):
    summary=[];contrib=[];structure=[];models=[]
    for h,scale in result['scales'].items():
        for branch,s in scale['summaries'].items():
            for metric,row in s['metrics'].items():summary.append({'hours':h,'k_branch':branch,'metric':metric,'k_values':s['k_values'],**row})
            for metric,rows in s['contributions'].items():
                for row in rows:contrib.append({'hours':h,'k_branch':branch,'metric':metric,**row})
        structure.extend({'hours':h,**row} for row in scale['structure_summary'])
        models.extend({'hours':h,**{k:v for k,v in m.items() if k not in ('eigenvectors','standardized_covariance')}} for m in scale['models'])
    q1_csv(res/'pca_summary.csv',summary);q1_csv(res/'pca_contributions.csv',contrib);q1_csv(res/'pca_structure.csv',structure);q1_csv(res/'pca_models.csv',models)


import tempfile
import shutil

_QF_REF_CACHE=None
_QF_MANIFEST=None
_QF_REUSED_FAMILIES=[]


def qf_read_json(path):
    p=Path(path)
    text=gzip.decompress(p.read_bytes()).decode('utf-8') if p.suffix=='.gz' else p.read_text(encoding='utf-8')
    return json.loads(text)


def qf_write_gzip(path,value):
    data=(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n').encode()
    Path(path).write_bytes(gzip.compress(data,compresslevel=9,mtime=0))


def qf_physical_fingerprint(features,values,seed,n_perm,n_boot,n_perm_boot):
    payload={'seed':seed,'counts':[n_perm,n_boot,n_perm_boot],'features':[]}
    for f in features:
        g=f['id'];scale=f['rank_scale'];record=values[g]
        convert=lambda row:{str(h):str(Fraction(row[str(h)],scale)) if row.get(str(h)) is not None else None for h in f['window_hours']}
        payload['features'].append({'id':g,'windows':f['window_hours'],'unit':f['unit'],'formula':f['formula'],'target':convert(record['target']),'seasonal_years':{y:{md:convert(row) for md,row in rows.items()} for y,rows in record['seasonal_years'].items()}})
    return hashlib.sha256(json.dumps(payload,sort_keys=True,ensure_ascii=False,separators=(',',':')).encode()).hexdigest()


def qf_run_ref(features,values,seed=20261024,*,label='B-aligned-corrected',n_perm=2000,n_boot=1000,n_perm_boot=500):
    if _QF_REF_CACHE is None:
        return run_ref(features,values,seed=seed,label=label,n_perm=n_perm,n_boot=n_boot,n_perm_boot=n_perm_boot)
    family='|'.join(f['id'] for f in features)
    expected=_QF_MANIFEST['modex_families'].get(family)
    if expected is None:raise ValueError('Frozen 参考方法 cache does not contain the exact ordered family')
    actual=qf_physical_fingerprint(features,values,seed,n_perm,n_boot,n_perm_boot)
    if actual!=expected['physical_window_sha256']:raise ValueError('参考方法 physical inputs/definitions/counts differ: full recomputation required')
    answer=copy.deepcopy(_QF_REF_CACHE['methods'][expected['method_key']])
    answer['_final_reused_modex']=True
    _QF_REUSED_FAMILIES.append({'method':expected['method_key'],'physical_window_sha256':actual,'new_draws':0})
    return answer


def qf_pca_inputs_from_analysis(analysis,cache):
    ids=[r['id'] for r in analysis['primary']['candidate_ranking']]
    if ids!=P3_FROZEN_V2['feature_ids']:
        raise ValueError('Reconstructed candidate pool differs from frozen Q1; do not silently change PCA variables')
    specs={f['id']:f for f in analysis['features']};windows={}
    for g in ids:
        scale=specs[g]['rank_scale'];rec=cache[g]
        convert=lambda row:{str(h):Fraction(row[str(h)],scale) if row.get(str(h)) is not None else None for h in specs[g]['window_hours']}
        windows[g]={'target':convert(rec['target']),'years':{str(y):{str(o):convert(rec['years'][str(y)][str(o)]) for o in range(-15,16)} for y in range(2021,2026)}}
    return p3_inputs(windows)


def qf_verify(results_dir):
    root=Path(results_dir);manifest=qf_read_json(root/'frozen_manifest.json')
    for name,expected in manifest['data_files'].items():
        p=root/name
        if not p.is_file() or q1_sha(p)!=expected['sha256']:raise ValueError('Packed frozen result missing or changed: '+name)
    a=qf_read_json(root/'analysis_results.json.gz');process=qf_read_json(root/'process_results.json.gz');pca=qf_read_json(root/'pca_results.json.gz')
    assert len(a['primary']['candidate_ranking'])==len(pca['feature_ids'])==13
    assert [r['id'] for r in a['primary']['candidate_ranking']]==pca['feature_ids']
    assert process['metadata']['primary_cutoff']=='2026-08-26T10:00:00+08:00'
    assert pca['support']['included_rows_total']==148 and pca['support']['same_endpoints_for_both_scales']
    print('FROZEN_RESULTS_VERIFIED: complete single-factor/process/comparison/sensitivity/PCA inputs; no analysis or 参考方法 executed')
    return a,process,pca


def qf_recompute(args):
    global _QF_REF_CACHE,_QF_MANIFEST
    started=dt.datetime.now(dt.timezone.utc).isoformat();out=Path(args.output_dir);res=out/'results';res.mkdir(parents=True,exist_ok=True)
    files=[Path(__file__), Path(__file__).resolve().parent / 'data/features.json', Path(__file__).resolve().parent / 'data/reference.json',Path(args.data_dir)/'run_summary_all.json',Path(args.data_dir)/'attachment01_02/variable_dictionary.csv']+[Path(args.data_dir)/'attachment01_02/station'/(t+'.csv.gz') for t in ('MET','SOIL','RADM','FLUX')]
    if args.reuse_ref:
        root=Path(args.reuse_ref);_QF_REF_CACHE,_,_=qf_verify(root);_QF_MANIFEST=qf_read_json(root/'frozen_manifest.json')
        for relative,expected in _QF_MANIFEST['frozen_inputs'].items():
            if q1_sha(Path(args.data_dir)/relative)!=expected['sha256']:raise ValueError('D0 version differs; frozen 参考方法 must not be reused: '+relative)
        files += [root/'frozen_manifest.json',root/'analysis_results.json.gz']
    before=[q1_fp(p) for p in files]

    with tempfile.TemporaryDirectory(prefix='q1-compute-',dir=out) as temp:
        stage=Path(temp);q1_execute(argparse.Namespace(data_dir=args.data_dir,output_dir=str(stage),v1_results=None,seed=20261024))
        oldres=stage/'results';analysis=qf_read_json(oldres/'analysis_results.json');process=qf_read_json(oldres/'process_v2.json')
        with gzip.open(oldres/'window_cache.json.gz','rt',encoding='utf-8') as stream:cache=json.load(stream)['H-END_C10_0.9']
        inputs,support_rows=qf_pca_inputs_from_analysis(analysis,cache)
        feature_meta=[{'id':f['id'],'label':f['label'],'unit':f['unit']} for f in P3_FROZEN_V2['features']]
        pca=build_pca_v3(inputs['feature_ids'],inputs['matrices_by_year'],inputs['targets_by_scale'],feature_meta=feature_meta)
        pca.update(support=inputs['support'],metadata={'version':'final complete Q1','cutoff':'2026-08-26T10:00:00+08:00','endpoint':'H-END','time_assumption':'hourly means and same channel cross-year; unconfirmed source facts remain conditional','V2_main_ranking_unchanged':True,'independent_2026_evaluation':False,'scientific_acceptance':'new execution pending independent review','input_provenance':{'kind':'complete_D0_Q1_path','Modex_reused':bool(args.reuse_ref)}},v2_baseline={k:analysis['primary'][k] for k in ('candidate_ranking','top_candidates','boundary_candidates','available_candidate_count')},v2_factor_evidence=[e for e in analysis['factor_evidence'] if e['id'] in inputs['feature_ids']])
        for scale in pca['scales'].values():
            for model in scale['models']:model['training_endpoint_ids_by_year']={y:inputs['support']['included_endpoint_ids_by_year'][y] for y in model['reference_years']}
            for conf in scale['configurations']:conf['calibration_endpoint_ids']=inputs['support']['included_endpoint_ids_by_year'][conf['calibration_year']]
        qf_write_gzip(res/'analysis_results.json.gz',analysis);qf_write_gzip(res/'process_results.json.gz',process);qf_write_gzip(res/'pca_results.json.gz',pca)
        q1_dump(res/'pca_inputs.json',inputs);q1_csv(res/'common_support.csv',support_rows);p3_tables(res,pca)
        q1_dump(res/'analysis_configuration.json',qf_read_json(oldres/'configuration.json'))
        shutil.copyfile(oldres/'window_cache.json.gz',res/'window_cache.json.gz')
        q1_dump(res/'pca_configuration.json',{'support':inputs['support'],'settings':pca['settings'],'feature_ids':pca['feature_ids']})
        for name in ('all_factors.csv','primary_windows.csv','single_sided_hypotheses.csv','experiments.csv','sensitivity.csv','loo_frequency.csv','method_comparison.csv','c_fixed7d.csv','extremes.csv','process_support.csv'):shutil.copyfile(oldres/name,res/name)
    after=[q1_fp(p) for p in files];assert before==after,'Input/source changed during execution'
    outputs=[p for p in sorted(res.iterdir()) if p.is_file() and p.name not in ('RUN_RECEIPT.json','FIGURE_RUN_RECEIPT.json')]
    q1_dump(res/'RUN_RECEIPT.json',{'record_type':'actual_final_integrated_Q1_process','started_utc':started,'finished_utc':dt.datetime.now(dt.timezone.utc).isoformat(),'argv':sys.argv,'inputs_before':before,'inputs_after':after,'inputs_unchanged':True,'outputs':[q1_fp(p) for p in outputs],'scope':'D0 feature/A/process/sensitivity/C/PCA recomputed; 参考方法 reused after exact matching' if args.reuse_ref else 'complete D0 Q1 including all 参考方法 computations','reused_Modex_families':_QF_REUSED_FAMILIES,'new_Modex_computation_executed':not bool(args.reuse_ref),'scientific_acceptance':'pending independent review of this execution'})
    print('FINAL_Q1_COMPLETE; 参考方法 '+('REUSED, no new draws' if args.reuse_ref else 'actually recomputed'))


def qf_main():
    here=Path(__file__).resolve().parent
    parser=argparse.ArgumentParser(description='Complete Q1: approved analysis/process/A-B-C/sensitivity and PCA, plus explicit frozen-result verification')
    parser.add_argument('--mode',choices=('verify','recompute'),default='verify')
    parser.add_argument('--results-dir',default=str(here/'results'))
    parser.add_argument('--data-dir',help='D0 cleaned_all directory, required only for recompute')
    parser.add_argument('--output-dir',default=str(here/'recomputed'))
    parser.add_argument('--reuse-ref',metavar='RESULTS_DIR',help='Reuse exact frozen 参考方法 families only after D0 and physical window checks; no new draws')
    args=parser.parse_args()
    if args.mode=='verify':qf_verify(args.results_dir)
    else:
        if not args.data_dir:parser.error('--data-dir is required for recompute')
        qf_recompute(args)


QX_PRECIP_IDS=('precip_sum','precip_max','rain_hours','api')
QX_NEW_IDS=tuple('d1_era5_'+g for g in QX_PRECIP_IDS)


def qx_precip_specs(baseline_features):
    index={f['id']:f for f in baseline_features};out=[]
    labels={'precip_sum':'ERA5格点累计总降水','precip_max':'ERA5格点最大小时总降水量','rain_hours':'ERA5格点有降水小时数（≥0.1 mm）','api':'ERA5格点前期总降水指数（K=0.9）'}
    formulas={'precip_sum':'sum(P_i)','precip_max':'max(P_i)/1 h','rain_hours':'sum(1{P_i>=0.1 mm}) * 1 h','api':'Cutoff-aligned 24h blocks; A_0=0 at each window start; A_d=P_d+0.9*A_(d-1); finite-window endpoint value'}
    for original in QX_PRECIP_IDS:
        f=copy.deepcopy(index[original]);f.update(id='d1_era5_'+original,name=labels[original],label=labels[original],original_factor_id=original,formula=formulas[original],window_hours=[24,72,168,360,720],windows_h=[24,72,168,360,720],representative_hours=168,computable=True,eligible_for_top5=True,status='product_defined_proxy',status_label='来源定义明确的格点代理',status_reason='Open-Meteo fixed ERA5 previous-hour total precipitation; independent D1, not restoration of station Precip',source_grid_scale=1,rank_scale=1,raw_columns=[{'table':'D1_ERA5','column':'precipitation_mm','qc_column':None,'source_qc_identity':'not_provided','support':'[interval_end-1h,interval_end)'}],conditional_assumptions=['ERA5 gridded total precipitation represents an environmental proxy at the returned grid point; no station/source-slope identity claimed','API finite-window initial state zero and cutoff-aligned 24h blocks, fixed before scoring'],limitations=['ERA5 total precipitation includes liquid-equivalent rain/snow; not a station rain gauge','No source QC flags supplied; completeness does not imply no error','No interpolation, zero-filling, elevation correction or blending into D0','0.1 mm is an analysis threshold, not instrument accuracy'],physical_family='era5_total_precipitation_memory' if original=='api' else 'era5_total_precipitation_amount_or_intensity',historical_precision={'encoding':'exact Decimal/Fraction; lossless integer common denominator','measurement_accuracy_claim':False})
        f['window_applicability']={str(h):'fixed and identically used for target/pseudo' for h in f['window_hours']}
        out.append(f)
        MECHANISM[f['id']]=MECHANISM[original];ORIGINAL_NAMES[f['id']]=labels[original]
    return out


def qx_read_precip(path):
    values={};products=set();grids=set();qc=set()
    with Path(path).open(encoding='utf-8-sig',newline='') as stream:
        for row in csv.DictReader(stream):
            t=dt.datetime.fromisoformat(row['interval_end']);start=dt.datetime.fromisoformat(row['interval_start'])
            if t.utcoffset()!=dt.timedelta(hours=8) or t.minute or t.second or t.microsecond or t-start!=dt.timedelta(hours=1):raise ValueError('D1 support must be exact aware UTC+8 preceding hour')
            if t in values:raise ValueError('Duplicate D1 interval')
            token=row['precipitation_mm'];p=None if token=='' else Fraction(decimal.Decimal(token))
            if p is not None and p<0:raise ValueError('Negative D1 precipitation')
            values[t]=p;products.add(row['source_product']);grids.add((row['grid_latitude'],row['grid_longitude'],row['grid_elevation_m']));qc.add(row['source_qc_status'])
    if len(products)!=1 or len(grids)!=1 or qc!={'not_provided'}:raise ValueError('D1 product/grid/QC identity inconsistent across target/history')
    return values,{'product':next(iter(products)),'returned_grid':list(next(iter(grids))),'source_qc_status':'not_provided','rows':len(values),'valid_rows':sum(v is not None for v in values.values()),'years':sorted({t.year for t in values}),'support':'[interval_end-1h,interval_end)','aggregation':'hourly total precipitation, rain/snow liquid equivalent','is_station_observation':False}


def qx_precip_window(series,end,hours,original):
    vals=[series.get(end-dt.timedelta(hours=i)) for i in range(hours-1,-1,-1)]
    if any(v is None for v in vals):return None
    if original=='precip_sum':return sum(vals,Fraction(0))
    if original=='precip_max':return max(vals)
    if original=='rain_hours':return Fraction(sum(v>=Fraction(1,10) for v in vals))
    if original=='api':
        state=Fraction(0)
        for i in range(0,hours,24):state=Fraction(9,10)*state+sum(vals[i:i+24],Fraction(0))
        return state
    raise ValueError(original)


def qx_precip_cache(series,specs,cutoff_hour):

    scale=1
    for value in series.values():
        if value is not None:scale=math.lcm(scale,value.denominator)
    raw=np.zeros(Q1_GRID_LENGTH,dtype=np.int64);valid=np.zeros(Q1_GRID_LENGTH,dtype=np.int64)
    for stamp,value in series.items():
        i=q1_hour(stamp.isoformat())
        if 0<=i<Q1_GRID_LENGTH and value is not None:
            n=value*scale;assert n.denominator==1
            raw[i]=n.numerator;valid[i]=1
    total=np.concatenate(([0],np.cumsum(raw)));count=np.concatenate(([0],np.cumsum(valid)))
    wet=np.concatenate(([0],np.cumsum((raw*10>=scale)&valid.astype(bool))))
    def all_windows(end):
        out={g:{} for g in QX_PRECIP_IDS}
        for hours in (24,72,168,360,720):
            lo,hi=end-hours+1,end+1;h=str(hours)
            if lo<0 or hi>Q1_GRID_LENGTH or count[hi]-count[lo]!=hours:
                for g in out:out[g][h]=None
                continue
            out['precip_sum'][h]=Fraction(int(total[hi]-total[lo]),scale)
            out['precip_max'][h]=Fraction(int(raw[lo:hi].max()),scale)
            out['rain_hours'][h]=Fraction(int(wet[hi]-wet[lo]))
            state=Fraction(0)
            for j in range(lo,hi,24):state=Fraction(9,10)*state+Fraction(int(total[j+24]-total[j]),scale)
            out['api'][h]=state
        return out
    first=dt.date(2021,5,31);last=dt.date(2021,10,15)
    monthdays=[(first+dt.timedelta(days=d)).strftime('%m-%d') for d in range((last-first).days+1)]
    target=all_windows(q1_hour(f'2026-08-26T{cutoff_hour:02d}:00:00+08:00'))
    annual={str(y):{md:all_windows(q1_hour(f'{y}-{md}T{cutoff_hour:02d}:00:00+08:00')) for md in monthdays} for y in A_DEFAULT_YEARS}
    values={}
    for f in specs:
        g=f['original_factor_id'];data={'target':target[g],'years':{},'seasonal_years':{}}
        for y in A_DEFAULT_YEARS:
            seasonal={md:annual[str(y)][md][g] for md in monthdays};data['seasonal_years'][str(y)]=seasonal
            data['years'][str(y)]={str(o):seasonal[(dt.date(y,8,26)+dt.timedelta(days=o)).strftime('%m-%d')] for o in range(-15,16)}
        values[f['id']]=data
    return values


def qx_encode_precip(caches,specs):
    for f in specs:
        rows=[]
        for branch in caches.values():
            data=branch[f['id']];rows.append(data['target']);rows.extend(r for yr in data['seasonal_years'].values() for r in yr.values())
        denominator=1
        for row in rows:
            for v in row.values():
                if v is not None:denominator=math.lcm(denominator,v.denominator)
        for row in rows:
            for h,v in row.items():
                if v is not None:row[h]=(v*denominator).numerator
        f['rank_scale']=denominator;f['exact_common_denominator']=denominator


def qx_c_precip(specs,values):
    rows=[]
    for f in specs:
        g=f['id'];data=values[g];target=data['target']['168'];annual={};pool=[]
        for year,records in data['years'].items():
            valid=[r['168'] for r in records.values() if r['168'] is not None]
            if len(valid)>=15 and target is not None:annual[year]={'n':len(valid),'cdf':float(a_midrank_fraction(sorted(valid),target))};pool.extend(valid)
        q=sum(r['cdf'] for r in annual.values())/len(annual) if len(annual)>=4 else None
        uninformative=bool(pool) and len(set(pool+[target]))==1
        rows.append({'id':g,'label':f['label'],'unit':f['unit'],'value':target/f['rank_scale'] if target is not None else None,'q':q,'A':abs(2*q-1) if q is not None else None,'annual':annual,'rank_eligible':not uninformative,'uninformative':uninformative,'reason':None if q is not None else 'target_or_15window_4year_support_NA','D1_support_rule':'complete precipitation window; no mean-gap relaxation'})
    return rows


def qx_pca_inputs(specs,values,primary):
    ids=[r['id'] for r in primary['candidate_ranking']];index={f['id']:f for f in specs};windows={}
    for g in ids:
        s=index[g]['rank_scale'];rec=values[g]
        convert=lambda row:{str(h):Fraction(row[str(h)],s) if row.get(str(h)) is not None else None for h in index[g]['window_hours']}
        windows[g]={'target':convert(rec['target']),'years':{str(y):{str(o):convert(rec['years'][str(y)][str(o)]) for o in range(-15,16)} for y in A_DEFAULT_YEARS}}
    global P3_FROZEN_V2
    old=P3_FROZEN_V2
    try:
        P3_FROZEN_V2={'feature_ids':ids,'features':[index[g] for g in ids]}
        inputs,rows=p3_inputs(windows)
    finally:P3_FROZEN_V2=old
    inputs['support']['rule']='All applicable windows of the full enhanced deduplicated candidate pool; same common endpoints at 7d/30d; no variable removed to salvage PCA'
    inputs['support']['note']='Support extension of approved V2 rule to D0+D1; no target-dependent endpoint selection'
    return inputs,rows


def qx_selftest():
    end=dt.datetime(2026,8,26,10,tzinfo=Q1_TZ)
    series={end-dt.timedelta(hours=i):Fraction(1) for i in range(72)}
    assert qx_precip_window(series,end,24,'precip_sum')==24
    assert qx_precip_window(series,end,24,'precip_max')==1
    assert qx_precip_window(series,end,24,'rain_hours')==24
    assert qx_precip_window(series,end,72,'api')==Fraction(24)*sum(Fraction(9,10)**i for i in range(3))

    series[end+dt.timedelta(hours=1)]=100000
    assert qx_precip_window(series,end,24,'precip_sum')==24
    missing=dict(series);del missing[end-dt.timedelta(hours=2)]
    assert qx_precip_window(missing,end,24,'precip_sum') is None
    zero={end-dt.timedelta(hours=i):Fraction(0) for i in range(24)}
    assert qx_precip_window(zero,end,24,'rain_hours')==0

    counter=[Fraction(0),Fraction(1),Fraction(3)];increments=[counter[i]-counter[i-1] for i in range(1,len(counter))]
    assert sum(increments)==3 and sum(counter)==4
    exact={'x':{g:{'target':{'24':Fraction(1,10)},'seasonal_years':{'2021':{'08-26':{'24':Fraction(3,10)}}}} for g in ['g']}}
    spec=[{'id':'g','rank_scale':1}];qx_encode_precip(exact,spec)
    assert Fraction(exact['x']['g']['target']['24'],spec[0]['rank_scale'])==Fraction(1,10)
    return {'status':'pass','checks':['hour_end_inclusion_and_next_hour_exclusion','complete_hour_support','sum_max_threshold_count','zero_is_valid_not_missing','counter_vs_increment_distinction','fixed_API_recursion_initial_state_and_decay','lossless_exact_encoding'],'scope':'author bounded checks, not independent approval'}


def qx_execute(args):
    started=dt.datetime.now(dt.timezone.utc).isoformat();clock0=time.monotonic();out=Path(args.output_dir);res=out/'results';res.mkdir(parents=True,exist_ok=True)
    if not args.d0_cache:
        if not args.data_dir:raise ValueError('enhance requires --data-dir to reconstruct D0 windows, or explicit --d0-cache')
        rebuilt=out/'reconstructed_d0'
        qx_rebuild_windows(argparse.Namespace(data_dir=args.data_dir,output_dir=str(rebuilt)))
        args.d0_cache=str(rebuilt/'results/window_cache.json.gz')
    source_files=[Path(__file__), Path(__file__).resolve().parent / 'data/features.json', Path(__file__).resolve().parent / 'data/reference.json',Path(args.d1_csv),Path(args.d0_cache),Path(args.baseline_results)/'analysis_results.json.gz',Path(args.baseline_results)/'c_fixed7d.csv']
    before=[q1_fp(p) for p in source_files];baseline=qf_read_json(Path(args.baseline_results)/'analysis_results.json.gz')
    features=copy.deepcopy(baseline['features']);new=qx_precip_specs(features);features+=new;active=[f for f in features if f.get('computable')]
    print('D1_READ_AND_EXACT_WINDOWS',flush=True)
    series,d1meta=qx_read_precip(args.d1_csv)
    write_d1_process(args.d1_csv,res)
    d1caches={hour:qx_precip_cache(series,new,hour) for hour in (10,0,9)};qx_encode_precip(d1caches,new)
    d0caches=qf_read_json(args.d0_cache)
    caches={key:{**data,**d1caches[int(key.split('_C')[1].split('_')[0])]} for key,data in d0caches.items()}
    values=caches['H-END_C10_0.9']
    print('ENHANCED_A_PRIMARY_AND_ALL_SENSITIVITIES',flush=True)
    primary=summarize_one(active,values)
    branches=[]
    for ep,ch,cov in [('H-END',10,.9),('H-START',10,.9),('H-END',0,.9),('H-START',0,.9),('H-END',10,.95)]:
        key=f'{ep}_C{ch}_{cov}';answer=primary if key=='H-END_C10_0.9' else summarize_one(active,caches[key])
        branches.append({'name':f'{ep}_C{ch:02d}_coverage{int(cov*100)}','endpoint':ep,'cutoff_hour':ch,'coverage':cov,'min_endpoints':25,'result':q1_trim_a(answer),'D1_time_policy':'actual product support; D0 endpoint assumption does not shift D1'})
    reduced=copy.deepcopy(active)
    for f in reduced:f['window_hours']=[h for h in f['window_hours'] if h!=720]
    branches.append({'name':'drop30d_all_factors','min_endpoints':25,'result':q1_trim_a(summarize_one(reduced,values))})
    for n in (20,15):branches.append({'name':f'all_factors_min{n}','min_endpoints':n,'result':q1_trim_a(summarize_one(active,values,min_endpoints=n))})
    for rho in (.85,.95):branches.append({'name':f'rho{rho}','min_endpoints':25,'result':q1_trim_a(summarize_one(active,values,redundancy_threshold=rho))})
    common,common_info=q1_common_support(active,values,primary)
    branches.append({'name':'common_history_support','min_endpoints':25,'common_support':common_info,'result':q1_trim_a(summarize_one(active,common))})
    print('ENHANCED_LEAVE_ONE_YEAR_OUT',flush=True);loo=leave_one_year_out(active,values)
    for k in loo['runs']:loo['runs'][k]=q1_trim_a(loo['runs'][k])
    singles=copy.deepcopy(new)
    for f in singles:f['direction']='high'
    single=summarize_one(singles,values)
    c_rows=copy.deepcopy(baseline['methods']['C_factors'])+qx_c_precip(new,values);c_index={r['id']:r for r in c_rows}
    specs={f['id']:f for f in features}
    common_ids=[r['id'] for r in primary['candidate_ranking'] if r['id'] in MECHANISM and c_index.get(r['id'],{}).get('A') is not None and c_index[r['id']].get('rank_eligible')]
    print('ENHANCED_ACTUAL_B_REFIT',len(common_ids),flush=True)
    b=v2b_fit(specs,values,common_ids,20261024,'B_D0_D1_common_actual_refit','enhanced_common_deduplicated_original_factor_analogues')
    b['metadata'].update(proxy_policy='Original mechanism grades inherited for separately named four ERA5 precipitation analogues; two height-difference proxies excluded because no original grades',D0_alignment='Unchanged exact D0 window cache plus separate product-defined D1 total-precipitation windows; no filling, clipping, recleaning or cross-source substitution',execution_in_enhancement=True,execution_in_v2=False,D1_role='Four original precipitation factors represented by a separately named fixed ERA5 gridded proxy; original grades transparently inherited',D0_reused_without_change=True)
    comparison=v2b_common_table(common_ids,specs,baseline['primary'],primary,b,c_index,scope='D0+D1 enhanced deduplicated shared original-factor analogues; no invented grades for height-difference proxies')
    for row in comparison:
        row.update(E1_scope='Unchanged D0 baseline ranking displayed only for old IDs; added D1 analogues have no baseline score',E2_scope='Enhanced D0+D1 A scores reranked in this actual common family',B_scope='Actual common-family entropy/TOPSIS refit with original prior mechanism grades inherited for D1 analogues',C_scope='Fixed7d baseline D0 plus complete-window D1 analogues, common-family reranking')
    print('ENHANCED_PCA_FIXED_7D_30D',flush=True)
    inputs,support_rows=qx_pca_inputs(active,values,primary)
    pca=build_pca_v3(inputs['feature_ids'],inputs['matrices_by_year'],inputs['targets_by_scale'],feature_meta=[{'id':g,'label':specs[g]['label'],'unit':specs[g]['unit']} for g in inputs['feature_ids']])
    pca.update(support=inputs['support'],metadata={'version':'D0+D1 complete enhanced pool','cutoff':'2026-08-26T10:00:00+08:00','D0_time_assumption':'H-END hourly mean working assumption','D1_time_definition':'documented previous hour total precipitation','independent_2026_evaluation':False,'scientific_acceptance':'pending independent review'})
    for scale in pca['scales'].values():
        for model in scale['models']:model['training_endpoint_ids_by_year']={y:inputs['support']['included_endpoint_ids_by_year'][y] for y in model['reference_years']}
        for conf in scale['configurations']:conf['calibration_endpoint_ids']=inputs['support']['included_endpoint_ids_by_year'][conf['calibration_year']]
    old_ids=[r['id'] for r in baseline['primary']['candidate_ranking']];old_ranks={r['id']:r for r in baseline['primary']['candidate_ranking']};new_ranks={r['id']:r for r in primary['candidate_ranking']}
    score_checks={g:baseline['primary']['factors'][g]['score_S_exact']==primary['factors'][g]['score_S_exact'] for g in old_ids}
    if not all(score_checks.values()):raise ValueError('Unchanged D0 marginal A scores changed')
    old_only_ranks=q1_rank_rows({g:primary['factors'][g]['score_S'] for g in old_ids});old_only={r['id']:r for r in old_only_ranks}
    changes=[{'id':g,'label':specs[g]['label'],'old_rank_min':old_ranks[g]['rank_min'],'enhanced_rank_min':new_ranks[g]['rank_min'],'common_old_only_rank_min':old_only[g]['rank_min'],'S_old':baseline['primary']['factors'][g]['score_S'],'S_enhanced':primary['factors'][g]['score_S'],'score_exact_unchanged':score_checks[g],'change_reason':'candidate addition only; marginal definition and own historical support unchanged'} for g in old_ids]
    result={'schema':'Q1_D0_D1_RECOMPUTED_20261030','features':features,'primary':primary,'branches':branches,'loo':loo,'single_sided_precipitation_hypotheses':single,'methods':{'B_common_actual_refit':b,'C_factors':c_rows,'B_D0_unchanged_reference':baseline['methods']['B_common_deduplicated_refit']},'comparisons':comparison,'pca':pca,'old_candidate_changes':changes,'D1_metadata':d1meta,'common_support':common_info,'D0_marginal_score_regression':score_checks,'scientific_acceptance':'pending independent Q1 total review','pre_score_operation_contract':'D1_OPERATION_CONTRACT.md'}
    qf_write_gzip(res/'enhanced_results.json.gz',result);qf_write_gzip(res/'enhanced_pca_results.json.gz',pca)
    qf_write_gzip(res/'enhanced_window_cache.json.gz',caches)
    q1_dump(res/'enhanced_pca_inputs.json',inputs);q1_csv(res/'enhanced_pca_common_support.csv',support_rows)
    pcadir=res/'enhanced_pca';pcadir.mkdir(exist_ok=True);p3_tables(pcadir,pca)
    q1_csv(res/'enhanced_method_comparison.csv',comparison);q1_csv(res/'enhanced_old_candidate_changes.csv',changes)
    q1_csv(res/'enhanced_sensitivity.csv',[{'branch':r['name'],'pool':r['result']['available_candidate_count'],'top':r['result']['top_candidates'],'boundary':r['result']['boundary_candidates']} for r in branches])
    q1_csv(res/'enhanced_loo_frequency.csv',[{'id':g,**r} for g,r in loo['selection_frequency'].items()])
    rank={r['id']:r for r in primary['candidate_ranking']};factor_rows=[];window_rows=[]
    for f in features:
        g=f['id'];r=primary['factors'].get(g,{});factor_rows.append({'id':g,'label':f['label'],'unit':f['unit'],'definition_status':f['status'],'score_status':r.get('status','definition_NA'),'S':r.get('score_S'),'rank_min':rank.get(g,{}).get('rank_min'),'rank_max':rank.get(g,{}).get('rank_max'),'top_candidate':g in primary['top_candidates'],'reason':r.get('reason') or f.get('status_reason'),'support':r.get('support')})
        for h,w in r.get('windows',{}).items():window_rows.append({'id':g,'label':f['label'],'hours':int(h),'value':w['target_value']/f['rank_scale'] if w['target_value'] is not None else None,'unit':f['unit'],**{k:v for k,v in w.items() if k!='target_value'}})
    q1_csv(res/'enhanced_all_factors.csv',factor_rows);q1_csv(res/'enhanced_primary_windows.csv',window_rows)
    checks={'bounded_precipitation_checks':qx_selftest(),'pca_reference_separation_and_contribution_sum':pca_v3_selftest(),'D0_marginal_scores_unchanged':score_checks,'target_history_same_product_and_grid':d1meta,'all_new_factor_year_endpoint_counts':{g:primary['factors'][g]['support'] for g in QX_NEW_IDS},'formal_new_B_counts':b['metadata'],'D1_no_QC_flags':'not_provided; never represented as QC0','API_rule_fixed_before_scoring':True,'direction':'main two-sided; inherited high-only hypothesis separate'}
    q1_dump(res/'ENHANCED_CHECKS.json',checks)
    after=[q1_fp(p) for p in source_files]
    if before!=after:raise ValueError('Input/source changed during enhanced execution')
    outputs=[p for p in sorted(res.rglob('*')) if p.is_file() and p.name not in ('RUN_RECEIPT.json','ENHANCED_RUN_RECEIPT.json') and 'baseline' not in p.parts]
    q1_dump(res/'ENHANCED_RUN_RECEIPT.json',{'record_type':'actual_D0_D1_enhancement','started_utc':started,'finished_utc':dt.datetime.now(dt.timezone.utc).isoformat(),'elapsed_seconds':time.monotonic()-clock0,'argv':sys.argv,'environment':{'python':sys.version,'numpy':np.__version__,'scipy':__import__('scipy').__version__},'inputs_before':before,'inputs_after':after,'inputs_unchanged':True,'outputs':[q1_fp(p) for p in outputs],'reused':'unchanged D0 windows and old separate B reference; baseline independently rerun by root with exact frozen cache checks','newly_computed':'D1 exact windows; full enhanced A/sensitivity/LOO; common B formal draws and bootstrap; C precipitation; full enhanced PCA','scientific_acceptance':'pending Q1 total review'})
    print('ENHANCED_COMPLETE',json.dumps({'top':primary['top_candidates'],'pool':len(primary['candidate_ranking']),'B_common':common_ids,'PCA_support':inputs['support']['n_by_year']},ensure_ascii=False),flush=True)


import argparse
import csv
import datetime as dt
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from pathlib import Path

D1_YEARS = (2021, 2022, 2023, 2024, 2025)
D1_PROBS = (("q05", Fraction(1,20)), ("q25", Fraction(1,4)),
            ("median", Fraction(1,2)), ("q75", Fraction(3,4)), ("q95", Fraction(19,20)))
D1_ZONE = dt.timezone(dt.timedelta(hours=8))


def d1_plain(obj):
    if isinstance(obj, dict):
        return {str(k):d1_plain(v) for k,v in obj.items()}
    if isinstance(obj, (tuple,list)):
        return [d1_plain(v) for v in obj]
    if isinstance(obj, (Fraction,Decimal)):
        return float(obj)
    return obj


def d1_equal_year_quantiles(groups):

    if len(groups) != 5 or any(not g for g in groups):
        return {name:None for name,_ in D1_PROBS}
    weighted = sorted((v,Fraction(1,5*len(g))) for g in groups for v in g)
    result = {}
    for name, p in D1_PROBS:
        weight = Fraction(0)
        for value, w in weighted:
            weight += w
            if weight >= p:
                result[name] = value
                break
    return result


def d1_load_precipitation(csv_path):
    records = {}
    grid = set()
    with Path(csv_path).open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            end = dt.datetime.fromisoformat(row["interval_end"])
            start = dt.datetime.fromisoformat(row["interval_start"])
            assert end.utcoffset() == dt.timedelta(hours=8)
            assert end.minute == end.second == end.microsecond == 0
            assert start == end-dt.timedelta(hours=1)
            assert end not in records, f"duplicate support end: {end}"
            text = row["precipitation_mm"].strip()
            value = Fraction(Decimal(text)) if text else None
            assert value is None or value >= 0
            assert row["source_product"] == "Open-Meteo archive / ERA5"
            assert row["source_qc_status"] == "not_provided"
            records[end] = value
            grid.add((float(row["grid_latitude"]),float(row["grid_longitude"]),float(row["grid_elevation_m"])))
    assert len(grid) == 1
    return records, next(iter(grid))


def d1_episodes(values, thresholds, ends, relative):
    flags = [None if v is None or q is None else v>q for v,q in zip(values,thresholds)]
    episodes, active = [], []
    def finish():
        if not active:
            return
        first,last = active[0],active[-1]
        episodes.append({"support_start":(ends[first]-dt.timedelta(hours=1)).isoformat(),
                         "support_end":ends[last].isoformat(),
                         "first_relative_hour":relative[first],"last_relative_hour":relative[last],
                         "duration_hours":len(active),
                         "sum_precipitation_mm":sum(values[k] for k in active),
                         "sum_excess_mm":sum(values[k]-thresholds[k] for k in active),
                         "max_hourly_precipitation_mm":max(values[k] for k in active),
                         "max_excess_mm":max(values[k]-thresholds[k] for k in active)})
        active.clear()
    for k,flag in enumerate(flags):
        if flag is True:
            active.append(k)
        else:
            finish()
    finish()
    return flags, episodes


def d1_period_summary(values, thresholds, ends, relative):
    available = [v for v in values if v is not None]
    flags,episodes = d1_episodes(values,thresholds,ends,relative)
    complete = len(available)==len(values)
    summary = {"expected_hours":len(values),"available_hours":len(available),
               "missing_hours":len(values)-len(available),"complete_support":complete,
               "support_start":(ends[0]-dt.timedelta(hours=1)).isoformat(),
               "support_end":ends[-1].isoformat(),
               "sum_precipitation_mm":sum(available) if complete else None,
               "available_sum_precipitation_mm":sum(available) if available else None,
               "max_hourly_precipitation_mm":max(available) if available else None,
               "wet_hours_ge_0_1mm":sum(v>=Fraction(1,10) for v in available),
               "zero_hours":sum(v==0 for v in available),
               "threshold_available_hours":sum(q is not None for q in thresholds),
               "assessed_hours":sum(f is not None for f in flags),
               "unassessed_hours":sum(f is None for f in flags),
               "high_p95_hours":sum(f is True for f in flags),
               "high_p95_episode_count":len(episodes),
               "longest_high_p95_episode_hours":max((e["duration_hours"] for e in episodes),default=0),
               "episodes":episodes,
               "count_note":"wet/zero/exceedance counts refer only to available/assessed hours; missing is never zero"}
    return summary


def build_d1_process(csv_path, cutoff=None):

    cutoff = cutoff or dt.datetime(2026,8,26,10,tzinfo=D1_ZONE)
    assert cutoff.utcoffset()==dt.timedelta(hours=8)
    assert cutoff.minute==cutoff.second==cutoff.microsecond==0
    records,grid = d1_load_precipitation(csv_path)
    relative = list(range(-719,1))
    ends = [cutoff+dt.timedelta(hours=h) for h in relative]
    values = [records.get(t) for t in ends]
    quantiles = {key:[] for key,_ in D1_PROBS}
    counts = [[] for _ in D1_YEARS]
    for end in ends:
        groups = []
        for yi,year in enumerate(D1_YEARS):
            group = []
            for day in range(-15,16):
                reference = end.replace(year=year)+dt.timedelta(days=day)
                assert reference.year in D1_YEARS
                v = records.get(reference)
                if v is not None:
                    group.append(v)
            counts[yi].append(len(group))
            groups.append(group)
        q = d1_equal_year_quantiles(groups)
        for key in quantiles:
            quantiles[key].append(q[key])
    thresholds = quantiles["q95"]
    flags,episodes = d1_episodes(values,thresholds,ends,relative)
    summaries = {name:d1_period_summary(values[-hours:],thresholds[-hours:],ends[-hours:],relative[-hours:])
                 for name,hours in (("30d",720),("7d",168),("72h",72))}
    assert len(ends)==720 and ends[0]-dt.timedelta(hours=1)==cutoff-dt.timedelta(days=30)
    assert ends[-1]==cutoff and all(t<=cutoff for t in ends)
    assert len(set(ends))==720
    assert all(all(n==31 for n in year) for year in counts), "incomplete D1 historical support"
    assert all(v is not None for v in values), "incomplete D1 target support"
    result = {"metadata":{"schema_version":"q1-d1-process-v1.0", "primary_cutoff":cutoff.isoformat(),
                          "timezone":"UTC+8", "source":"Open-Meteo archive API / fixed ERA5",
                          "requested_point":{"latitude":28.86,"longitude":85.29},
                          "returned_grid":{"latitude":grid[0],"longitude":grid[1],"elevation_m":grid[2]},
                          "cell_selection":"nearest","elevation_request":"nan; downscaling disabled",
                          "support":"API precipitation at t is total over preceding hour [t-1h,t)",
                          "unit":"mm","source_qc_status":"not_provided",
                          "role":"D1 independent grid proxy, not D0 station or source-slope observation",
                          "source_csv":str(csv_path),
                          "source_sha256":hashlib.sha256(Path(csv_path).read_bytes()).hexdigest(),
                          "reference_years":list(D1_YEARS),"seasonal_days_each_side":15,
                          "quantile_rule":"Per relative hour: each of five years has weight 1/5; each eligible sample in year y has weight 1/(5*n_y); inverse ECDF inf{x:F(x)>=p}, no interpolation; all five years required.",
                          "threshold_rule":"strictly greater than historical q95; historical years only, never target values",
                          "trajectory_note":"Historical same-calendar ±15d trajectories overlap; not independent events or a confidence interval.",
                          "future_target_support_used":False,
                          "limitations":"Reanalysis grid precipitation, terrain-smoothed spatial support; no causal attribution or no-rain-trigger inference."},
              "time_relative_hours":relative,
              "support_start_iso":[(t-dt.timedelta(hours=1)).isoformat() for t in ends],
              "support_end_iso":[t.isoformat() for t in ends],
              "values":values,
              **{f"background_{key}":q for key,q in quantiles.items()},
              "history_n_by_year":counts,
              "history_years_count":[sum(c[k]>0 for c in counts) for k in range(720)],
              "high_p95_exceedance":flags,
              "high_p95_signed_excess_mm":[v-q if v is not None and q is not None else None for v,q in zip(values,thresholds)],
              "high_p95_episodes":episodes,"summaries":summaries}
    return d1_plain(result)


def d1_process_selftest():
    q=d1_equal_year_quantiles([[Fraction(0)],[Fraction(10)]*100,[Fraction(20)],[Fraction(30)],[Fraction(40)]])
    assert q["median"]==20, "annual weighting accidentally pooled"
    assert d1_equal_year_quantiles([[],[1],[1],[1],[1]])["q95"] is None
    ends=[dt.datetime(2026,8,26,7+i,tzinfo=D1_ZONE) for i in range(4)]
    flags,episodes=d1_episodes([Fraction(1),Fraction(2),None,Fraction(3)],[Fraction(1)]*4,ends,[-3,-2,-1,0])
    assert flags==[False,True,None,True] and len(episodes)==2
    s=d1_period_summary([Fraction(0),None],[Fraction(0),Fraction(0)],ends[:2],[-1,0])
    assert s["sum_precipitation_mm"] is None and s["available_sum_precipitation_mm"]==0
    assert s["zero_hours"]==1 and s["missing_hours"]==1 and s["high_p95_hours"]==0
    return {"unequal_annual_n_equal_weights":"pass","one_empty_year_quantiles_na":"pass",
            "strict_threshold_and_missing_breaks_episodes":"pass","zero_distinct_from_missing":"pass"}


def write_d1_process(csv_path, output_dir):
    output_dir=Path(output_dir); output_dir.mkdir(parents=True,exist_ok=True)
    result=build_d1_process(csv_path)
    (output_dir/"d1_process.json").write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False),encoding="utf-8")
    fields=["relative_hour","support_start","support_end","precipitation_mm","background_q05",
            "background_q25","background_median","background_q75","background_q95",
            *[f"history_n_{year}" for year in D1_YEARS],"history_years_count","high_p95_exceedance",
            "high_p95_signed_excess_mm","source_product","source_qc_status"]
    with (output_dir/"process_source_data.csv").open("w",encoding="utf-8",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader()
        for k,rel in enumerate(result["time_relative_hours"]):
            writer.writerow({"relative_hour":rel,"support_start":result["support_start_iso"][k],
                             "support_end":result["support_end_iso"][k],"precipitation_mm":result["values"][k],
                             **{f"background_{key}":result[f"background_{key}"][k] for key,_ in D1_PROBS},
                             **{f"history_n_{year}":result["history_n_by_year"][yi][k] for yi,year in enumerate(D1_YEARS)},
                             "history_years_count":result["history_years_count"][k],
                             "high_p95_exceedance":result["high_p95_exceedance"][k],
                             "high_p95_signed_excess_mm":result["high_p95_signed_excess_mm"][k],
                             "source_product":"Open-Meteo archive / ERA5","source_qc_status":"not_provided"})
    (output_dir/"d1_process_checks.json").write_text(json.dumps({"selftest":d1_process_selftest(),
                    "actual_target_hours":len(result["values"]),"actual_target_missing_hours":sum(v is None for v in result["values"]),
                    "historical_count_per_year_per_relative_hour":31,"historical_years":5,
                    "actual_cutoff_and_source_intervals":"pass","not_independent_review":True},ensure_ascii=False,indent=2),encoding="utf-8")
    return result


def qx_rebuild_windows(args):

    started=dt.datetime.now(dt.timezone.utc).isoformat();out=Path(args.output_dir)/'results';out.mkdir(parents=True,exist_ok=True)
    base=Path(args.data_dir);files=[Path(__file__), Path(__file__).resolve().parent / 'data/features.json', Path(__file__).resolve().parent / 'data/reference.json',base/'run_summary_all.json',base/'attachment01_02/variable_dictionary.csv']+[base/'attachment01_02/station'/(t+'.csv.gz') for t in ('MET','SOIL','RADM','FLUX')]
    before=[q1_fp(p) for p in files];dictionary=q1_features();features=dictionary['factors'];active=[f for f in features if f.get('computable')]
    arrays,counts=q1_read_d0(base,features)
    for f in active:
        _,_,places,_=q2_inputs(arrays,f);f['source_grid_scale']=10**places;f['rank_scale']=10**places
    caches={};ledger=[]
    for ep,ch,cov in [('H-END',10,.9),('H-START',10,.9),('H-END',0,.9),('H-START',0,.9),('H-END',10,.95),('H-END',9,.9)]:
        print('REBUILD_D0_EXACT_WINDOWS',ep,ch,cov,flush=True);caches[(ep,ch,cov)]=q1_cache(arrays,active,ep,ch,cov,ledger)
        ledger.clear()
    q1_encode_exact_cache(caches,features)
    qf_write_gzip(out/'window_cache.json.gz',{f'{ep}_C{ch}_{cov}':data for (ep,ch,cov),data in caches.items()})
    q1_dump(out/'feature_dictionary.json',dictionary)
    after=[q1_fp(p) for p in files];assert before==after
    q1_dump(out/'WINDOW_REBUILD_RECEIPT.json',{'record_type':'actual_D0_exact_windows_reconstruction','started_utc':started,'finished_utc':dt.datetime.now(dt.timezone.utc).isoformat(),'inputs_before':before,'inputs_after':after,'inputs_unchanged':True,'outputs':[q1_fp(out/'window_cache.json.gz'),q1_fp(out/'feature_dictionary.json')],'data_counts':counts,'analysis_or_B_sampling_executed':False,'scope':'D0 read/exact aggregation only; does not approve scientific results'})
    print('D0_WINDOWS_REBUILT',str(out/'window_cache.json.gz'),flush=True)


def sol():
    mode=sys.argv[sys.argv.index('--mode')+1] if '--mode' in sys.argv else 'verify'
    if mode not in ('enhance','rebuild-windows'):return qf_main()
    parser=argparse.ArgumentParser(description='Explicit independent ERA5 D1 precipitation enhancement; approved Q1 mathematics retained')
    parser.add_argument('--mode',choices=['enhance','rebuild-windows'],required=True);parser.add_argument('--d1-csv');parser.add_argument('--d0-cache');parser.add_argument('--data-dir');parser.add_argument('--baseline-results');parser.add_argument('--output-dir',required=True)
    args=parser.parse_args()
    if mode=='rebuild-windows':
        if not args.data_dir:parser.error('rebuild-windows requires --data-dir')
        qx_rebuild_windows(args)
    else:
        if not args.d1_csv or not args.baseline_results:parser.error('enhance requires --d1-csv and --baseline-results')
        if not args.d0_cache and not args.data_dir:parser.error('enhance requires --data-dir or --d0-cache')
        qx_execute(args)


if __name__=='__main__':sol()
