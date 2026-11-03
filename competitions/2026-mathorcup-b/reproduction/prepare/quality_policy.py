from dataclasses import dataclass
from datetime import datetime, timezone
from math import isfinite, isnan
from numbers import Real


@dataclass(frozen=True)
class PhysicalBounds:
    lower: float = float("-inf")
    upper: float = float("inf")
    lower_inclusive: bool = True
    upper_inclusive: bool = True

    def __post_init__(self):
        if isnan(self.lower) or isnan(self.upper) or self.lower > self.upper:
            raise ValueError("Invalid physical bounds")
        if self.lower == self.upper and not (self.lower_inclusive and self.upper_inclusive):
            raise ValueError("Empty physical domain")

    def accepts(self, value):
        return ((value >= self.lower if self.lower_inclusive else value > self.lower)
                and (value <= self.upper if self.upper_inclusive else value < self.upper))


@dataclass(frozen=True)
class QualityPolicy:
    no_source_qc: str = "reject"
    retain_flagged_codes: tuple = ()
    sentinel_tokens: tuple = (9999.9,)

    def __post_init__(self):

        object.__setattr__(self, "retain_flagged_codes", tuple(self.retain_flagged_codes))
        object.__setattr__(self, "sentinel_tokens", tuple(self.sentinel_tokens))
        if self.no_source_qc not in {"reject", "retain_flagged"}:
            raise ValueError("Unknown no-source-QC policy")
        if any(type(c) is not int or c not in {1, 4, 9} for c in self.retain_flagged_codes):
            raise ValueError("Only QC1/4/9 may be explicitly retained")
        if any(isinstance(x, bool) or not isinstance(x, Real) or not isfinite(x)
               for x in self.sentinel_tokens):
            raise ValueError("Sentinel tokens must be explicit finite numbers")


def _aware(value):
    return isinstance(value, datetime) and value.tzinfo is not None and value.utcoffset() is not None


def assess_cell(raw_value, *, source_qc_presence, source_qc=None,
                physical_bounds=None, policy=QualityPolicy(),
                interval_start=None, interval_end=None, aggregation=None,
                available_at=None, mode="retrospective", as_of=None,
                empirical_flags=()):


    if source_qc_presence not in {"not_provided", "provided"}:
        raise ValueError("source_qc_presence must be explicit")
    if source_qc_presence == "not_provided" and source_qc is not None:
        raise ValueError("Do not fabricate a code for a column without source QC")
    if mode not in {"retrospective", "causal"}:
        raise ValueError("Unknown mode")
    if mode == "causal" and not _aware(as_of):
        raise ValueError("Causal mode requires timezone-aware as_of")
    if mode == "retrospective" and as_of is not None:
        raise ValueError("as_of applies only in causal mode")
    if isinstance(empirical_flags, str):
        raise ValueError("empirical_flags must be a collection of nonempty strings")
    flags = tuple(empirical_flags)
    if any(not isinstance(x, str) or not x for x in flags):
        raise ValueError("empirical_flags must be a collection of nonempty strings")
    reasons = []
    numeric = None
    if isinstance(raw_value, bool):
        reasons.append("nonnumeric")
    else:
        try:
            numeric = float(raw_value)
        except (TypeError, ValueError, OverflowError):
            reasons.append("nonnumeric")
    structural = numeric is not None and isfinite(numeric)
    if numeric is not None and not isfinite(numeric):
        reasons.append("nonfinite")
    if structural and numeric in policy.sentinel_tokens:
        structural = False
        reasons.append("declared_sentinel")

    qc_known = (isinstance(source_qc, Real) and not isinstance(source_qc, bool)
                and isfinite(source_qc) and source_qc in {0, 1, 2, 4, 8, 9})
    if source_qc_presence == "not_provided":
        source_status = "source_quality_unverified"
        source_allowed = policy.no_source_qc == "retain_flagged"
    elif source_qc is None:
        source_status, source_allowed = "provided_code_missing", False
    elif not qc_known:
        source_status, source_allowed = "provided_code_unknown", False
    else:
        source_status = {0: "source_flag_0_correct", 1: "source_flag_1_suspicious",
                         2: "source_flag_2_wrong", 4: "source_flag_4_imputed",
                         8: "source_flag_8_missing", 9: "source_flag_9_uncontrolled"}[source_qc]
        source_allowed = (source_qc not in {2, 8}
                          and (source_qc == 0 or source_qc in policy.retain_flagged_codes))
    if not source_allowed:
        reasons.append(source_status)

    if physical_bounds is not None and not isinstance(physical_bounds, PhysicalBounds):
        raise TypeError("physical_bounds must be an explicit PhysicalBounds")
    physical = ("not_evaluated" if not structural or physical_bounds is None else
                "pass" if physical_bounds.accepts(numeric) else "fail")
    if physical != "pass":
        reasons.append("physical_" + physical)

    aware_support = _aware(interval_start) and _aware(interval_end)

    start_utc = interval_start.astimezone(timezone.utc) if aware_support else None
    end_utc = interval_end.astimezone(timezone.utc) if aware_support else None
    time_known = bool(aware_support and (
        (aggregation == "instantaneous" and start_utc == end_utc)
        or (aggregation in {"mean", "sum"} and start_utc < end_utc)))
    time_status = "declared_support_known" if time_known else "unknown_or_invalid_support"
    if not time_known:
        reasons.append(time_status)
    release_known = _aware(available_at)
    release_utc = available_at.astimezone(timezone.utc) if release_known else None
    cutoff_utc = as_of.astimezone(timezone.utc) if mode == "causal" else None
    release_consistent = release_known and (not time_known or release_utc >= end_utc)
    availability_status = ("unknown" if not release_known else
                           "invalid_before_support_end" if not release_consistent else "known")
    causal_available = bool(time_known and release_consistent and mode == "causal"
                            and end_utc <= cutoff_utc and release_utc <= cutoff_utc)
    if mode == "causal" and not causal_available:
        reasons.append("causal_availability_not_established")
    if release_known and not release_consistent:
        reasons.append(availability_status)
    candidate = structural and physical == "pass" and source_allowed
    usable = bool(candidate and time_known and (not release_known or release_consistent)
                  and (mode == "retrospective" or causal_available))
    return dict(raw_value=raw_value, source_qc=source_qc,
                source_qc_presence=source_qc_presence, source_quality_status=source_status,
                source_policy_eligible=bool(source_allowed),
                structural_status="pass" if structural else "fail", physical_status=physical,
                time_support_status=time_status, interval_start=interval_start, interval_end=interval_end,
                aggregation=aggregation, availability_status=availability_status,
                available_at=available_at, causal_available=causal_available, mode=mode,
                research_candidate=bool(candidate), research_usable=usable,
                research_value=numeric if usable else None, empirical_flags=flags,
                empirical_status="suspect" if flags else "not_flagged",
                source_imputed=bool(qc_known and source_qc == 4), imputed_here=False,
                clipped_here=False, reasons=tuple(reasons))
