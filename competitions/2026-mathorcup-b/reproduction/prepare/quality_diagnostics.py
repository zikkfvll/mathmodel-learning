import json
import math
from numbers import Integral

import pandas as pd


def _scalar(value):

    if value is None or pd.isna(value):
        return None
    if hasattr(value, "item"):
        value = value.item()
    return value


def _label(frame, position):
    value = _scalar(frame.at[position, "timestamp_label"])
    return str(value) if value is not None else None


def _row(frame, position):
    return _scalar(frame.at[position, "source_row"])


def _json(details):


    def safe(value):
        if isinstance(value, dict):
            return {str(key): safe(item) for key, item in value.items()}
        if isinstance(value, (tuple, list)):
            return [safe(item) for item in value]
        value = _scalar(value)
        if isinstance(value, float) and not math.isfinite(value):
            return str(value)
        return value

    return json.dumps(safe(details), ensure_ascii=False, separators=(",", ":"), allow_nan=False)


def diagnose(data: pd.DataFrame, measurement_columns: list[str], *,
             dataset: str, step_hours: int) -> list[dict]:


    if isinstance(step_hours, bool) or not isinstance(step_hours, Integral) or step_hours <= 0:
        raise ValueError("step_hours must be a positive integer")
    if data.empty:
        return []
    columns = list(dict.fromkeys(measurement_columns))
    required = ["timestamp_label", "source_file", "source_sheet", "source_row", *columns]
    missing = [column for column in required if column not in data.columns]
    if missing:
        raise ValueError(f"Missing diagnostic input columns: {missing}")
    if not data.columns.is_unique:
        raise ValueError("Diagnostic input column names must be unique")

    output = []
    step = pd.Timedelta(hours=int(step_hours))
    kind = dataset.upper()
    for (source_file, source_sheet), source in data.groupby(
            ["source_file", "source_sheet"], sort=False, dropna=False):
        frame = source.reset_index(drop=True)


        times = pd.to_datetime(frame["timestamp_label"], format="ISO8601", errors="coerce", utc=True)
        numeric = frame[columns].apply(pd.to_numeric, errors="coerce")
        finite = numeric.notna() & numeric.abs().lt(float("inf"))
        sheet_detail = {"source_sheet": _scalar(source_sheet)}

        def emit(check, variable, count, start=None, end=None,
                 before=None, after=None, **details):
            output.append({
                "dataset": dataset,
                "source_file": _scalar(source_file),
                "check": check,
                "variable": variable,
                "count": int(count),
                "start": _label(frame, start) if start is not None else None,
                "end": _label(frame, end) if end is not None else None,
                "value_before": _scalar(before),
                "value_after": _scalar(after),
                "detail": _json({**sheet_detail, **details}),
            })

        for column in columns:
            usable = finite[column] & times.notna()
            observed = numeric.loc[usable, column]
            eligible = times.loc[usable].diff().eq(step)
            changes = observed.diff().loc[eligible]
            if changes.empty:
                emit("max_adjacent_change", column, 0,
                     reason="no_adjacent_valid_pair_at_declared_step", step_hours=int(step_hours))
            else:
                end = int(changes.abs().idxmax())
                end_offset = observed.index.get_loc(end)
                start = int(observed.index[end_offset - 1])
                emit("max_adjacent_change", column, len(changes), start, end,
                     float(observed.at[start]), float(observed.at[end]),
                     source_rows=[_row(frame, start), _row(frame, end)],
                     signed_change=float(changes.at[end]), absolute_change=float(abs(changes.at[end])),
                     step_hours=int(step_hours))


            continuation = (usable & usable.shift(1, fill_value=False)
                            & times.diff().eq(step) & numeric[column].eq(numeric[column].shift(1)))
            run_ids = (~continuation).cumsum()
            runs = numeric.loc[usable, column].groupby(run_ids.loc[usable], sort=False)
            sizes = runs.size()
            if sizes.empty:
                emit("longest_exact_constant_run", column, 0, reason="no_valid_timed_value")
            else:
                best_run = sizes.idxmax()
                positions = run_ids.index[usable & run_ids.eq(best_run)]
                start, end = int(positions[0]), int(positions[-1])
                value = float(numeric.at[start, column])
                emit("longest_exact_constant_run", column, int(sizes.at[best_run]),
                     start, end, value, value,
                     source_rows=[_row(frame, start), _row(frame, end)],
                     equality="exact_float_equality", step_hours=int(step_hours))

        def layer_spread(prefix, group_name):
            selected = [column for column in columns if column.lower().startswith(prefix)]
            if len(selected) < 2:
                return
            values = numeric[selected].where(finite[selected])
            simultaneous = values.notna().sum(axis=1).ge(2) & times.notna()
            spreads = (values.max(axis=1) - values.min(axis=1)).loc[simultaneous]
            if spreads.empty:
                emit("max_simultaneous_layer_spread", group_name, 0,
                     reason="no_time_with_two_valid_layers", total_layers=len(selected))
                return
            position = int(spreads.idxmax())
            valid_layers = values.loc[position].dropna()
            low, high = valid_layers.idxmin(), valid_layers.idxmax()
            emit("max_simultaneous_layer_spread", group_name, len(spreads),
                 position, position, float(valid_layers.at[low]), float(valid_layers.at[high]),
                 source_row=_row(frame, position), spread=float(spreads.at[position]),
                 low_layer=low, high_layer=high, valid_layers=len(valid_layers), total_layers=len(selected),
                 layer_values={column: float(value) for column, value in valid_layers.items()},
                 interpretation="observation_only_not_an_error_test")

        def range_observation(selected, name, predicate):
            if not selected:
                return
            mask = predicate(numeric[selected]) & finite[selected]
            cell_count = int(mask.sum().sum())
            details = {"count_unit": "observed_cells", "distinct_source_rows": int(mask.any(axis=1).sum()),
                       "counts_by_variable": {column: int(mask[column].sum()) for column in selected},
                       "interpretation": "observation_only_no_value_removed"}
            if not cell_count:
                emit(name, ";".join(selected), 0, **details)
                return
            position = int(mask.any(axis=1).loc[lambda series: series].index[0])
            column = next(column for column in selected if mask.at[position, column])
            emit(name, ";".join(selected), cell_count, position, position,
                 float(numeric.at[position, column]), source_row=_row(frame, position),
                 representative_variable=column, representative="first_in_source_order", **details)

        if kind == "MET":
            layer_spread("ta_", "Ta_layers")
            range_observation([c for c in columns if c.upper().startswith(("RH_", "RH("))],
                              "observed_RH_outside_0_100", lambda values: (values < 0) | (values > 100))
            range_observation([c for c in columns if c.upper().startswith(("WS_", "WS("))],
                              "observed_negative_WS", lambda values: values < 0)
            range_observation([c for c in columns if c.upper().startswith(("WD_", "WD("))],
                              "observed_WD_outside_0_360", lambda values: (values < 0) | (values > 360))
        elif kind == "SOIL":
            layer_spread("soilt_", "soilT_layers")
            layer_spread("vwc_", "VWC_layers")
        elif kind == "RADM":
            radiation = [c for c in columns if c.split("(", 1)[0].lower() in {"rsd", "rsu", "rld", "rlu"}]
            range_observation(radiation, "observed_negative_radiation", lambda values: values < 0)
    return output
