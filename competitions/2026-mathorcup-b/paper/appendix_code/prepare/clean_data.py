import argparse
from contextlib import ExitStack, contextmanager
from decimal import Decimal, InvalidOperation
import collections
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import math
import pathlib
import re
import sys
import xml.etree.ElementTree as ET
import zipfile
import pandas as pd
import quality_policy as API
from quality_diagnostics import diagnose

BASE = pathlib.Path(__file__).resolve().parent
RULES = json.loads((BASE / "quality_rules.json").read_text(encoding="utf-8"))
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def physical(variable, station=True):
    if not station: return "none_daily_unit_pending", None
    rules = RULES["physical_rules"]
    if variable.startswith(("Ta_", "soilT_")) or variable == "T_SONIC(deg C)": key = "temperature_C"
    elif variable in rules["pressure"]["columns"]: key = "pressure"
    elif variable.startswith("VWC_") and variable.endswith("(m3/m3)"): key = "vwc"
    else: return "none_source_definition_pending", None
    r = rules[key]
    return key, API.PhysicalBounds(lower=r.get("lower", -math.inf), upper=r.get("upper", math.inf),
        lower_inclusive=r.get("lower_inclusive", True), upper_inclusive=r.get("upper_inclusive", True))


NUMBER_GRAMMAR = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z", re.ASCII)
ASCII_SPACE = " \t\r\n\v\f"


def parse_number(raw, station=True):

    if not isinstance(raw, str):
        return None, "nonnumeric"
    token = raw.strip(ASCII_SPACE)
    if not token:
        return None, "empty_token"
    if re.fullmatch(r"[+-]?(?:inf(?:inity)?|s?nan[0-9]*)", token, re.IGNORECASE | re.ASCII):
        return None, "nonfinite"
    if not NUMBER_GRAMMAR.fullmatch(token):
        return None, "nonnumeric"
    try:
        exact = Decimal(token)
    except InvalidOperation:
        return None, "decimal_conversion_invalid"
    if not exact.is_finite():
        return None, "nonfinite"
    if station and exact == Decimal("9999.9"):
        return None, "declared_sentinel"
    try:
        value = float(exact)
    except (ValueError, OverflowError):
        return None, "numeric_conversion_nonfinite"
    if not math.isfinite(value):
        return None, "numeric_conversion_nonfinite"
    if exact != 0 and value == 0:
        return None, "numeric_conversion_underflow"
    return value, None


def assess(raw, variable, presence, qc_raw, label_status="unique", station=True):
    qc = None
    if qc_raw.strip():
        try: qc = float(qc_raw)
        except ValueError: qc = qc_raw
    rule, bounds = physical(variable, station)
    val, parse_reason = parse_number(raw, station)

    out = API.assess_cell(val if parse_reason is None else math.nan,
        source_qc_presence=presence, source_qc=qc, physical_bounds=bounds,
        policy=API.QualityPolicy(no_source_qc=RULES["no_source_qc"],
            retain_flagged_codes=tuple(c for c in RULES["main_allowed_qc_codes"] if c != 0),
            sentinel_tokens=()))
    out["raw_value"] = raw
    out["structural_status"] = "pass" if parse_reason is None else "fail"
    reasons = ([parse_reason] if parse_reason else []) + [r for r in out["reasons"] if r not in {"nonfinite","nonnumeric","declared_sentinel"}]
    flags = []
    if station and variable == "PA(kPa)":
        flags.append("pa_definition_reference_pending")
        reasons.append("pa_definition_reference_pending")
    if val is not None and parse_reason is None:
        if val < 0 and variable == "Precip(mm)": flags.append("precip_nonnegative_definition_pending")
        if val < 0 and variable == "USTAR(m/s)": flags.append("ustar_nonnegative_definition_pending")
        if val < 0 and variable in {"Rsd(W/m2)","Rsu(W/m2)","Rld(W/m2)","Rlu(W/m2)"}:
            flags.append("radiation_direction_noise_pending")
    if label_status != "unique": reasons.append("time_label_"+label_status)
    structural = parse_reason is None
    source_candidate = structural and out["source_policy_eligible"]
    candidate = source_candidate and out["physical_status"] != "fail" and label_status == "unique"
    screened = candidate and out["physical_status"] == "pass"
    if not structural: part = "quarantine_structure"
    elif not out["source_policy_eligible"]: part = "quarantine_source_qc"
    elif out["physical_status"] == "fail": part = "quarantine_physical"
    elif label_status != "unique": part = "hold_time_label_conflict"
    elif out["physical_status"] != "pass": part = "retain_physical_pending"
    else: part = "retain_screened_candidate"
    out.update(physical_rule=rule, conditional_flags=";".join(flags), primary_disposition=part, empirical_status="not_evaluated",
        reasons=";".join(dict.fromkeys(reasons)), time_label_status=label_status,
        source_candidate_value=val if source_candidate else None,
        candidate_value=val if candidate else None, screened_candidate_value=val if screened else None,
        research_candidate=bool(screened))
    return out


def parse_label(raw):
    try: return dt.datetime.strptime(raw.replace("/", "-"), "%Y-%m-%d %H:%M")
    except ValueError: return None


def label_statuses(labels, records, year, hourly):
    counts = collections.Counter(x for x in labels if x is not None)
    statuses = []
    prev = None
    for i, t in enumerate(labels):
        why = []
        if t is None: why.append("unparseable")
        else:
            if counts[t] > 1: why.append("duplicate")
            if t.year != year: why.append("year_mismatch")
            if t.minute or t.second or t.microsecond or (not hourly and t.hour): why.append("off_grid")
            if prev and t < prev: why.append("out_of_order")
            prev = t
        statuses.append("+".join(why) or "unique")
    conflicts = []
    for t, n in counts.items():
        if n > 1:
            ix = [i for i, x in enumerate(labels) if x == t]
            conflicts.append({"timestamp_label":t.isoformat(),"rows":[i+2 for i in ix],"count":n,
                              "identical_full_records":len({tuple(map(str,records[i])) for i in ix}) == 1})
    return statuses, conflicts


def excel_rows(path):

    import openpyxl
    decoded = openpyxl.load_workbook(path, read_only=True, data_only=False)
    try:
        with zipfile.ZipFile(path) as z:
            strings = []
            if "xl/sharedStrings.xml" in z.namelist():
                strings = ["".join(n.itertext()) for n in ET.fromstring(z.read("xl/sharedStrings.xml"))]

            relns = "{http://schemas.openxmlformats.org/package/2006/relationships}"
            rels = {a.attrib["Id"]: a.attrib["Target"] for a in ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))}
            names = ET.fromstring(z.read("xl/workbook.xml")).find("m:sheets", NS)
            for s in names:
                rid = s.attrib["{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"]
                target = rels[rid]
                target = target.lstrip("/") if target.startswith("/") else "xl/"+target
                sheet = s.attrib["name"]
                xml = ET.fromstring(z.read(target))
                raw_rows = []
                for row in xml.findall("m:sheetData/m:row", NS):
                    cells = {}
                    for c in row:
                        ref = c.attrib["r"]
                        typ = c.attrib.get("t", "n")
                        v = c.find("m:v", NS)
                        token = v.text if v is not None and v.text is not None else ""
                        if typ == "s": value = strings[int(token)]
                        elif typ == "inlineStr": value = "".join(c.find("m:is", NS).itertext())
                        else: value = token
                        formula = c.find("m:f", NS)
                        cells[re.sub(r"[0-9]", "", ref)] = (value,typ,formula.text if formula is not None else "",ref)
                    raw_rows.append((int(row.attrib["r"]), cells))
                yield sheet, raw_rows, list(decoded[sheet].values)
    finally: decoded.close()


LOCATION = ["source_file", "source_sheet", "source_row", "timestamp_raw", "timestamp_label"]
FLAG_FIELDS = ["dataset", "source_file", "source_sheet", "source_row", "source_column",
               "cell_reference", "timestamp_label", "variable", "raw_value", "source_qc",
               "source_quality_status", "status", "reasons"]


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


@contextmanager
def gz_writer(path, fields):

    with path.open("wb") as binary, gzip.GzipFile(
            filename="", fileobj=binary, mode="wb", mtime=0, compresslevel=6) as gz:
        with io.TextIOWrapper(gz, encoding="utf-8", newline="") as text:
            writer = csv.DictWriter(text, fieldnames=fields)
            writer.writeheader()
            yield writer


def write_csv(path, rows, fields=None):
    fields = fields or list(rows[0])
    with path.open("w", encoding="utf-8-sig", newline="") as text:
        writer = csv.DictWriter(text, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def discover(root):

    found = []
    for number, pattern, station in [(1, "*.csv", True), (2, "*.xlsx", False)]:
        folders = [p for p in root.iterdir() if p.is_dir() and
                   p.name.startswith((f"（{number}）", f"({number})"))]
        if len(folders) != 1:
            raise ValueError(f"Expected one attachment ({number}) folder directly inside {root}")
        paths = sorted(p for p in folders[0].rglob(pattern)
                       if not p.name.startswith(("._", "~$")) and "__MACOSX" not in p.parts)
        if not paths:
            raise ValueError(f"No {pattern} inputs inside {folders[0]}")
        for path in paths:
            if station and not re.fullmatch(r"(MET|RADM|SOIL|FLUX)_Gyirong_20[0-9]{2}\.csv", path.name):
                raise ValueError(f"Unrecognized station CSV (not silently skipped): {path}")
            year_match = re.search(r"(20[0-9]{2})(?:\.csv|\.xlsx)$", path.name)
            if not year_match:
                raise ValueError(f"No source year in filename: {path.name}")
            found.append(dict(path=path, station=station, year=int(year_match[1]),
                              dataset=path.name.split("_")[0] if station else "daily_temperature"))
    return found


def read_sheets(source):

    path, station = source["path"], source["station"]
    if station:
        with path.open(encoding="utf-8-sig", newline="") as text:
            reader = csv.reader(text)
            header, data = next(reader), list(reader)
        if len(set(header)) != len(header) or any(len(row) != len(header) for row in data):
            raise ValueError(f"Duplicate header or inconsistent row width: {path.name}")
        columns = [j for j, name in enumerate(header) if j and not name.startswith("QC_")]
        qcmap = {}
        for j in columns:
            wanted = "qc_" + header[j].split("(")[0].lower()
            hits = [k for k, name in enumerate(header) if name.lower() == wanted]
            if len(hits) > 1:
                raise ValueError(f"Ambiguous QC column: {path.name}, {header[j]}")
            qcmap[j] = hits[0] if hits else None
        if {k for k, name in enumerate(header) if name.startswith("QC_")} != set(
                k for k in qcmap.values() if k is not None):
            raise ValueError(f"Unmapped QC column: {path.name}")
        yield "", header, data, columns, qcmap, [
            (i + 2, parse_label(row[0]), {}) for i, row in enumerate(data)]
    else:
        for sheet, raw, decoded in excel_rows(path):
            if not raw or raw[0][0] != 1:
                raise ValueError(f"Missing row-1 header: {path.name}/{sheet}")
            letters = list(raw[0][1])
            header = [raw[0][1][letter][0] for letter in letters]
            if header[0] != "时间" or len(set(header)) != len(header):
                raise ValueError(f"Unexpected daily header: {path.name}/{sheet}")
            data, extra = [], []
            for rn, cells in raw[1:]:
                row = [cells.get(c, ("", "missing", "", f"{c}{rn}"))[0] for c in letters]
                data.append(row)
                value = decoded[rn - 1][0]
                label = value if isinstance(value, dt.datetime) else None
                extra.append((rn, label, {j: cells.get(c, ("", "missing", "", f"{c}{rn}"))
                                          for j, c in enumerate(letters)}))
            yield sheet, header, data, list(range(1, len(header))), {
                j: None for j in range(1, len(header))}, extra


def disposition(a):

    part = a["primary_disposition"]
    if part == "quarantine_structure":
        return "missing" if a["reasons"].split(";")[0] in {
            "empty_token", "declared_sentinel", "nonfinite"} else "excluded_parse"
    if part == "quarantine_source_qc":
        return {"source_flag_8_missing": "missing", "source_flag_2_wrong": "excluded_source_error"}.get(
            a["source_quality_status"], "excluded_qc_policy")
    if part == "quarantine_physical":
        return "excluded_physical"
    if part == "hold_time_label_conflict":
        return "excluded_time"
    if a["physical_status"] != "pass" or a["source_qc_presence"] == "not_provided" or a["conditional_flags"]:
        return "retained_flagged"
    return "retained"


def flag_reasons(a):
    reasons = [x for x in a["reasons"].split(";")
               if x and x != "unknown_or_invalid_support"]
    if a["source_qc_presence"] == "not_provided":
        reasons.append("source_qc_not_provided")
    reasons.extend(a["conditional_flags"].split(";") if a["conditional_flags"] else [])
    return ";".join(dict.fromkeys(reasons))


def main(input_root, output_dir):
    root, out = input_root.resolve(), output_dir.resolve()
    inputs = discover(root)
    if out.exists():
        raise ValueError(f"Output already exists; choose a new directory to preserve prior results: {out}")
    if out == root or root in out.parents:
        raise ValueError("Put outputs outside the original attachment directory")
    started = dt.datetime.now(dt.timezone.utc).isoformat()
    before = {s["path"]: sha(s["path"]) for s in inputs}

    label_counts = collections.defaultdict(collections.Counter)
    schemas = {}
    for source in inputs:
        for sheet, header, data, columns, qcmap, extra in read_sheets(source):
            dataset = source["dataset"]
            if dataset in schemas and schemas[dataset] != header:
                raise ValueError(f"Schema drift in {source['path'].name}; explicit mapping needed")
            schemas[dataset] = header
            label_counts[dataset].update(t for _, t, _ in extra if t is not None)
    out.mkdir(parents=True)
    (out / "station").mkdir()
    (out / "daily").mkdir()
    inventory, summaries, dictionary, diagnostics, time_issues, gaps = [], [], {}, [], [], []
    total_status = collections.Counter()
    totals = collections.Counter()
    with ExitStack() as stack:
        flags = stack.enter_context(gz_writer(out / "quality_flags.csv.gz", FLAG_FIELDS))
        writers = {}
        for dataset, header in schemas.items():
            target = (out / "daily/temperature.csv.gz" if dataset == "daily_temperature"
                      else out / f"station/{dataset}.csv.gz")
            writers[dataset] = stack.enter_context(gz_writer(target, LOCATION + header[1:]))
        for source in inputs:
            path, station, year, dataset = (source[k] for k in ["path", "station", "year", "dataset"])
            for sheet, header, data, columns, qcmap, extra in read_sheets(source):
                labels = [e[1] for e in extra]
                statuses, _ = label_statuses(labels, data, year, station)
                for i, label in enumerate(labels):
                    if label is not None and label_counts[dataset][label] > 1 and "duplicate" not in statuses[i]:
                        statuses[i] = "cross_file_duplicate" if statuses[i] == "unique" else statuses[i] + "+cross_file_duplicate"
                stats = {j: dict(status=collections.Counter(), qc=collections.Counter(),
                                 reasons=collections.Counter(), present=0, physical_pass=0,
                                 first="", last="", minimum=None, maximum=None) for j in columns}
                cleaned = []
                for i, row in enumerate(data):
                    rn, label, cell_info = extra[i]
                    label_text = (label.isoformat() + ("+08:00" if station else "")) if label else ""
                    result = dict(source_file=path.name, source_sheet=sheet, source_row=rn,
                                  timestamp_raw=row[0], timestamp_label=label_text)
                    if statuses[i] != "unique":
                        time_issues.append(dict(dataset=dataset, source_file=path.name, source_sheet=sheet,
                                                source_row=rn, timestamp_raw=row[0], status=statuses[i]))
                    for j in columns:
                        qi = qcmap[j]
                        qraw = row[qi] if qi is not None else ""
                        a = assess(row[j], header[j], "provided" if qi is not None else "not_provided",
                                   qraw, statuses[i], station)
                        ref = str(j + 1)
                        if not station:
                            _, typ, formula, ref = cell_info[j]
                            if formula or typ in {"b", "e"}:
                                a.update(candidate_value=None, primary_disposition="quarantine_structure",
                                         structural_status="fail", reasons="formula_not_evaluated" if formula else "nonnumeric_xml_type")
                        value, status = a["candidate_value"], disposition(a)
                        result[header[j]] = value
                        if qi is not None:
                            result[header[qi]] = qraw
                        s = stats[j]
                        s["status"][status] += 1
                        s["qc"][qraw if qi is not None else "not_provided"] += 1
                        s["present"] += value is not None
                        s["physical_pass"] += value is not None and a["physical_status"] == "pass"
                        if value is not None:
                            s["first"] = s["first"] or label_text
                            s["last"] = label_text
                            s["minimum"] = value if s["minimum"] is None else min(s["minimum"], value)
                            s["maximum"] = value if s["maximum"] is None else max(s["maximum"], value)
                        if status != "retained":
                            reason = flag_reasons(a)
                            s["reasons"].update(reason.split(";"))
                            flags.writerow(dict(dataset=dataset, source_file=path.name, source_sheet=sheet,
                                source_row=rn, source_column=j + 1, cell_reference=ref,
                                timestamp_label=label_text, variable=header[j], raw_value=row[j], source_qc=qraw,
                                source_quality_status=a["source_quality_status"], status=status, reasons=reason))
                        key = dataset, header[j]
                        unit_match = re.search(r"\((.*)\)$", header[j]) if station else None
                        dictionary[key] = dict(dataset=dataset, variable=header[j],
                            unit=unit_match[1] if unit_match else "unknown_in_local_workbook",
                            unit_basis="source_header" if unit_match else "not_declared",
                            qc_column=header[qi] if qi is not None else "not_provided",
                            physical_rule=physical(header[j], station)[0],
                            timezone="UTC+08:00" if station else "not_declared",
                            time_support="label_only; interval/aggregation_unconfirmed" if station else "date_only; day_boundary_unconfirmed")
                    writers[dataset].writerow(result)
                    cleaned.append(result)
                valid_labels = sorted(set(t for t in labels if t is not None))
                step = dt.timedelta(hours=1) if station else dt.timedelta(days=1)
                for left, right in zip(valid_labels, valid_labels[1:]):
                    if right - left > step:
                        gaps.append(dict(dataset=dataset, source_file=path.name, source_sheet=sheet,
                                         start=(left + step).isoformat(), end=(right - step).isoformat(),
                                         nominal_missing_labels=int((right - left) / step) - 1,
                                         basis="between_existing_source_labels; no rows inserted"))
                for j in columns:
                    s = stats[j]
                    if sum(s["status"].values()) != len(data):
                        raise RuntimeError("Cell counts do not close")
                    summaries.append(dict(dataset=dataset, source_file=path.name, source_sheet=sheet,
                        variable=header[j], cells=len(data), analysis_present=s["present"],
                        analysis_missing=len(data) - s["present"],
                        necessary_physical_screen_present=s["physical_pass"],
                        first_analysis_label=s["first"], last_analysis_label=s["last"],
                        min_analysis_value=s["minimum"], max_analysis_value=s["maximum"],
                        statuses=json.dumps(s["status"], ensure_ascii=False, sort_keys=True),
                        source_qc_counts=json.dumps(s["qc"], ensure_ascii=False, sort_keys=True),
                        reason_counts=json.dumps(s["reasons"], ensure_ascii=False, sort_keys=True)))
                    total_status.update(s["status"])
                    totals.update(cells=len(data), analysis_present=s["present"], physical_screen_present=s["physical_pass"])
                inventory.append(dict(dataset=dataset, source_file=path.name, source_sheet=sheet,
                    input_path=path.relative_to(root).as_posix(), sha256=before[path],
                    encoding="utf-8-sig strict" if station else "OOXML raw numeric tokens + decoded date labels",
                    rows=len(data), measurement_columns=len(columns),
                    first_label=valid_labels[0].isoformat() if valid_labels else "",
                    last_label=valid_labels[-1].isoformat() if valid_labels else "",
                    time_status_counts=json.dumps(collections.Counter(statuses), ensure_ascii=False),
                    output_path=f"station/{dataset}.csv.gz" if station else "daily/temperature.csv.gz"))
                diagnostics.extend(diagnose(pd.DataFrame(cleaned), [header[j] for j in columns],
                                            dataset=dataset, step_hours=1 if station else 24))
                print(f"{path.name}: {len(data):,} rows, {len(columns)} measurement columns", flush=True)
    if not all(sha(s["path"]) == before[s["path"]] for s in inputs):
        raise RuntimeError("A source input changed during cleaning")
    write_csv(out / "input_inventory.csv", inventory)
    write_csv(out / "variable_summary.csv", summaries)
    write_csv(out / "variable_dictionary.csv", list(dictionary.values()))
    write_csv(out / "internal_diagnostics.csv", diagnostics,
              ["dataset", "source_file", "check", "variable", "count", "start", "end", "value_before", "value_after", "detail"])
    write_csv(out / "time_issues.csv", time_issues,
              ["dataset", "source_file", "source_sheet", "source_row", "timestamp_raw", "status"])
    write_csv(out / "record_gaps.csv", gaps,
              ["dataset", "source_file", "source_sheet", "start", "end", "nominal_missing_labels", "basis"])
    summary = dict(started_utc=started, finished_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
        scope="official attachments 1 and 2 only", input_files=len(inputs),
        station_files=sum(s["station"] for s in inputs), daily_files=sum(not s["station"] for s in inputs),
        station_rows=sum(x["rows"] for x in inventory if x["dataset"] != "daily_temperature"),
        daily_rows=sum(x["rows"] for x in inventory if x["dataset"] == "daily_temperature"),
        **totals, statuses=dict(total_status), source_inputs_unchanged=True,
        imputed_here=0, clipped_here=0, unit_conversions=0, inserted_observation_rows=0,
        model_or_factor_computation=False, time_issue_rows=len(time_issues), record_gap_intervals=len(gaps),
        code_sha256={p.name: sha(p) for p in [BASE / "clean_data.py", BASE / "quality_policy.py",
                                            BASE / "quality_diagnostics.py", BASE / "quality_rules.json"]},
        environment=dict(python=sys.version.split()[0], pandas=pd.__version__,
                         openpyxl=__import__("openpyxl").__version__),
        output_sha256={p.relative_to(out).as_posix(): sha(p) for p in sorted(out.rglob("*")) if p.is_file()},
        limitations=["Retained numbers have passed the declared cleaning policy, not every possible quality check.",
                     "No-source-QC and physical-definition-pending values are retained and flagged.",
                     "Station time labels are UTC+8; interval endpoints/aggregation remain unconfirmed.",
                     "Daily units/day boundary/spatial extraction are not declared locally; no conversion.",
                     "PA is not treated as MET pressure; precipitation is not accumulated here."])
    (out / "run_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: summary[k] for k in ["input_files", "station_rows", "daily_rows", "cells", "analysis_present", "statuses"]},
                     ensure_ascii=False, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=pathlib.Path, required=True, help="Official 附件 directory containing （1） and （2）")
    parser.add_argument("--output-dir", type=pathlib.Path, required=True, help="New output directory outside original attachments")
    args = parser.parse_args()
    main(args.input_root, args.output_dir)
