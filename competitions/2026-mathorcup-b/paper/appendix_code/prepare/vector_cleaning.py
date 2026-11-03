from __future__ import annotations

import codecs
import csv
import gzip
import hashlib
import json
import re
import warnings
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd
import pyogrio
import shapely


_ATTACHMENT_RE = re.compile(r"^[（(]\s*([368])\s*[）)]")
_ATTACHMENTS = (3, 6, 8)
_OUTPUT_NAMES = {
    "JilongFlood_Buffer_200km": "flood_buffer_200km",
    "RGI_V6_JilongFlood_Buffer_200km": "rgi_v6_glaciers",
    "RGI_V7_JilongFlood_Buffer_200km": "rgi_v7_glaciers",
    "西藏日喀则市吉隆县泥石流灾害区历史灾点数据": "historical_hazards",
    "西藏日喀则市吉隆县第七次人口普查数据": "census_population",
}
_ID_COLUMNS = {
    "rgiid", "rgi_id", "glimsid", "glims_id", "code", "orig_fid",
    "区划码", "县级码",
}
_RESERVED = {
    "source_row", "source_fid", "source_file", "qc_source_presence",
    "qc_source_encoding", "qc_geometry_status",
}


def _sha(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def _relative(path: Path, root: Path) -> str:
    return path.relative_to(root).as_posix()


def _write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def _write_csv(path: Path, frame: pd.DataFrame) -> None:


    frame.to_csv(
        path, index=False, encoding="utf-8", na_rep="", lineterminator="\n",
        compression={"method": "gzip", "mtime": 0},
        date_format="%Y-%m-%dT%H:%M:%S",
    )


def _dbf_layout(path: Path) -> dict[str, Any]:
    data = path.read_bytes()
    if len(data) < 33:
        raise ValueError("DBF header is truncated")
    count = int.from_bytes(data[4:8], "little")
    header_size = int.from_bytes(data[8:10], "little")
    record_size = int.from_bytes(data[10:12], "little")
    if header_size < 33 or record_size < 1:
        raise ValueError("DBF header/record length is invalid")
    if len(data) < header_size + count * record_size:
        raise ValueError("DBF declared records exceed available source bytes")
    fields: list[dict[str, Any]] = []
    offset = 1
    position = 32
    while position < header_size and data[position] != 0x0D:
        descriptor = data[position:position + 32]
        if len(descriptor) != 32:
            raise ValueError("DBF field descriptor is truncated")
        field = {
            "name_bytes": descriptor[:11].split(b"\0", 1)[0],
            "type": chr(descriptor[11]), "width": int(descriptor[16]),
            "decimals": int(descriptor[17]), "record_offset": offset,
        }
        fields.append(field)
        offset += field["width"]
        position += 32
    if position >= header_size or data[position] != 0x0D or offset != record_size:
        raise ValueError("DBF field layout does not match declared record length")
    flags = [data[header_size + row * record_size] for row in range(count)]
    if any(flag not in (0x20, 0x2A) for flag in flags):
        raise ValueError("Unrecognised DBF deletion flag; refusing silent row loss")
    return {
        "bytes": data, "version": data[0], "record_count": count,
        "header_size": header_size, "record_size": record_size,
        "language_driver_id": data[29], "fields": fields, "flags": flags,
    }


def _dbf_text_cells(dbf: dict[str, Any]):
    for field in dbf["fields"]:
        yield None, field["name_bytes"]
    for row in range(dbf["record_count"]):
        start = dbf["header_size"] + row * dbf["record_size"]
        for field in dbf["fields"]:
            if field["type"] == "C":
                offset = start + field["record_offset"]
                yield row, dbf["bytes"][offset:offset + field["width"]].rstrip(b" \0")


def _decode_probe(dbf: dict[str, Any], encoding: str) -> dict[str, Any]:
    failures = names = values = 0
    for row, raw in _dbf_text_cells(dbf):
        if row is None:
            names += 1
        else:
            values += 1
        try:
            raw.decode(encoding, errors="strict")
        except UnicodeDecodeError:
            failures += 1
    return {
        "encoding": encoding, "strict_decode_errors": failures,
        "field_names_checked": names, "character_values_checked": values,
        "all_text_bytes_decodable": failures == 0,
    }


def _normalise_encoding(label: str) -> str:
    aliases = {"65001": "utf-8", "936": "gbk", "cp936": "gbk"}
    return codecs.lookup(aliases.get(label.strip().lower(), label.strip())).name


def _encoding_decision(path: Path, dbf: dict[str, Any]) -> dict[str, Any]:
    cpg = path.with_suffix(".cpg")
    label = cpg.read_text(encoding="utf-8-sig").strip() if cpg.exists() else None
    probes = [_decode_probe(dbf, name) for name in ("utf-8", "gbk", "gb18030")]
    evidence = [
        f"DBF language driver byte: 0x{dbf['language_driver_id']:02X}",
        "Strict decoding checks include every DBF field name and character cell.",
        "GDAL reports its Unicode transfer encoding, which is not proof of source byte encoding.",
    ]
    selected = None
    reason = None
    if label:
        try:
            candidate = _normalise_encoding(label)
            probe = _decode_probe(dbf, candidate)
            if probe["all_text_bytes_decodable"]:
                selected, reason = candidate, "source_cpg_and_strict_byte_check"
            else:
                evidence.append("Source CPG conflicts with DBF text bytes; no automatic override.")
        except LookupError:
            evidence.append("Source CPG label is not recognised by this Python codec registry.")
    elif dbf["language_driver_id"] == 0x4D:

        if _decode_probe(dbf, "gbk")["all_text_bytes_decodable"]:
            selected, reason = "gbk", "dbf_ldid_0x4d_cp936_and_strict_byte_check"
            evidence.append("DBF LDID 0x4D identifies Chinese CP936/GBK.")
    elif all(all(byte < 128 for byte in raw) for _, raw in _dbf_text_cells(dbf)):
        selected, reason = "ascii", "all_field_names_and_character_bytes_are_ascii"
    unresolved = selected is None
    return {
        "source_cpg": label, "source_cpg_present": cpg.exists(),
        "dbf_language_driver_id_hex": f"0x{dbf['language_driver_id']:02X}",
        "selected_source_encoding": selected,
        "reader_encoding": selected if selected else "ISO-8859-1",
        "status": "unresolved_raw_bytes" if unresolved else "source_supported",
        "decision_basis": reason, "evidence": evidence, "candidate_probes": probes,
        "fallback_policy": (
            "Undecided bytes use reversible ISO-8859-1 transport only; they are not interpreted as human text. "
            "Raw DBF records and candidate decodings are retained. Confirm source encoding before text use."
            if unresolved else None
        ),
        "recommended_options_if_unresolved": [
            "Retain raw bytes and candidate decodings until the producer supplies encoding metadata.",
            "Confirm a candidate against independent source labels before selecting it in a later version.",
        ] if unresolved else [],
    }


def _xml_evidence(path: Path, field_names: list[str]) -> dict[str, Any]:
    xml_path = path.with_suffix(".shp.xml")
    if not xml_path.exists():
        return {"present": False}
    data = xml_path.read_bytes()
    try:
        root = ET.fromstring(data)
        all_text = "\n".join(text for text in root.itertext() if text.strip())
        matches = [name for name in field_names if any(ord(c) > 127 for c in name) and name in all_text]
        declared = re.search(br'<\?xml[^>]*encoding=[\'\"]([^\'\"]+)', data[:256])
        return {
            "present": True, "sha256": _sha(xml_path), "xml_parse": "pass",
            "xml_document_encoding": declared.group(1).decode("ascii") if declared else None,
            "matched_non_ascii_dbf_field_names": matches,
            "arcgis_process_entries": sum(element.tag.rsplit("}", 1)[-1] == "Process" for element in root.iter()),
            "interpretation": "XML document encoding is independent of DBF attribute encoding. Processing history corroborates field labels only; it is not a current attribute codebook.",
        }
    except (ET.ParseError, UnicodeError) as error:
        return {"present": True, "sha256": _sha(xml_path), "xml_parse": "pending", "error": str(error)}


def _write_raw_dbf(path: Path, dbf: dict[str, Any]) -> None:
    with path.open("wb") as binary:
        with gzip.GzipFile(filename="", mode="wb", fileobj=binary, mtime=0) as compressed:
            import io
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as stream:
                writer = csv.writer(stream, lineterminator="\n")
                writer.writerow(["source_row", "source_fid", "dbf_deleted", "record_hex"])
                for row, flag in enumerate(dbf["flags"]):
                    start = dbf["header_size"] + row * dbf["record_size"]
                    writer.writerow([row + 1, row, flag == 0x2A, dbf["bytes"][start:start + dbf["record_size"]].hex()])


def _write_candidate_decodings(path: Path, dbf: dict[str, Any]) -> None:
    rows = []
    for index, field in enumerate(dbf["fields"]):
        name_bytes = field["name_bytes"]
        for encoding in ("utf-8", "gbk", "gb18030"):
            try:
                text, state = name_bytes.decode(encoding, errors="strict"), "decoded_candidate_only"
            except UnicodeDecodeError:
                text, state = "", "strict_decode_failed"
            rows.append({"source_row": 0, "field_index": index, "encoding": encoding, "raw_hex": name_bytes.hex(), "decoded_candidate": text, "status": state})
        if field["type"] != "C":
            continue
        for row in range(dbf["record_count"]):
            start = dbf["header_size"] + row * dbf["record_size"] + field["record_offset"]
            raw = dbf["bytes"][start:start + field["width"]].rstrip(b" \0")
            for encoding in ("utf-8", "gbk", "gb18030"):
                try:
                    text, state = raw.decode(encoding, errors="strict"), "decoded_candidate_only"
                except UnicodeDecodeError:
                    text, state = "", "strict_decode_failed"
                rows.append({"source_row": row + 1, "field_index": index, "encoding": encoding, "raw_hex": raw.hex(), "decoded_candidate": text, "status": state})
    _write_csv(path, pd.DataFrame(rows, columns=["source_row", "field_index", "encoding", "raw_hex", "decoded_candidate", "status"]))


def _geometry_quality(frame: gpd.GeoDataFrame) -> pd.DataFrame:
    geometries = frame.geometry.array
    missing = np.asarray(shapely.is_missing(geometries))
    empty = np.asarray(shapely.is_empty(geometries))
    valid = np.asarray(shapely.is_valid(geometries))
    status = np.full(len(frame), "valid", dtype=object)
    status[~valid] = "invalid"
    status[empty] = "empty"
    status[missing] = "missing"
    wkbs = shapely.to_wkb(geometries)
    return pd.DataFrame({
        "source_row": frame["source_row"].to_numpy(), "source_fid": frame["source_fid"].to_numpy(),
        "geometry_status": status, "geometry_type": frame.geometry.geom_type.to_numpy(),
        "validity_reason": shapely.is_valid_reason(geometries),
        "has_z": shapely.has_z(geometries), "has_m": shapely.has_m(geometries),
        "coordinate_dimension": shapely.get_coordinate_dimension(geometries),
        "geometry_wkb_sha256": [hashlib.sha256(wkb).hexdigest() if wkb is not None else None for wkb in wkbs],
        "source_qc_presence": "not_provided",
        "qc_interpretation": "structural diagnostics only; source scientific QC is not supplied",
    })


def _attribute_quality(frame: gpd.GeoDataFrame, columns: list[str], dbf: dict[str, Any], encoding: str):
    layout = {field["name_bytes"].decode(encoding): field for field in dbf["fields"]}
    profile = []
    duplicates = []
    issues = []
    for column in columns:
        series = frame[column]
        field = layout[column]
        blanks = series.map(lambda value: isinstance(value, str) and not value.strip())
        profile.append({
            "field": column, "decoded_dtype": str(series.dtype), "dbf_type": field["type"],
            "dbf_width": field["width"], "dbf_decimals": field["decimals"],
            "rows": len(frame), "null_count": int(series.isna().sum()),
            "blank_string_count": int(blanks.sum()), "non_null_distinct_count": int(series.nunique(dropna=True)),
            "source_qc_presence": "not_provided", "missing_value_policy": "preserved; no null-to-zero or undeclared sentinel substitution",
        })
        if column.casefold() in _ID_COLUMNS:
            duplicated = series.notna() & series.duplicated(keep=False)
            for position in np.flatnonzero(duplicated.to_numpy()):
                duplicates.append({
                    "source_row": int(frame.iloc[position]["source_row"]), "source_fid": int(frame.iloc[position]["source_fid"]),
                    "field": column, "value": str(series.iloc[position]), "action": "retained_without_deduplication",
                    "key_status": "identifier_candidate; uniqueness is not assumed",
                })

        sentinel_like = series.astype("string").isin(["-9999", "-9999.0", "-9999999", "9999.9"])
        for position in np.flatnonzero(sentinel_like.to_numpy(dtype=bool, na_value=False)):
            issues.append({
                "source_row": int(frame.iloc[position]["source_row"]), "source_fid": int(frame.iloc[position]["source_fid"]),
                "field": column, "observed_value": str(series.iloc[position]),
                "observation": "placeholder_like_token_uninterpreted", "action": "retained; field-specific source definition required before recoding",
            })
    return (
        pd.DataFrame(profile),
        pd.DataFrame(duplicates, columns=["source_row", "source_fid", "field", "value", "action", "key_status"]),
        pd.DataFrame(issues, columns=["source_row", "source_fid", "field", "observed_value", "observation", "action"]),
    )


def _verify_layer(expected: gpd.GeoDataFrame, path: Path, layer: str) -> dict[str, Any]:
    actual = pyogrio.read_dataframe(path, layer=layer)
    if len(expected) != len(actual):
        raise ValueError(f"Readback row-count mismatch for {layer}")
    if expected.crs != actual.crs:
        raise ValueError(f"Readback CRS mismatch for {layer}")
    if list(expected.columns) != list(actual.columns):
        raise ValueError(f"Readback attribute/geometry column mismatch for {layer}")
    pd.testing.assert_frame_equal(
        expected.drop(columns=expected.geometry.name).reset_index(drop=True),
        actual.drop(columns=actual.geometry.name).reset_index(drop=True),
        check_dtype=False, check_exact=True,
    )
    before, after = shapely.to_wkb(expected.geometry.array), shapely.to_wkb(actual.geometry.array)
    if not all(left == right for left, right in zip(before, after)):
        raise ValueError(f"Readback WKB mismatch (coordinates/order/dimension changed) for {layer}")
    return {
        "status": "pass", "rows": len(actual), "all_attribute_values_equal": True,
        "row_order_and_stable_ids_equal": True, "crs_equivalent": True,
        "geometry_wkb_exact": True, "coordinates_rings_parts_and_z_unchanged": True,
        "dtype_policy": "GeoPackage integer/date storage widths can differ; all decoded values and missingness compared exactly",
    }


def _clean_layer(path: Path, attachment: int, input_root: Path, output_root: Path) -> dict[str, Any]:
    source_file = _relative(path, input_root)
    slug = _OUTPUT_NAMES.get(path.stem, "layer_" + hashlib.sha256(source_file.encode()).hexdigest()[:12])
    target = output_root / f"attachment{attachment:02d}" / "vectors"
    target.mkdir(parents=True, exist_ok=True)
    required = [path, path.with_suffix(".shx"), path.with_suffix(".dbf"), path.with_suffix(".prj")]
    missing = [item.suffix for item in required if not item.is_file()]
    if missing:
        raise ValueError(f"Required shapefile components are missing: {missing}")
    dbf = _dbf_layout(path.with_suffix(".dbf"))
    encoding = _encoding_decision(path, dbf)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        info = pyogrio.read_info(path, encoding=encoding["reader_encoding"])
        frame = pyogrio.read_dataframe(path, encoding=encoding["reader_encoding"], fid_as_index=True)
    source_columns = [column for column in frame.columns if column != frame.geometry.name]
    if any(column in _RESERVED for column in source_columns):
        raise ValueError("Source attributes collide with reserved provenance columns; explicit name mapping is required")
    if len(source_columns) != len(set(source_columns)):
        raise ValueError("Duplicate decoded field names are ambiguous; raw DBF bytes remain authoritative")
    fields_encoding = encoding["selected_source_encoding"] or "latin-1"
    if source_columns != [field["name_bytes"].decode(fields_encoding) for field in dbf["fields"]]:
        raise ValueError("Decoded DBF field order does not match the GIS reader")
    active_fids = [row for row, flag in enumerate(dbf["flags"]) if flag == 0x20]
    if list(frame.index) != active_fids:
        raise ValueError("GIS reader omitted/reordered source DBF records; refusing silent row loss")
    if len(active_fids) != dbf["record_count"]:
        raise ValueError("Deleted DBF records require an explicit retained-record import; no source rows are silently dropped")
    frame.insert(0, "source_fid", np.asarray(frame.index, dtype="int64"))
    frame.insert(0, "source_row", frame["source_fid"].to_numpy() + 1)
    frame.insert(2, "source_file", source_file)
    frame["qc_source_presence"] = "not_provided"
    frame["qc_source_encoding"] = encoding["status"]
    quality = _geometry_quality(frame)
    frame["qc_geometry_status"] = quality["geometry_status"].to_numpy()


    geometry_name = frame.geometry.name
    frame = frame[[column for column in frame.columns if column != geometry_name] + [geometry_name]]
    frame = frame.reset_index(drop=True)
    profile, duplicate_ids, issues = _attribute_quality(frame, source_columns, dbf, fields_encoding)
    analysis = frame.loc[frame["qc_geometry_status"].eq("valid")].copy()
    files = {
        "geopackage": target / f"{slug}.gpkg",
        "attributes": target / f"{slug}_attributes.csv.gz",
        "geometry_quality": target / f"{slug}_geometry_quality.csv.gz",
        "attribute_profile": target / f"{slug}_attribute_profile.csv.gz",
        "duplicate_ids": target / f"{slug}_duplicate_ids.csv.gz",
        "attribute_observations": target / f"{slug}_attribute_observations.csv.gz",
        "raw_dbf_records": target / f"{slug}_raw_dbf_records.csv.gz",
        "metadata": target / f"{slug}_metadata.json",
    }
    if any(item.exists() for item in files.values()):
        raise FileExistsError(f"Refusing to overwrite an existing vector layer output: {slug}")
    pyogrio.write_dataframe(frame, files["geopackage"], layer="all_features", driver="GPKG", promote_to_multi=False)
    pyogrio.write_dataframe(analysis, files["geopackage"], layer="analysis_valid_geometry", driver="GPKG", promote_to_multi=False, append=False)
    checks = {
        "all_features": _verify_layer(frame, files["geopackage"], "all_features"),
        "analysis_valid_geometry": _verify_layer(analysis, files["geopackage"], "analysis_valid_geometry"),
    }
    _write_csv(files["attributes"], frame.drop(columns=frame.geometry.name))
    _write_csv(files["geometry_quality"], quality)
    _write_csv(files["attribute_profile"], profile)
    _write_csv(files["duplicate_ids"], duplicate_ids)
    _write_csv(files["attribute_observations"], issues)
    _write_raw_dbf(files["raw_dbf_records"], dbf)
    if encoding["status"] == "unresolved_raw_bytes":
        files["candidate_decodings"] = target / f"{slug}_candidate_decodings.csv.gz"
        _write_candidate_decodings(files["candidate_decodings"], dbf)
    source_components = []
    for component in sorted(path.parent.iterdir()):
        if component.is_file() and (component.stem == path.stem or component.name.startswith(path.name + ".")):
            source_components.append({"path": _relative(component, input_root), "bytes": component.stat().st_size, "sha256": _sha(component)})
    geometry_counts = {str(key): int(value) for key, value in quality["geometry_status"].value_counts().items()}
    metadata = {
        "attachment": attachment, "source_file": source_file, "output_name": slug,
        "source_components": source_components, "source_rows": dbf["record_count"],
        "retained_rows": len(frame), "analysis_geometry_rows": len(analysis),
        "source_fid_definition": "GDAL Shapefile FID, zero-based physical DBF/SHP record position",
        "source_row_definition": "one-based physical DBF/SHP record position = source_fid + 1",
        "source_attribute_columns": source_columns, "source_qc_presence": "not_provided",
        "source_qc_note": "No dedicated scientific QC field is documented in the supplied vector attributes; inventory status/classification attributes are retained without reinterpretation.",
        "source_crs_wkt": frame.crs.to_wkt() if frame.crs else None,
        "source_crs_authority": list(frame.crs.to_authority()) if frame.crs and frame.crs.to_authority() else None,
        "source_prj_sha256": _sha(path.with_suffix(".prj")),
        "crs_status": "source_metadata_preserved" if frame.crs else "unknown_not_guessed",
        "source_geometry_type": str(info["geometry_type"]), "geometry_status_counts": geometry_counts,
        "geometry_types": {str(key): int(value) for key, value in frame.geometry.geom_type.value_counts(dropna=False).items()},
        "geometry_with_z": int(quality["has_z"].sum()), "geometry_with_m": int(quality["has_m"].sum()),
        "encoding": encoding, "xml_evidence": _xml_evidence(path, source_columns),
        "dbf_layout": {
            "version": dbf["version"], "header_size": dbf["header_size"], "record_size": dbf["record_size"],
            "deleted_record_count": dbf["flags"].count(0x2A),
            "fields": [{**{key: value for key, value in field.items() if key != "name_bytes"}, "name": field["name_bytes"].decode(fields_encoding), "name_bytes_hex": field["name_bytes"].hex()} for field in dbf["fields"]],
        },
        "policies": {
            "all_features": "all source rows and decoded attributes; original geometry including invalid geometry retained",
            "analysis_valid_geometry": "only structural geometry screening; excludes missing/empty/invalid geometry without repairing or splitting",
            "geometry_repair": "none", "coordinate_transformation": "none", "deduplication": "none",
            "missing_values": "no fill, zero substitution, global sentinel recoding or undocumented missingness inference",
            "units": "source attributes retained; no field units are guessed",
            "csv": "UTF-8 compressed CSV, decoded nulls exported as empty fields; profiles and raw DBF bytes preserve interpretation evidence",
            "raw_dbf": "record_hex includes the full deletion marker and fixed-width record bytes; use dbf_layout to decode without losing whitespace",
        },
        "duplicate_identifier_rows": len(duplicate_ids), "attribute_observations": len(issues),
        "reader_warnings": [str(item.message).replace(str(path), source_file) for item in caught],
        "readback_verification": checks,
        "engine_versions": {"geopandas": gpd.__version__, "pyogrio": pyogrio.__version__, "gdal": pyogrio.__gdal_version_string__, "shapely": shapely.__version__, "geos": shapely.geos_version_string, "pandas": pd.__version__},
        "outputs": {key: _relative(value, output_root) for key, value in files.items()},
        "remaining_semantic_questions": [
            "Attribute units, categorical codes and observation/interpretation date meanings require the field-specific source documentation before analytical use."
        ],
    }
    if attachment == 6:
        metadata["remaining_semantic_questions"].append("The source is Polygon Z, despite being called historical disaster points. Resolution/imagery/river fields contain potentially inconsistent semantics; values are retained exactly.")
    if attachment == 8:
        metadata["remaining_semantic_questions"].append("合计 and 人口数 are distinct source fields; population/household scope is not equated or overwritten.")
    _write_json(files["metadata"], metadata)
    return {
        "attachment": attachment, "source_file": source_file, "name": slug,
        "status": "complete" if encoding["status"] == "source_supported" else "pending_encoding",
        "rows": len(frame), "analysis_geometry_rows": len(analysis),
        "geometry_status_counts": geometry_counts, "source_encoding": encoding["selected_source_encoding"],
        "source_qc_presence": "not_provided", "duplicate_identifier_rows": len(duplicate_ids),
        "readback_verification": checks,
        "outputs": metadata["outputs"],
        "output_fingerprints": {key: {"sha256": _sha(value), "bytes": value.stat().st_size} for key, value in files.items()},
    }


def clean_vectors(input_root: Path, output_root: Path) -> dict[str, Any]:


    input_root, output_root = Path(input_root).resolve(), Path(output_root).resolve()
    if not input_root.is_dir():
        raise FileNotFoundError("input_root must be an existing local attachment directory")
    if output_root == input_root or input_root.is_relative_to(output_root) or output_root.is_relative_to(input_root):
        raise ValueError("Input and output roots must be separate directories, with neither containing the other")
    discovered = []
    for path in sorted(input_root.rglob("*.shp")):
        match = next((match for part in path.relative_to(input_root).parts[:-1] if (match := _ATTACHMENT_RE.match(part))), None)
        if match:
            discovered.append((int(match.group(1)), path))
    layers, failures = [], []
    for attachment, path in discovered:
        try:
            layers.append(_clean_layer(path, attachment, input_root, output_root))
        except Exception as error:

            failures.append({"attachment": attachment, "source_file": _relative(path, input_root), "status": "pending", "error_type": type(error).__name__, "error": str(error).replace(str(input_root), "<input_root>").replace(str(output_root), "<output_root>")})
    counts = {number: sum(attachment == number for attachment, _ in discovered) for number in _ATTACHMENTS}
    inventory = {f"attachment{number:02d}": {"minimum_expected_layers": 1, "found_layers": counts[number], "source_layers": [_relative(path, input_root) for attachment, path in discovered if attachment == number], "status": "pass" if counts[number] >= 1 else "pending"} for number in _ATTACHMENTS}
    complete = not failures and all(record["status"] == "pass" for record in inventory.values()) and all(layer["status"] == "complete" for layer in layers)
    result = {
        "status": "complete" if complete else "pending", "scope": "attachments 03, 06, 08 vector basic cleaning only",
        "source_qc_presence": "not_provided", "inventory": inventory,
        "layers": layers, "failures": failures,
        "total_retained_rows": sum(layer["rows"] for layer in layers),
        "no_geometry_repair_or_coordinate_transform": True,
        "independent_review": "not_performed_by_this_executor",
        "inventory_path": "quality/vector_inventory.json",
    }

    json.dumps(result, ensure_ascii=False, allow_nan=False)
    inventory_path = output_root / result["inventory_path"]
    inventory_path.parent.mkdir(parents=True, exist_ok=True)
    if inventory_path.exists():
        raise FileExistsError("Refusing to overwrite quality/vector_inventory.json")
    _write_json(inventory_path, result)
    return result
