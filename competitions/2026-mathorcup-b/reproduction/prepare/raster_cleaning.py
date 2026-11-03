from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
from typing import Any, Iterator

import numpy as np
import rasterio
from rasterio.shutil import copy as copy_raster
from rasterio.windows import Window


ATTACHMENTS = (3, 5, 7)
WINDOW_SIZE = 512


LAYOUT_TAG_DOMAINS = {"IMAGE_STRUCTURE", "DERIVED_SUBDATASETS"}


def _json_value(value: Any) -> Any:
    if isinstance(value, (np.integer, np.floating)):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return "NaN" if math.isnan(value) else ("Infinity" if value > 0 else "-Infinity")
    if isinstance(value, dict):
        return {str(k): _json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(v) for v in value]
    return value


def _dump_json(path: Path, data: Any) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(_json_value(data), handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write("\n")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _windows(width: int, height: int) -> Iterator[Window]:

    for row in range(0, height, WINDOW_SIZE):
        for col in range(0, width, WINDOW_SIZE):
            yield Window(col, row, min(WINDOW_SIZE, width - col), min(WINDOW_SIZE, height - row))


def _tags(src: Any, band: int = 0) -> dict[str, dict[str, str]]:
    tags = {"default": src.tags(band)}
    for namespace in src.tag_namespaces(band):
        if namespace not in LAYOUT_TAG_DOMAINS:
            tags[namespace] = src.tags(band, ns=namespace)
    return tags


def _metadata(src: Any) -> dict[str, Any]:
    crs_wkt = src.crs.to_wkt() if src.crs else None
    return {
        "driver": src.driver,
        "width": src.width,
        "height": src.height,
        "count": src.count,
        "dtypes": list(src.dtypes),
        "crs": str(src.crs) if src.crs else None,
        "crs_wkt": crs_wkt,
        "epsg": src.crs.to_epsg() if src.crs else None,
        "transform": list(src.transform),
        "resolution": list(src.res),
        "bounds": list(src.bounds),
        "nodata": _json_value(src.nodatavals),
        "scales": list(src.scales),
        "offsets": list(src.offsets),
        "units": list(src.units),
        "descriptions": list(src.descriptions),
        "color_interpretation": [item.name for item in src.colorinterp],
        "mask_flags": [[item.name for item in band] for band in src.mask_flag_enums],
        "tags": _tags(src),
        "band_tags": [_tags(src, band) for band in src.indexes],
    }


def _grid_signature(metadata: dict[str, Any]) -> str:
    grid = {key: metadata[key] for key in ("crs_wkt", "transform", "width", "height")}
    grid["pixel_registration"] = metadata["tags"]["default"].get("AREA_OR_POINT", "unspecified")
    return hashlib.sha256(json.dumps(grid, sort_keys=True, ensure_ascii=False).encode("utf-8")).hexdigest()


def _copy_metadata(src: Any, dst: Any) -> None:
    dst.scales = src.scales
    dst.offsets = src.offsets
    dst.units = src.units
    dst.descriptions = src.descriptions
    dst.colorinterp = src.colorinterp
    for band in (0, *src.indexes):
        for namespace, tags in _tags(src, band).items():
            if tags:
                dst.update_tags(band, ns=None if namespace == "default" else namespace, **tags)


def _flag(code: str, detail: str, severity: str = "pending") -> dict[str, str]:
    return {"code": code, "severity": severity, "detail": detail}


def _metadata_flags(metadata: dict[str, Any], number: int) -> list[dict[str, str]]:
    flags = []
    if not metadata["crs"]:
        flags.append(_flag("crs_unconfirmed", "Source has no declared CRS; coordinates were not relabelled."))
    if any(value is None for value in metadata["nodata"]):
        flags.append(_flag("nodata_not_declared", "No new sentinel was inferred; source masks and nonfinite values determine validity."))
    if number == 3 and any(not unit for unit in metadata["units"]):
        flags.append(_flag("dhdt_units_unconfirmed", "Provider document identifies elevation-change rate/error but does not state their numeric units; source scale and offset are retained."))
    if number == 5:
        flags.append(_flag("land_cover_class_dictionary_unconfirmed", "Observed integer codes are retained without assigning class names; no provider class dictionary was found in this attachment."))
    if number == 7:
        if any(not unit for unit in metadata["units"]):
            flags.append(_flag("elevation_units_unconfirmed", "Elevation band units are missing; no conversion was applied."))

        flags.append(_flag("vertical_datum_unconfirmed", "A confirmed vertical datum is not supplied by the raster CRS; consult provider documentation before using absolute elevations across sources."))
    transform = metadata["transform"]
    if transform[1] != 0 or transform[3] != 0:
        flags.append(_flag("rotated_grid_preserved", "Source affine rotation/shear is retained.", "info"))
    return flags


def _source_document_notes(folder: Path, number: int) -> list[dict[str, str]]:

    notes = []
    if number == 3:
        for path in sorted(folder.rglob("*.docx")):


            notes.append({"path": path.relative_to(folder).as_posix(), "sha256": _sha256(path), "role": "provider_document; rate/error numeric unit remains unconfirmed"})
    if number == 7:
        for path in sorted(folder.rglob("*.tif.xml")):
            note = {"path": path.relative_to(folder).as_posix(), "sha256": _sha256(path), "role": "provider_spatial_metadata"}
            try:
                text = " ".join(ET.parse(path).getroot().itertext())
                if re.search(r"VDATUM\s*\[\s*['\"]unknown['\"]", text, flags=re.I):
                    note["vertical_datum_evidence"] = "provider XML explicitly declares VDATUM unknown"
            except (ET.ParseError, UnicodeError, OSError) as exc:
                note["read_status"] = f"pending: {type(exc).__name__}"
            notes.append(note)
    return notes


def _new_band_stats(band: int, src: Any) -> dict[str, Any]:
    return {
        "band": band, "dtype": src.dtypes[band - 1],
        "scale": src.scales[band - 1], "offset": src.offsets[band - 1],
        "unit": src.units[band - 1], "sample_count": 0,
        "source_masked_count": 0, "finite_source_masked_count": 0,
        "nonfinite_count": 0, "nonfinite_source_unmasked_count": 0,
        "analysis_valid_count": 0, "analysis_invalid_count": 0,
        "valid_raw_min": None, "valid_raw_max": None,
        "valid_scaled_min": None, "valid_scaled_max": None,
        "valid_negative_count": 0,
    }


def _update_stats(stats: dict[str, Any], values: np.ndarray, source_mask: np.ndarray,
                  analysis_valid: np.ndarray, class_counts: Counter | None) -> None:
    finite = np.isfinite(values)
    source_masked = source_mask == 0
    valid_values = values[analysis_valid]
    stats["sample_count"] += int(values.size)
    stats["source_masked_count"] += int(np.count_nonzero(source_masked))
    stats["finite_source_masked_count"] += int(np.count_nonzero(source_masked & finite))
    stats["nonfinite_count"] += int(np.count_nonzero(~finite))
    stats["nonfinite_source_unmasked_count"] += int(np.count_nonzero(~source_masked & ~finite))
    count = int(valid_values.size)
    stats["analysis_valid_count"] += count
    stats["analysis_invalid_count"] += int(values.size) - count
    if count:
        low, high = valid_values.min().item(), valid_values.max().item()
        stats["valid_raw_min"] = low if stats["valid_raw_min"] is None else min(low, stats["valid_raw_min"])
        stats["valid_raw_max"] = high if stats["valid_raw_max"] is None else max(high, stats["valid_raw_max"])
        stats["valid_negative_count"] += int(np.count_nonzero(valid_values < 0))
        if class_counts is not None:
            unique, counts = np.unique(valid_values, return_counts=True)
            class_counts.update({str(value.item()): int(n) for value, n in zip(unique, counts)})


def _finalize_stats(stats: dict[str, Any], class_counts: Counter | None) -> None:
    if stats["valid_raw_min"] is not None:
        limits = [float(stats[key]) * stats["scale"] + stats["offset"] for key in ("valid_raw_min", "valid_raw_max")]
        stats["valid_scaled_min"], stats["valid_scaled_max"] = min(limits), max(limits)
    if class_counts is not None:
        stats["observed_valid_code_counts"] = dict(sorted(class_counts.items(), key=lambda item: float(item[0])))


def _prepare_one(source: Path, folder: Path, number: int, input_root: Path,
                 output_root: Path) -> dict[str, Any]:
    relative = source.relative_to(folder)
    destination = output_root / f"attachment{number:02d}" / "rasters" / relative
    validity_path = destination.with_name(destination.stem + ".analysis_valid_mask.tif")
    destination.parent.mkdir(parents=True, exist_ok=True)
    source_hash = _sha256(source)
    with rasterio.open(source) as src:
        source_metadata = _metadata(src)
        source_dependencies = []
        for filename in src.files:
            dependency = Path(filename).resolve()
            if dependency != source and dependency.is_file():
                source_dependencies.append({
                    "path": dependency.relative_to(input_root).as_posix(),
                    "sha256": _sha256(dependency),
                    "role": "source GDAL companion; archived separately by caller",
                })
        flags = _metadata_flags(source_metadata, number)


        with rasterio.Env(GDAL_TIFF_INTERNAL_MASK=True):
            copy_raster(src, destination, driver="GTiff", strict=True,
                        compress="DEFLATE", tiled=True, blockxsize=256, blockysize=256)
        with rasterio.open(destination, "r+") as dst:
            _copy_metadata(src, dst)
        mask_profile = {
            "driver": "GTiff", "width": src.width, "height": src.height,
            "count": src.count, "dtype": "uint8", "crs": src.crs,
            "transform": src.transform, "nodata": 0,
            "compress": "DEFLATE", "tiled": True, "blockxsize": 256, "blockysize": 256,
        }
        band_stats = [_new_band_stats(band, src) for band in src.indexes]
        class_counts = [Counter() if number == 5 and np.issubdtype(np.dtype(src.dtypes[band - 1]), np.integer) else None for band in src.indexes]
        checked_windows = 0
        with rasterio.open(destination) as dst, rasterio.open(validity_path, "w", **mask_profile) as mask:
            copied_metadata = _metadata(dst)


            compare_keys = [key for key in source_metadata if key != "mask_flags"]
            mismatch = [key for key in compare_keys if source_metadata[key] != copied_metadata[key]]
            if mismatch:
                raise ValueError(f"Raster metadata readback differs ({mismatch}): {source}")
            mask.update_tags(
                mask_meaning="0=invalid/missing; 1=analysis-valid",
                validity_rule="source_band_mask_nonzero AND isfinite(source_raw_value)",
                source_raster=destination.relative_to(output_root).as_posix(),
                source_sha256=source_hash,
                numeric_policy="Raw source tokens, nodata, scales and offsets are preserved in source_raster; no imputation or inferred sentinel.",
            )
            mask.descriptions = tuple(f"analysis_valid_band_{band}" for band in src.indexes)
            for band in src.indexes:
                for window in _windows(src.width, src.height):
                    values = src.read(band, window=window)
                    readback = dst.read(band, window=window)
                    source_mask = src.read_masks(band, window=window)
                    readback_mask = dst.read_masks(band, window=window)
                    if not np.array_equal(values, readback, equal_nan=True):
                        raise ValueError(f"Raster values differ on readback: {source}, band {band}, {window}")
                    if not np.array_equal(source_mask, readback_mask):
                        raise ValueError(f"Source mask differs on readback: {source}, band {band}, {window}")
                    valid = (source_mask != 0) & np.isfinite(values)
                    mask.write(valid.astype("uint8"), band, window=window)
                    _update_stats(band_stats[band - 1], values, source_mask, valid, class_counts[band - 1])
                    checked_windows += 1


        with rasterio.open(validity_path) as mask:
            if mask.crs != src.crs or mask.transform != src.transform or mask.shape != src.shape or mask.count != src.count:
                raise ValueError(f"Analysis mask grid differs: {validity_path}")
            for band in src.indexes:
                for window in _windows(src.width, src.height):
                    values = src.read(band, window=window)
                    expected = ((src.read_masks(band, window=window) != 0) & np.isfinite(values)).astype("uint8")
                    if not np.array_equal(mask.read(band, window=window), expected):
                        raise ValueError(f"Analysis mask readback differs: {validity_path}, band {band}, {window}")
        for stats, counts in zip(band_stats, class_counts):
            _finalize_stats(stats, counts)
            if stats["nonfinite_source_unmasked_count"]:
                flags.append(_flag("nonfinite_unmasked_samples", f"Band {stats['band']}: {stats['nonfinite_source_unmasked_count']} source-unmasked NaN/Inf samples are excluded by analysis mask only.", "warning"))
            if not stats["analysis_valid_count"]:
                flags.append(_flag("no_analysis_valid_samples", f"Band {stats['band']} has no finite unmasked samples.", "warning"))
            if "dhdt_err_UTM" in source.parts and stats["valid_negative_count"]:
                flags.append(_flag("negative_error_values_retained", f"Band {stats['band']}: {stats['valid_negative_count']} finite unmasked negative error values require provider review; no values were removed.", "warning"))
    if _sha256(source) != source_hash:
        raise RuntimeError(f"Source changed during raster preparation: {source}")
    return {
        "attachment": f"attachment{number:02d}",
        "source": source.relative_to(input_root).as_posix(),
        "source_sha256": source_hash,
        "source_companion_files": source_dependencies,
        "output": destination.relative_to(output_root).as_posix(),
        "output_sha256": _sha256(destination),
        "analysis_valid_mask": validity_path.relative_to(output_root).as_posix(),
        "analysis_valid_mask_sha256": _sha256(validity_path),
        "metadata": source_metadata,
        "output_mask_flags": copied_metadata["mask_flags"],
        "grid_signature": _grid_signature(source_metadata),
        "band_statistics": band_stats,
        "flags": flags,
        "verification": {
            "metadata": "pass", "all_values": "pass", "source_masks": "pass",
            "analysis_mask_values": "pass", "analysis_mask_grid": "pass",
            "band_windows_checked": checked_windows, "maximum_window_shape": [WINDOW_SIZE, WINDOW_SIZE],
            "scope": "Exact source/output equality and mask rule; this is not independent scientific validation.",
        },
    }


def _attachment_folder(input_root: Path, number: int) -> Path:
    pattern = re.compile(rf"^[（(]0?{number}[）)]")
    found = [path for path in input_root.iterdir() if path.is_dir() and pattern.match(path.name)]
    if len(found) != 1:
        raise ValueError(f"Expected exactly one attachment {number} directory under {input_root}; found {len(found)}")
    return found[0]


def clean_rasters(input_root: Path, output_root: Path) -> dict[str, Any]:


    input_root, output_root = Path(input_root).resolve(), Path(output_root).resolve()
    if not input_root.is_dir() or not output_root.is_dir():
        raise ValueError("Input and output roots must already exist as directories")
    if input_root == output_root or input_root in output_root.parents:
        raise ValueError("Output root must be outside the source attachment tree")
    folders = {number: _attachment_folder(input_root, number) for number in ATTACHMENTS}
    sources = {number: sorted(path for path in folder.rglob("*") if path.is_file() and path.suffix.lower() in {".tif", ".tiff"}) for number, folder in folders.items()}
    inventory_json = output_root / "quality" / "raster_inventory.json"
    inventory_csv = output_root / "quality" / "raster_inventory.csv"
    targets = [inventory_json, inventory_csv]
    for number, paths in sources.items():
        if not paths:
            raise ValueError(f"No source TIFFs found in attachment {number}: {folders[number]}")
        for source in paths:
            if source.is_symlink():
                raise ValueError(f"Symlink source is not accepted: {source}")
            destination = output_root / f"attachment{number:02d}" / "rasters" / source.relative_to(folders[number])
            targets.extend([destination, destination.with_name(destination.stem + ".analysis_valid_mask.tif")])
    if len(set(targets)) != len(targets):
        raise ValueError("Raster output paths collide with one another")
    for target in targets:
        if target.exists():
            raise FileExistsError(f"Refusing to overwrite raster output: {target}")


        if not target.resolve().is_relative_to(output_root):
            raise ValueError(f"Raster output escapes output root: {target}")
    inventory_json.parent.mkdir(parents=True, exist_ok=True)
    records = []
    document_notes = {}
    for number in ATTACHMENTS:
        document_notes[f"attachment{number:02d}"] = _source_document_notes(folders[number], number)
        for source in sources[number]:
            records.append(_prepare_one(source, folders[number], number, input_root, output_root))
    grid_groups: dict[str, list[str]] = {}
    for record in records:
        grid_groups.setdefault(record["grid_signature"], []).append(record["output"])
    inventory = {
        "schema_version": "1.0", "scope": "basic_lossless_raster_preparation",
        "policy": {
            "raw_values": "preserved per band, including finite source NoData tokens and nonfinite values",
            "source_masks": "preserved and fully compared after writing",
            "analysis_valid_mask": "uint8 per band: 1 iff source mask is nonzero and source raw numeric value is finite; 0 means missing/invalid",
            "spatial_operations": "none: no reprojection, resampling, mosaic, crop or forced alignment",
            "numeric_operations": "none: no imputation, clipping, inferred sentinel, unit conversion or class remapping",
            "statistics": "computed only on finite source-unmasked values; raw and scale/offset-derived ranges are separate",
            "tags": "source default/band and substantive namespace tags preserved; layout and generated DERIVED_SUBDATASETS domains describe new files and are not copied",
            "grid_signatures": "exact CRS WKT, affine, dimensions and pixel registration; different signatures must not be multiplied by array index without a separate alignment decision; tiny floating differences are retained",
            "metadata_pending": "missing units, class definitions and vertical datum remain pending locally; numerical preservation may pass independently",
        },
        "counts": {f"attachment{number:02d}": len(sources[number]) for number in ATTACHMENTS},
        "provider_document_notes": document_notes,
        "grid_groups": grid_groups,
        "records": records,
    }
    _dump_json(inventory_json, inventory)
    with inventory_csv.open("x", encoding="utf-8-sig", newline="") as handle:
        columns = ["attachment", "source", "output", "analysis_valid_mask", "grid_signature", "crs", "width", "height", "count", "dtype", "resolution", "bounds", "affine", "nodata", "units", "analysis_valid_count", "source_masked_count", "nonfinite_count", "flags", "readback_verification"]
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for record in records:
            metadata = record["metadata"]
            row = {key: record[key] for key in ("attachment", "source", "output", "analysis_valid_mask", "grid_signature")}
            row.update({key: metadata[key] for key in ("crs", "width", "height", "count", "resolution", "bounds", "nodata", "units")})
            row.update(dtype=metadata["dtypes"], affine=metadata["transform"], flags=record["flags"], readback_verification="pass")
            for key in ("analysis_valid_count", "source_masked_count", "nonfinite_count"):
                row[key] = [band[key] for band in record["band_statistics"]]
            writer.writerow({key: json.dumps(_json_value(value), ensure_ascii=False, allow_nan=False) if isinstance(value, (dict, list)) else value for key, value in row.items()})
    outputs = [record[key] for record in records for key in ("output", "analysis_valid_mask")]
    outputs.extend([inventory_json.relative_to(output_root).as_posix(), inventory_csv.relative_to(output_root).as_posix()])
    summary = {
        "status": "completed_with_metadata_pending" if any(flag["severity"] == "pending" for record in records for flag in record["flags"]) else "completed",
        "raster_count": len(records), "analysis_mask_count": len(records),
        "counts_by_attachment": inventory["counts"], "grid_signature_count": len(grid_groups),
        "readback_verification": "pass", "independent_scientific_validation": "not_performed",
        "pending_metadata_flag_count": sum(flag["severity"] == "pending" for record in records for flag in record["flags"]),
        "warning_flag_count": sum(flag["severity"] == "warning" for record in records for flag in record["flags"]),
        "inventory_json": inventory_json.relative_to(output_root).as_posix(),
        "inventory_csv": inventory_csv.relative_to(output_root).as_posix(),
        "outputs": outputs,
    }

    summary.update(
        primary_files=len(records),
        band_cells=sum(band["sample_count"] for record in records for band in record["band_statistics"]),
        analysis_valid_cells=sum(band["analysis_valid_count"] for record in records for band in record["band_statistics"]),
        output_paths=outputs,
        inventory_path=inventory_json.relative_to(output_root).as_posix(),
    )
    return summary
