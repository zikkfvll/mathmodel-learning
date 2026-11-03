from pathlib import Path
import csv
import gzip
import io
import json
import shutil

from netCDF4 import Dataset
import numpy as np
from pyproj import CRS


def _jsonable(value):
    if isinstance(value, dict):
        return {k: _jsonable(v) for k, v in value.items()}
    if isinstance(value, (tuple, list, np.ndarray)):
        return [_jsonable(v) for v in value]
    if isinstance(value, np.generic):
        return _jsonable(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return str(value)
    return value


def _decoded(variable):
    variable.set_auto_maskandscale(True)
    masked = np.ma.asarray(variable[:], dtype=np.float64)
    values = np.asarray(masked.filled(np.nan))
    source_missing = np.ma.getmaskarray(masked) | np.isnan(values)
    nonfinite = ~np.isfinite(values) & ~source_missing
    return values, source_missing, nonfinite


def _physical(name, units, values):

    failed = np.zeros(values.shape, dtype=bool)
    rule = "none; retain_without_inventing_bounds"
    if name == "slope" and units == "degree":
        failed = (values < 0) | (values > 90)
        rule = "terrain_slope_degrees_in_0_90"
    elif name == "thickness" and units == "m":
        failed = values < 0
        rule = "ice_thickness_metres_nonnegative"
    elif name in {"v", "v_amp", "vx_amp", "vy_amp"} and units == "meter/year":
        failed = values < 0
        rule = "speed_or_amplitude_nonnegative; signed_vx_vy_exempt"
    return failed & np.isfinite(values), rule


def clean_netcdf(input_root: Path, output_root: Path) -> dict:
    folders = [p for p in input_root.iterdir() if p.is_dir() and p.name.startswith(("（4）", "(4)"))]
    if len(folders) != 1:
        raise ValueError("Expected one attachment 4 directory")
    files = sorted(p for p in folders[0].rglob("*.nc") if not p.name.startswith("._"))
    if not files:
        raise ValueError("No attachment 4 NetCDF files")
    destination = output_root / "attachment04"
    destination.mkdir(parents=True, exist_ok=False)
    quality = output_root / "quality"
    quality.mkdir(exist_ok=True)
    reports = []
    outputs = []
    for source in files:
        target = destination / source.name
        table = destination / (source.stem + "_grid_cells.csv.gz")
        flags_file = destination / (source.stem + "_cell_quality.csv.gz")
        shutil.copy2(source, target)
        grids, statuses, variable_reports = {}, {}, []
        updates = {}
        with Dataset(source) as dataset:
            dimensions = {name: dict(size=len(dim), unlimited=dim.isunlimited())
                          for name, dim in dataset.dimensions.items()}
            global_attrs = {name: dataset.getncattr(name) for name in dataset.ncattrs()}
            if not {"x", "y", "spatial_ref"}.issubset(dataset.variables):
                raise ValueError(f"No explicit x/y/CRS in {source.name}; cannot create a spatial table")
            x, xm, xn = _decoded(dataset["x"])
            y, ym, yn = _decoded(dataset["y"])
            if xm.any() or ym.any() or xn.any() or yn.any():
                raise ValueError("Invalid coordinate labels; raw NetCDF is retained, coordinate alignment requires review")
            reference = dataset["spatial_ref"]
            wkt = getattr(reference, "crs_wkt", getattr(reference, "spatial_ref", None))
            crs = CRS.from_wkt(wkt) if wkt else None
            transform_text = getattr(reference, "GeoTransform", "")
            transform = [float(v) for v in transform_text.split()] if transform_text else None
            coordinate_check = "not_declared"
            if transform:
                if len(transform) != 6:
                    raise ValueError("Invalid source GeoTransform")
                a, b, c, d, e, f = transform
                if c or e:
                    raise ValueError("Rotated transform cannot be inferred from 1D coordinates")
                expected_x = a + (np.arange(len(x)) + 0.5) * b
                expected_y = d + (np.arange(len(y)) + 0.5) * f
                if not np.allclose(x, expected_x, rtol=0, atol=1e-7) or not np.allclose(y, expected_y, rtol=0, atol=1e-7):
                    raise ValueError("Source coordinates conflict with source transform; no automatic shift")
                coordinate_check = "source_coordinates_equal_transform_pixel_centres"
            if len(np.unique(x)) != len(x) or len(np.unique(y)) != len(y):
                raise ValueError("Duplicate coordinates require an explicit reconciliation policy")
            for name, variable in dataset.variables.items():
                if variable.dtype.kind not in "iuf":
                    raise ValueError(f"Unsupported variable type (not silently skipped): {name}")
                values, missing, nonfinite = _decoded(variable)
                units = getattr(variable, "units", "not_declared")
                failed, rule = _physical(name, units, values)
                status = np.ones(values.shape, dtype=np.uint8)
                status[missing] = 0
                status[nonfinite] = 2
                status[failed] = 3
                analysis = values.copy()
                analysis[status != 1] = np.nan
                if np.any(nonfinite | failed):
                    updates[name] = nonfinite | failed
                if variable.dimensions == ("y", "x"):
                    grids[name], statuses[name] = analysis, status
                valid = analysis[np.isfinite(analysis)]
                variable_reports.append(dict(name=name, dimensions=list(variable.dimensions), shape=list(values.shape),
                    source_dtype=str(variable.dtype), units=units,
                    attrs={k: variable.getncattr(k) for k in variable.ncattrs()},
                    source_missing=int(missing.sum()), nonfinite_excluded=int(nonfinite.sum()),
                    physical_excluded=int(failed.sum()), retained=int(np.isfinite(analysis).sum()),
                    min=float(valid.min()) if valid.size else None, max=float(valid.max()) if valid.size else None,
                    physical_rule=rule, source_qc="not_provided"))
            velocity_check = None
            if {"v", "vx", "vy"}.issubset(grids):
                common = np.isfinite(grids["v"]) & np.isfinite(grids["vx"]) & np.isfinite(grids["vy"])
                discrepancy = np.abs(grids["v"][common] - np.hypot(grids["vx"][common], grids["vy"][common]))
                velocity_check = dict(paired_cells=int(common.sum()),
                    max_absolute_difference=float(discrepancy.max()) if discrepancy.size else None,
                    action="observation_only; bilinear-source-processing-can-change-identity; no values altered")
        if updates:
            with Dataset(target, "r+") as dataset:
                for name, exclusion in updates.items():
                    variable = dataset[name]
                    if "_FillValue" not in variable.ncattrs():
                        raise ValueError(f"No declared fill value for excluding {name}; cannot invent one")
                    variable.set_auto_maskandscale(False)
                    raw = np.asarray(variable[:]).copy()
                    raw[exclusion] = variable.getncattr("_FillValue")
                    variable[:] = raw

        with Dataset(source) as original, Dataset(target) as cleaned:
            if original.ncattrs() != cleaned.ncattrs() or set(original.variables) != set(cleaned.variables):
                raise RuntimeError("NetCDF metadata/schema changed unexpectedly")
            for key in original.ncattrs():
                if _jsonable(original.getncattr(key)) != _jsonable(cleaned.getncattr(key)):
                    raise RuntimeError(f"Changed global attribute: {key}")
            for name, variable in original.variables.items():
                copy = cleaned[name]
                if variable.dimensions != copy.dimensions or variable.dtype != copy.dtype:
                    raise RuntimeError(f"Changed variable schema: {name}")
                if _jsonable({k: variable.getncattr(k) for k in variable.ncattrs()}) != _jsonable({k: copy.getncattr(k) for k in copy.ncattrs()}):
                    raise RuntimeError(f"Changed variable attributes: {name}")
                variable.set_auto_maskandscale(False)
                copy.set_auto_maskandscale(False)
                raw, new = np.asarray(variable[:]), np.asarray(copy[:])
                keep = ~updates[name] if name in updates else np.ones(raw.shape, dtype=bool)
                if not np.array_equal(raw[keep], new[keep], equal_nan=True):
                    raise RuntimeError(f"Changed retained raw values: {name}")
                if name in grids:
                    decoded, missing, nonfinite = _decoded(copy)
                    if not np.array_equal(np.where(np.isfinite(decoded), decoded, np.nan), grids[name], equal_nan=True):
                        raise RuntimeError(f"NetCDF/table analysis value mismatch: {name}")
        for filename, fields, flag_table in [(table, ["source_row", "source_column", "x", "y", *grids], False),
                                              (flags_file, ["source_row", "source_column", "x", "y", "variable", "status"], True)]:
            with filename.open("wb") as binary, gzip.GzipFile(filename="", fileobj=binary, mode="wb", mtime=0) as gz:
                with io.TextIOWrapper(gz, encoding="utf-8", newline="") as text:
                    writer = csv.DictWriter(text, fieldnames=fields)
                    writer.writeheader()
                    for row in range(len(y)):
                        for col in range(len(x)):
                            location = dict(source_row=row + 1, source_column=col + 1, x=float(x[col]), y=float(y[row]))
                            if flag_table:
                                for name, status in statuses.items():
                                    if status[row, col] != 1:
                                        writer.writerow(dict(location, variable=name, status={0: "source_missing", 2: "excluded_nonfinite", 3: "excluded_physical"}[int(status[row, col])]))
                            else:
                                writer.writerow(dict(location, **{name: float(values[row, col]) if np.isfinite(values[row, col]) else None
                                                                  for name, values in grids.items()}))
        report = dict(source_path=source.relative_to(input_root).as_posix(),
            output_path=target.relative_to(output_root).as_posix(), grid_table=table.relative_to(output_root).as_posix(),
            quality_table=flags_file.relative_to(output_root).as_posix(), dimensions=dimensions, global_attrs=global_attrs,
            crs=crs.to_string() if crs else None, crs_wkt=wkt, geotransform=transform,
            coordinate_check=coordinate_check, has_time_dimension="time" in dimensions,
            time_note="No time dimension: calendar/time decoding not applicable; velocities are source-described 2014–2024 climatology.",
            variables=variable_reports, grid_rows=len(x) * len(y), velocity_consistency_observation=velocity_check,
            vertical_datum="not_declared; do not assume compatibility with attachment07 DEM", decoded_once=True,
            nonfinite_or_physical_exclusions=sum(int(m.sum()) for m in updates.values()),
            no_resampling=True, no_imputation=True)
        reports.append(_jsonable(report))
        outputs.extend([p.relative_to(output_root).as_posix() for p in [target, table, flags_file]])
        print(f"NetCDF {source.name}: {len(x) * len(y):,} grid cells; coordinates/units/masks retained", flush=True)
    (quality / "netcdf_inventory.json").write_text(json.dumps(reports, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return dict(files=len(files), grid_rows=sum(r["grid_rows"] for r in reports), outputs=outputs,
                exclusions=sum(r["nonfinite_or_physical_exclusions"] for r in reports),
                inventory="quality/netcdf_inventory.json")
