import argparse
import csv
import datetime as dt
import hashlib
import importlib.metadata
import json
from pathlib import Path
import shutil
import sys

import clean_data
from netcdf_cleaning import clean_netcdf
from raster_cleaning import clean_rasters
from vector_cleaning import clean_vectors


BASE = Path(__file__).resolve().parent
CODE_FILES = ["clean_all.py", "clean_data.py", "quality_policy.py", "quality_diagnostics.py",
              "quality_rules.json", "netcdf_cleaning.py", "raster_cleaning.py", "vector_cleaning.py"]


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def inventory(input_root):
    folders = {}
    for number in range(1, 9):
        matches = [p for p in input_root.iterdir() if p.is_dir() and p.name.startswith((f"（{number}）", f"({number})"))]
        if len(matches) != 1:
            raise ValueError(f"Expected one folder for attachment {number}, found {len(matches)}")
        folders[number] = matches[0]
    records = []
    for number, folder in folders.items():
        for path in sorted(p for p in folder.rglob("*") if p.is_file()):
            suffix = path.suffix.lower()
            primary = {1: {".csv": "station_csv"}, 2: {".xlsx": "daily_xlsx"},
                       3: {".tif": "glacier_raster", ".shp": "glacier_vector"},
                       4: {".nc": "glacier_netcdf"}, 5: {".tif": "landcover_raster"},
                       6: {".shp": "historical_disaster_vector"}, 7: {".tif": "dem_raster"},
                       8: {".shp": "population_vector"}}
            if path.name.startswith(("._", ".DS_Store", "~$")) or "__MACOSX" in path.parts:
                role, handling = "operating_system_metadata", "listed_only_not_data"
            elif suffix == ".lock":
                role, handling = "software_lock", "listed_only_not_data"
            elif suffix in primary[number]:
                role, handling = primary[number][suffix], "cleaned"
            elif suffix in {".doc", ".docx", ".pdf", ".jpg", ".png"}:
                role, handling = "reference_document_or_image", "copied_unchanged_to_support"
            elif suffix in {".dbf", ".shx", ".prj", ".cpg", ".sbn", ".sbx", ".xml", ".tfw", ".ovr"}:
                role, handling = "source_companion", "copied_unchanged_to_support"
            else:
                raise ValueError(f"Unclassified attachment file (not silently skipped): {path}")
            archive = handling == "copied_unchanged_to_support" or suffix == ".shp"
            support = (Path("support") / f"attachment{number:02d}" / path.relative_to(folder)).as_posix() if archive else ""
            records.append(dict(attachment=number, input_path=path.relative_to(input_root).as_posix(),
                                bytes=path.stat().st_size, sha256=sha(path), role=role, handling=handling,
                                support_copy=support))
    return records


def main(input_root, output_dir):
    input_root, out = input_root.resolve(), output_dir.resolve()
    if not input_root.is_dir():
        raise ValueError("--input-root must be the official 附件 directory")
    if out.exists():
        raise ValueError("Output exists; choose a new directory to preserve previous deliveries")
    if input_root == out or input_root in out.parents:
        raise ValueError("Output must be outside the original attachment directory")
    started = dt.datetime.now(dt.timezone.utc).isoformat()
    sources = inventory(input_root)
    code_before = {name: sha(BASE / name) for name in CODE_FILES}
    out.mkdir(parents=True)
    (out / "quality").mkdir()
    with (out / "input_inventory_all.csv").open("w", encoding="utf-8-sig", newline="") as text:
        writer = csv.DictWriter(text, fieldnames=list(sources[0]))
        writer.writeheader()
        writer.writerows(sources)

    clean_data.main(input_root, out / "attachment01_02")
    tabular = json.loads((out / "attachment01_02/run_summary.json").read_text(encoding="utf-8"))
    rasters = clean_rasters(input_root, out)
    vectors = clean_vectors(input_root, out)
    netcdf = clean_netcdf(input_root, out)
    for source in sources:
        if source["support_copy"]:
            target = out / source["support_copy"]
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(input_root / source["input_path"], target)
            if sha(target) != source["sha256"]:
                raise RuntimeError(f"Support file changed: {target.name}")
    if any(sha(input_root / source["input_path"]) != source["sha256"] for source in sources):
        raise RuntimeError("An original attachment changed during cleaning")
    if any(sha(BASE / name) != digest for name, digest in code_before.items()):
        raise RuntimeError("Cleaner code changed during execution")
    versions = {name: importlib.metadata.version(name) for name in [
        "numpy", "pandas", "openpyxl", "rasterio", "geopandas", "pyogrio", "pyproj", "shapely", "netCDF4"]}
    summary = dict(scope="all eight official attachments; basic cleaning only", started_utc=started,
        finished_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
        input_files=len(sources), primary_files=sum(s["handling"] == "cleaned" for s in sources),
        copied_support_files=sum(bool(s["support_copy"]) for s in sources),
        ignored_non_data_files=[s["input_path"] for s in sources if s["handling"] == "listed_only_not_data"],
        counts_by_attachment={str(i): dict(
            files=sum(s["attachment"] == i for s in sources),
            primary_files=sum(s["attachment"] == i and s["handling"] == "cleaned" for s in sources)) for i in range(1, 9)},
        tabular=dict(source_files=tabular["input_files"], station_rows=tabular["station_rows"], daily_rows=tabular["daily_rows"],
                     measurement_cells=tabular["cells"], retained_values=tabular["analysis_present"], statuses=tabular["statuses"]),
        rasters=rasters, vectors=vectors, netcdf=netcdf,
        source_inputs_unchanged=True, code_sha256=code_before,
        environment=dict(python=sys.version.split()[0], **versions),
        imputation=False, clipping=False, reprojection=False, resampling=False,
        model_training=False, new_observations=False,
        limitations=[
            "Station interval/aggregation and daily temperature units/day boundary remain unconfirmed.",
            "Raster native grids/CRS/NoData are preserved; different grids are not automatically aligned.",
            "Glacier dhdt units, land-cover class meanings and vertical datums remain partially unconfirmed.",
            "Original invalid vector geometries remain in all_features layers and are isolated in analysis_valid_geometry layers.",
            "Historical-disaster geometries are polygons; no centroids or artificial dates/labels are created.",
            "Attachment04 velocity fields are climatological, with no time dimension.",
            "Population attribute meanings/totals are preserved; no guessed reallocation or normalization."],
        output_sha256={p.relative_to(out).as_posix(): sha(p) for p in sorted(out.rglob("*")) if p.is_file()})
    dump(out / "run_summary_all.json", summary)
    print(json.dumps({k: summary[k] for k in ["scope", "input_files", "primary_files", "copied_support_files", "tabular"]},
                     ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True, help="Official 附件 folder with （1）…（8）")
    parser.add_argument("--output-dir", type=Path, required=True, help="New output folder outside originals")
    args = parser.parse_args()
    main(args.input_root, args.output_dir)
