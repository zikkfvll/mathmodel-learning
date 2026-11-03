from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

BASE = Path(__file__).resolve().parent
ENDPOINT = "https://archive-api.open-meteo.com/v1/archive"
DOC_URL = "https://open-meteo.com/en/docs/historical-weather-api"
ZONE = dt.timezone(dt.timedelta(hours=8))


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fetch(url: str) -> tuple[bytes, dict]:
    request = Request(url, headers={"User-Agent": "MathorCup-Q1-research/1.0"})
    with urlopen(request, timeout=55) as response:
        body = response.read()
        metadata = {"http_status": response.status,
                    "retrieved_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                    "response_headers": dict(response.headers)}
    return body, metadata


def parse_year(year: int, obj: dict) -> tuple[list[dict], dict]:
    if obj.get("error"):
        raise ValueError(f"API error: {obj}")
    assert obj["utc_offset_seconds"] == 28800
    assert obj["timezone"] == "Asia/Shanghai"
    assert obj["hourly_units"]["precipitation"] == "mm"
    assert obj["hourly_units"]["rain"] == "mm"
    times = obj["hourly"]["time"]
    precip = obj["hourly"]["precipitation"]
    rain = obj["hourly"]["rain"]
    assert len(times) == len(precip) == len(rain)
    start = dt.datetime(year, 5, 1, tzinfo=ZONE)
    last = dt.datetime(year, 10, 16, 23, tzinfo=ZONE)
    expected = int((last - start).total_seconds() / 3600) + 1
    assert len(times) == expected, (year, len(times), expected)
    rows = []
    for i, (label, p, r) in enumerate(zip(times, precip, rain)):
        end = dt.datetime.fromisoformat(label).replace(tzinfo=ZONE)
        assert end == start + dt.timedelta(hours=i), (year, i, label)
        for value in (p, r):
            assert value is not None and math.isfinite(value) and value >= 0, (label, p, r)

        rows.append({"interval_start": (end-dt.timedelta(hours=1)).isoformat(),
                     "interval_end": end.isoformat(),
                     "precipitation_mm": p, "rain_mm": r,
                     "source_product": "Open-Meteo archive / ERA5",
                     "source_qc_status": "not_provided",
                     "grid_latitude": obj["latitude"],
                     "grid_longitude": obj["longitude"],
                     "grid_elevation_m": obj["elevation"]})
    checks = {"year": year, "status": "pass", "requested_first_label": start.isoformat(),
              "requested_last_label": last.isoformat(), "expected_count": expected,
              "actual_count": len(rows), "missing_precipitation": 0, "missing_rain": 0,
              "nonfinite_values": 0, "negative_values": 0, "duplicate_labels": 0,
              "hourly_continuity": True, "sum_precipitation_mm": round(sum(precip), 10),
              "sum_rain_mm": round(sum(rain), 10),
              "max_precipitation_mm": max(precip), "wet_hours_gt_zero": sum(p>0 for p in precip),
              "zero_hours": sum(p==0 for p in precip),
              "returned_grid": {k:obj[k] for k in ["latitude","longitude","elevation"]},
              "timezone": obj["timezone"], "utc_offset_seconds": obj["utc_offset_seconds"],
              "hourly_units": obj["hourly_units"]}
    return rows, checks


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true", help="refetch existing raw years")
    args = parser.parse_args()
    BASE.mkdir(exist_ok=True, parents=True)
    rows, checks, requests = [], [], []
    for year in range(2021, 2027):
        params = dict(latitude=28.86, longitude=85.29, start_date=f"{year}-05-01",
                      end_date=f"{year}-10-16", hourly="precipitation,rain", models="era5",
                      cell_selection="nearest", elevation="nan", timezone="Asia/Shanghai",
                      precipitation_unit="mm", timeformat="iso8601")
        url = ENDPOINT + "?" + urlencode(params)
        raw = BASE / f"era5_hourly_{year}.json"
        receipt = BASE / f"request_{year}.json"
        if args.refresh or not raw.exists():
            body, metadata = fetch(url)
            raw.write_bytes(body)
            receipt.write_text(json.dumps({"url":url,"params":params,**metadata},
                                          ensure_ascii=False, indent=2), encoding="utf-8")
        request_record = json.loads(receipt.read_text(encoding="utf-8"))
        assert request_record["url"] == url
        year_rows, year_checks = parse_year(year, json.loads(raw.read_text(encoding="utf-8")))
        year_checks["raw_sha256"] = sha(raw)
        rows.extend(year_rows)
        checks.append(year_checks)
        requests.append(request_record)
        print(f"{year}: {len(year_rows)} hourly rows checked", flush=True)
    grids = {(c["returned_grid"]["latitude"], c["returned_grid"]["longitude"],
              c["returned_grid"]["elevation"]) for c in checks}
    assert len(grids) == 1, grids
    out = BASE / "era5_precipitation_hourly_2021_2026.csv"
    with out.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    report = {"status":"pass", "script_sha256":sha(Path(__file__)),
              "csv_sha256":sha(out), "rows":len(rows), "same_product_all_years":True,
              "same_returned_grid_all_years":True, "year_checks":checks,
              "source_qc_status":"not_provided", "requests":requests,
              "support":"Each API label t is preceding-hour sum, assigned [t-1h,t).",
              "role":"D1 independent grid proxy; never inserted into D0 station missing values."}
    (BASE / "DOWNLOAD_CHECKS.json").write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
    doc = BASE / "open_meteo_historical_weather_api.html"
    if args.refresh or not doc.exists():
        body, metadata = fetch(DOC_URL)
        doc.write_bytes(body)
        (BASE / "documentation_request.json").write_text(
            json.dumps({"url":DOC_URL, **metadata},ensure_ascii=False,indent=2),encoding="utf-8")


if __name__ == "__main__":
    main()
