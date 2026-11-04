#!/usr/bin/env python3
"""Generate ms1-family-map/data.json from an ms1-families CSV (the
isoforma ladder table: spacing/member_count/RT window/masses per family).
Pass-through: the chart derives spacing, RT centre and MS2 support itself."""
import argparse
import csv
import json
from pathlib import Path


def main() -> None:
    here = Path(__file__).parent
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--src", type=Path, default=here / "ms1-families.csv")
    ap.add_argument("--out", type=Path, default=here / "data.json")
    args = ap.parse_args()

    with args.src.open(encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise SystemExit(f"no rows in {args.src}")

    payload = {
        "data": rows,
        "meta": {
            "source": args.src.name,
            "transform": "CSV rows passed through verbatim; chart derives RT centre = (common_rt_start + common_rt_end)/2",
            "n": len(rows),
        },
    }
    args.out.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
    print(f"{len(rows)} isotope families -> {args.out}")


if __name__ == "__main__":
    main()
