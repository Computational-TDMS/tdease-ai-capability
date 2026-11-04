#!/usr/bin/env python3
"""Generate mgf-accumulated/data.json from an MGF exported by
`tdease-ms spectrum --mode scans --scans ... --out-format mgf`.

Pass-through adapter: peaks are kept verbatim (m/z 4 dp, intensity 3 dp) — the
viewer re-accumulates the checked scans live with the quantify.py
_sum_peak_lists contract (ppm chain-merge, intensity-weighted m/z), so scan
selection stays interactive. TITLE's `activation=` (tdease-ms ≥ 2026-08-30)
groups ETD vs HCD so complementary series are never blindly summed.

Default --mgf is the bundled sample targeted family: precursor 631.749
z = 10+, 6 ETD + 5 HCD scans."""
import argparse
import json
import re
from pathlib import Path

ACTIVATION_SHORT = {
    "electron transfer dissociation": "ETD",
    "beam-type collision-induced dissociation": "HCD",
    "higher-energy collisional dissociation": "HCD",
    "collision-induced dissociation": "CID",
    "electron capture dissociation": "ECD",
}


def parse_mgf(path: Path) -> list[dict]:
    scans, cur = [], None
    for line in path.read_text(encoding="utf-8").splitlines():
        t = line.strip()
        if t.upper() == "BEGIN IONS":
            cur = {"scan": None, "rt": None, "activation": None,
                   "pepmass": None, "charge": None, "peaks": []}
        elif t.upper() == "END IONS":
            if cur and cur["peaks"]:
                scans.append(cur)
            cur = None
        elif cur is None:
            continue
        elif t.startswith("TITLE="):
            m = re.search(r"scan=(\d+)", t)
            if m:
                cur["scan"] = int(m.group(1))
            if "activation=" in t:
                cur["activation"] = t.split("activation=", 1)[1].strip()
        elif t.startswith("PEPMASS="):
            cur["pepmass"] = float(t.split("=", 1)[1].split()[0])
        elif t.startswith("CHARGE="):
            raw = t.split("=", 1)[1]
            sign = -1 if raw.rstrip().endswith("-") else 1
            cur["charge"] = sign * int(raw.rstrip("+-"))
        elif t.startswith("RTINSECONDS="):
            cur["rt"] = round(float(t.split("=", 1)[1]) / 60.0, 2)
        else:
            parts = t.split()
            if len(parts) == 2:
                cur["peaks"].append([round(float(parts[0]), 4), round(float(parts[1]), 3)])
    return scans


def main() -> None:
    here = Path(__file__).parent
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mgf", type=Path, default=here / "scans-631z10.mgf")
    ap.add_argument("--tolerance-ppm", type=float, default=10.0,
                    help="chain-merge tolerance the viewer uses for accumulation")
    ap.add_argument("--out", type=Path, default=here / "data.json")
    args = ap.parse_args()

    scans = parse_mgf(args.mgf)
    if not scans:
        raise SystemExit(f"no BEGIN IONS blocks with peaks in {args.mgf}")

    pepmass = next((s["pepmass"] for s in scans if s["pepmass"] is not None), None)
    charge = next((s["charge"] for s in scans if s["charge"] is not None), None)
    for s in scans:
        s["activation"] = ACTIVATION_SHORT.get(s["activation"] or "", s["activation"] or "?")

    groups = {}
    for s in scans:
        groups.setdefault(s["activation"], []).append(s["scan"])
    tally = " + ".join(f"{len(v)} {k}" for k, v in sorted(groups.items()))

    # annotations fixture: flag the strongest peak and its isotope
    # neighbourhood (dynamic — survives data refreshes)
    agg = {}
    for s in scans:
        for mz, inten in s["peaks"]:
            agg[mz] = agg.get(mz, 0.0) + inten
    strongest = max(agg, key=agg.get)

    payload = {
        "data": {
            "precursor": {"mz": pepmass, "charge": charge},
            "tolerancePpm": args.tolerance_ppm,
            "scans": [{"scan": s["scan"], "rt": s["rt"], "activation": s["activation"],
                       "peaks": s["peaks"]} for s in scans],
            "annotations": [
                {"kind": "callout", "mz": strongest, "label": "strongest peak"},
                {"kind": "band", "from": round(strongest - 8, 4), "to": round(strongest + 8, 4),
                 "label": "isotope cluster"},
            ],
        },
        "meta": {
            "title": f"Accumulated MS2 — precursor {pepmass:.4f} z={abs(charge)}{'+' if charge >= 0 else '-'}",
            "subtitle": (f"{tally} scans · grey = individual scans (display-capped), accent = accumulated "
                         f"at {args.tolerance_ppm:g} ppm chain-merge, intensity-weighted m/z · "
                         "ETD/HCD groups toggle independently"),
            "source": f"tdease-ms spectrum --mode scans --out-format mgf · {args.mgf.name}",
            "transform": "MGF parsed verbatim (scan/RT/activation/peaks); accumulation runs live in the viewer over the checked scans",
            "thresholds": {"tolerance_ppm": args.tolerance_ppm},
            "n": f"{len(scans)} scans ({tally})",
        },
    }
    args.out.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
    n_peaks = sum(len(s["peaks"]) for s in scans)
    print(f"{len(scans)} scans ({tally}), {n_peaks} peaks -> {args.out} ({args.out.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
