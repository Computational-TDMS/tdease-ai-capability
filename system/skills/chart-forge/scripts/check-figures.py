# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""Self-containment check for figures/reports handed to the user.

Verifies what chart-forge's build.py and the toppicr report guarantee by
construction, so an agent does not re-derive the checks by hand per figure:

  1. no unreplaced template tokens (/*__Xxx__*/ leftovers);
  2. no external network references (http(s) in src/href/url(); data: URIs are fine);
  3. a data payload exists (window.DATA / const DATA assignment).

One bounded JSON verdict on stdout — per file and overall. Directories are
scanned for *.html (non-recursive). Exit code 1 when any file fails.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

TOKEN_LEFTOVER = re.compile(r"/\*__[A-Za-z0-9]+__\*/")
EXTERNAL_REF = re.compile(r"""(?:src|href)\s*=\s*["']https?://[^"']*["']""", re.IGNORECASE)
EXTERNAL_FETCH = re.compile(r"""(?:fetch|import)\s*\(\s*["']https?://""")
DATA_PAYLOAD = re.compile(r"(?:window\.DATA\s*=|(?:const|let|var)\s+DATA\s*=)")


def check_file(path: Path) -> dict:
    html = path.read_text(encoding="utf-8", errors="replace")
    external = sorted(set(m.group(0) for m in EXTERNAL_REF.finditer(html)))
    external += sorted(set(m.group(0) for m in EXTERNAL_FETCH.finditer(html)))
    verdict = {
        "file": str(path),
        "bytes": len(html),
        "unreplacedTokens": sorted(set(TOKEN_LEFTOVER.findall(html))),
        "externalRefs": external,
        "dataPayload": bool(DATA_PAYLOAD.search(html)),
    }
    verdict["ok"] = not verdict["unreplacedTokens"] and not external and verdict["dataPayload"]
    return verdict


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("targets", nargs="+",
                        help="HTML files, or directories to scan for *.html")
    args = parser.parse_args()

    files: list[Path] = []
    for target in args.targets:
        path = Path(target)
        if path.is_dir():
            files.extend(sorted(path.glob("*.html")))
        else:
            files.append(path)
    if not files:
        raise SystemExit(f"no HTML files found under: {' '.join(args.targets)}")

    results = [check_file(path) for path in files]
    payload = {
        "checked": len(results),
        "ok": all(result["ok"] for result in results),
        "results": results,
    }
    print(json.dumps(payload, separators=(",", ":")))
    return 0 if payload["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
