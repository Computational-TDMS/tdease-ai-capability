# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""TDP OpenSearch capability command router."""

import importlib
from pathlib import Path
import sys

_scripts = Path(__file__).resolve().parent / "scripts"
if str(_scripts) not in sys.path:
    sys.path.insert(0, str(_scripts))
_skills = Path(__file__).resolve().parent.parent
if str(_skills) not in sys.path:
    sys.path.insert(0, str(_skills))

GROUPS = {
    "inspect": "census_inspect",
    "ingest": "ingest",
    "families": "family_discovery",
    "figures": "figures",
}

USAGE = f"""usage: {Path(sys.argv[0]).name} <family> [args...]

tools:   {' '.join(sorted(GROUPS))}"""


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        print(USAGE)
        return 2
    head, rest = argv[0], argv[1:]
    if head in ("-h", "--help"):
        print(USAGE)
        return 0
    module = GROUPS.get(head)
    if module is None:
        print(f"unknown family {head!r}\n\n{USAGE}", file=sys.stderr)
        return 2
    sys.argv = [sys.argv[0]] + rest
    return importlib.import_module(module).main()


if __name__ == "__main__":
    raise SystemExit(main())
