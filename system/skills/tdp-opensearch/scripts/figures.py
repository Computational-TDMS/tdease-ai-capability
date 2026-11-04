"""Census figures from OpenSearch CSV/JSON artifacts — chart-forge adapter.

Ruled 2026-08-31: this script computes and writes `{data, meta}` contracts,
then calls chart-forge `build.py --example <template>` for the self-contained
interactive HTML (the runtime's only rendering stack). No plotly, no CDN;
PNG/SVG are the figure page's client-side export buttons.

``report`` is the pipeline default: ONE self-contained HTML carrying the small
discovery view set (census + families, with focused panels opt-in). The
per-figure modes remain for special plotting paths, each with its own
parameters; they must not become the default loop of the pipeline.

Standard interface (figure modes): ``--out`` takes the deliverable path
(``--plot`` is an alias); the deliverable is one self-contained interactive
HTML plus a ``.data.json`` contract sidecar and a ``.manifest.json``; PNG/SVG
export lives in the figure page's Save buttons. Stdout is one JSON envelope
({figure, plot, written}); ``report`` emits ({report, figures, sections}).

Templates (system/skills/chart-forge/examples/):
  mass-shifts (alias: census) / shift-spectrum -> opensearch-census
      (sticks + evidence marks)
  grading-levels            -> grading-levels        (Smith grading histogram)
  featuremap              -> ms1-family-map (one row per isotope family)
  proteoform              -> proteoform-map       (packed lane bars + boxes)
"""

import argparse

from figures_report import report_main
from figures_census import census_main, shift_spectrum_main
from figures_featuremap import featuremap_main
from figures_proteoform import proteoform_main


def main():
    root = argparse.ArgumentParser(add_help=False, description=__doc__)
    root.add_argument("mode", nargs="?",
                      choices=["report", "mass-shifts", "census",
                               "shift-spectrum", "featuremap", "proteoform"])
    args, rest = root.parse_known_args()
    if args.mode is None:
        if rest and rest[0] in ("-h", "--help"):
            root.print_help()
            return 0
        root.print_usage()
        return 2
    return {
        "report": report_main,
        "mass-shifts": census_main,
        "census": census_main,  # legacy alias; the name is now mass-shifts
        "shift-spectrum": shift_spectrum_main,
        "featuremap": featuremap_main,
        "proteoform": proteoform_main,
    }[args.mode](rest)


if __name__ == "__main__":
    raise SystemExit(main())
