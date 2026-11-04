"""Theory-pair mode: TV-distance / apex comparison of two theoretical
isotope distributions per panel (kernel-free)."""

import json
from pathlib import Path

from _shared.manifest import write_manifest
from envelope_common import (
    FigureError,
    _compare_distributions,
    _load_theory_distribution,
    _render,
)


def run(args) -> int:
    if args.pair_label and len(args.pair_label) != len(args.theory_pair):
        raise FigureError("PAIR_LABEL_MISMATCH: provide one --pair-label per --theory-pair")
    labels = args.pair_label or [f"pair {index + 1}" for index in range(len(args.theory_pair))]
    comparisons = []
    for label, (left_path, right_path) in zip(labels, args.theory_pair):
        comparisons.append(_compare_distributions(
            _load_theory_distribution(left_path),
            _load_theory_distribution(right_path),
            label,
        ))
    plot_path = args.plot or (Path(args.theory_pair[0][0]).stem + "-theory-comparison.png")
    panels = [{
        "title": item["label"],
        "observed": {"x": [], "y": []},
        "series": [
            {"label": args.left_label, "x": list(range(len(item["left"]))), "y": item["left"]},
            {"label": args.right_label, "x": list(range(len(item["right"]))), "y": item["right"]},
        ],
        "apex": {"observed": None,
                 "series": [item["left_apex_index"], item["right_apex_index"]]},
    } for item in comparisons]
    payload = {"data": {"panels": panels,
                        "xlabel": "Isotope index (M+n)",
                        "ylabel": "Normalized probability"},
               "meta": {"title": args.label or "Theoretical isotope distributions",
                        "n": f"{len(panels)} pair panels"}}
    written = _render(plot_path, payload)
    output_comparisons = [{
        "label": item["label"],
        "tv_distance": round(item["tv_distance"], 6),
        "left_apex_index": item["left_apex_index"],
        "right_apex_index": item["right_apex_index"],
    } for item in comparisons]
    write_manifest(
        plot_path, figure="theory-comparison", script=Path(__file__).name,
        params={
            "left_label": args.left_label, "right_label": args.right_label,
            "pair_labels": labels,
        },
        inputs=[path for pair in args.theory_pair for path in pair],
        caption={"comparisons": output_comparisons},
        plot_data=output_comparisons,
    )
    print(json.dumps({
        "figure": "theory-comparison",
        "plot": str(Path(plot_path).with_suffix(".html").resolve()),
        "renderer": "chart-forge d3",
        "comparisons": output_comparisons,
        "written": written,
    }, ensure_ascii=False))
    return 0
