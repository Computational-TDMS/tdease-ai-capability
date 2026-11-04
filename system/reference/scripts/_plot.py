"""Plotting for analysis scripts, split by chart type, not by skill.

- Envelope / theoretical-vs-experimental stick plots use **matplotlib** -- a
  standard chart a general library draws well. matplotlib is imported only inside
  the plotting step (`_render_matplotlib`), so a script run with `--no-plot`, or
  one that draws only hand-written SVG, stays standard-library only and starts
  fast; the numeric analysis never imports it.
- Error heatmaps and sequence-fragment ladders stay **hand-written SVG**
  (`sequence_coverage_svg`): they are domain-specific 2D layouts that a general
  library renders badly, and SVG is text, so such a chart diffs, greps and embeds
  in a report.

matplotlib is the only renderer for `Plot` figures; if it is unavailable a
figure request fails loudly rather than emitting a degraded hand-written SVG.
`--no-plot` still lets the numeric analysis run without matplotlib.

This module is duplicated across skill packages on purpose. A skill is
independently distributable; a shared import across packages would break that.

Every chart also carries its data: `Plot.series_data()` returns exactly what was
drawn, so a figure can be re-plotted elsewhere without re-running the analysis.
"""

from __future__ import annotations

import html
from dataclasses import dataclass, field

# Colour-blind safe, readable on white. Index 0 is the primary series.
PALETTE = ("#1f6feb", "#d1242f", "#1a7f37", "#9a6700", "#8250df", "#57606a")
GRID = "#d0d7de"
AXIS = "#57606a"
TEXT = "#1f2328"
MUTED = "#8c959f"


@dataclass
class Series:
    kind: str  # stem | line | bar | point
    label: str
    points: list  # [(x, y), ...]
    color: str
    annotations: list = field(default_factory=list)  # [(x, y, text), ...]


class Plot:
    """One 2D chart with a linear x axis and a linear y axis."""

    def __init__(
        self,
        title: str = "",
        x_label: str = "",
        y_label: str = "",
        width: int = 900,
        height: int = 380,
        subtitle: str = "",
    ):
        self.title = title
        self.subtitle = subtitle
        self.x_label = x_label
        self.y_label = y_label
        self.width = width
        self.height = height
        self.series: list[Series] = []
        self.markers: list = []  # vertical reference lines: (x, label, color)
        self.pad_left, self.pad_right = 74, 22
        self.pad_top = 56 if subtitle else 44
        self.pad_bottom = 54

    # ---------------------------------------------------------------- series
    def add(self, kind: str, points, label: str = "", color: str | None = None, annotations=None) -> None:
        if not points:
            return
        self.series.append(Series(
            kind=kind,
            label=label,
            points=[(float(x), float(y)) for x, y in points],
            color=color or PALETTE[len(self.series) % len(PALETTE)],
            annotations=list(annotations or []),
        ))

    def stem(self, points, label="", color=None, annotations=None) -> None:
        self.add("stem", points, label, color, annotations)

    def line(self, points, label="", color=None) -> None:
        self.add("line", points, label, color)

    def bar(self, points, label="", color=None, annotations=None) -> None:
        self.add("bar", points, label, color, annotations)

    def point(self, points, label="", color=None, annotations=None) -> None:
        self.add("point", points, label, color, annotations)

    def marker(self, x: float, label: str = "", color: str = MUTED) -> None:
        self.markers.append((float(x), label, color))

    def series_data(self) -> list:
        """The drawn data, so the figure can be reproduced without the analysis."""
        return [
            {
                "kind": item.kind,
                "label": item.label,
                "color": item.color,
                "points": [{"x": x, "y": y} for x, y in item.points],
                "annotations": [{"x": x, "y": y, "text": text} for x, y, text in item.annotations],
            }
            for item in self.series
        ]

    # ------------------------------------------------------------------ math
    def _bounds(self):
        xs = [x for item in self.series for x, _ in item.points]
        ys = [y for item in self.series for _, y in item.points]
        xs += [x for x, _, _ in self.markers]
        if not xs or not ys:
            return 0.0, 1.0, 0.0, 1.0
        x_min, x_max = min(xs), max(xs)
        y_min, y_max = min(0.0, min(ys)), max(ys)
        if x_max - x_min < 1e-12:
            x_min, x_max = x_min - 1.0, x_max + 1.0
        if y_max - y_min < 1e-12:
            y_max = y_min + 1.0
        span = x_max - x_min
        return x_min - span * 0.02, x_max + span * 0.02, y_min, y_max * 1.08

    def _scaler(self):
        x_min, x_max, y_min, y_max = self._bounds()
        left, right = self.pad_left, self.width - self.pad_right
        top, bottom = self.pad_top, self.height - self.pad_bottom

        def sx(value: float) -> float:
            return left + (value - x_min) / (x_max - x_min) * (right - left)

        def sy(value: float) -> float:
            return bottom - (value - y_min) / (y_max - y_min) * (bottom - top)

        return sx, sy, (x_min, x_max, y_min, y_max)

    @staticmethod
    def _ticks(low: float, high: float, count: int = 6) -> list:
        span = high - low
        if span <= 0:
            return [low]
        raw = span / count
        magnitude = 10 ** _floor_log10(raw)
        for step in (1, 2, 2.5, 5, 10):
            if raw <= step * magnitude:
                nice = step * magnitude
                break
        else:
            nice = 10 * magnitude
        start = nice * int(low / nice) + (nice if low > 0 and low % nice else 0)
        values, value = [], start
        while value <= high + nice * 0.001 and len(values) < 40:
            values.append(round(value, 10))
            value += nice
        return values

    # ---------------------------------------------------------------- render
    def render(self) -> str:
        sx, sy, (x_min, x_max, y_min, y_max) = self._scaler()
        left, right = self.pad_left, self.width - self.pad_right
        top, bottom = self.pad_top, self.height - self.pad_bottom
        out = [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{self.width}" height="{self.height}" '
            f'viewBox="0 0 {self.width} {self.height}" font-family="ui-sans-serif,Segoe UI,Helvetica,Arial,sans-serif">',
            f'<rect width="{self.width}" height="{self.height}" fill="#ffffff"/>',
        ]
        if self.title:
            out.append(_text(self.pad_left, 26, self.title, size=15, weight="600"))
        if self.subtitle:
            out.append(_text(self.pad_left, 43, self.subtitle, size=11, fill=MUTED))

        for value in self._ticks(y_min, y_max):
            y = sy(value)
            if not top - 1 <= y <= bottom + 1:
                continue
            out.append(f'<line x1="{left}" y1="{y:.1f}" x2="{right}" y2="{y:.1f}" stroke="{GRID}" stroke-width="1"/>')
            out.append(_text(left - 8, y + 4, _compact(value), size=10, fill=AXIS, anchor="end"))
        for value in self._ticks(x_min, x_max):
            x = sx(value)
            if not left - 1 <= x <= right + 1:
                continue
            out.append(f'<line x1="{x:.1f}" y1="{bottom}" x2="{x:.1f}" y2="{bottom + 4}" stroke="{AXIS}" stroke-width="1"/>')
            out.append(_text(x, bottom + 18, _compact(value), size=10, fill=AXIS, anchor="middle"))

        out.append(f'<line x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}" stroke="{AXIS}" stroke-width="1.2"/>')
        out.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{bottom}" stroke="{AXIS}" stroke-width="1.2"/>')
        if self.x_label:
            out.append(_text((left + right) / 2, self.height - 12, self.x_label, size=11, fill=AXIS, anchor="middle"))
        if self.y_label:
            out.append(
                f'<text x="16" y="{(top + bottom) / 2:.0f}" font-size="11" fill="{AXIS}" text-anchor="middle" '
                f'transform="rotate(-90 16 {(top + bottom) / 2:.0f})">{html.escape(self.y_label)}</text>'
            )

        for x, label, color in self.markers:
            position = sx(x)
            out.append(
                f'<line x1="{position:.1f}" y1="{top}" x2="{position:.1f}" y2="{bottom}" stroke="{color}" '
                f'stroke-width="1" stroke-dasharray="4 3"/>'
            )
            if label:
                out.append(_text(position + 3, top + 11, label, size=10, fill=color))

        baseline = sy(max(0.0, y_min))
        for item in self.series:
            out.extend(self._draw(item, sx, sy, baseline))

        labelled = [item for item in self.series if item.label]
        for index, item in enumerate(labelled):
            y = top + 12 + index * 15
            out.append(f'<rect x="{right - 152}" y="{y - 8}" width="10" height="10" fill="{item.color}"/>')
            out.append(_text(right - 138, y + 1, item.label, size=11, fill=TEXT))
        out.append("</svg>")
        return "\n".join(out)

    def _draw(self, item: Series, sx, sy, baseline) -> list:
        out = []
        if item.kind == "stem":
            for x, y in item.points:
                out.append(
                    f'<line x1="{sx(x):.1f}" y1="{baseline:.1f}" x2="{sx(x):.1f}" y2="{sy(y):.1f}" '
                    f'stroke="{item.color}" stroke-width="1.2"/>'
                )
        elif item.kind == "bar":
            count = max(1, len(item.points))
            width = max(2.0, (self.width - self.pad_left - self.pad_right) / (count * 1.6))
            for x, y in item.points:
                out.append(
                    f'<rect x="{sx(x) - width / 2:.1f}" y="{sy(y):.1f}" width="{width:.1f}" '
                    f'height="{max(0.0, baseline - sy(y)):.1f}" fill="{item.color}" opacity="0.85"/>'
                )
        elif item.kind == "line":
            path = " ".join(f"{'M' if index == 0 else 'L'}{sx(x):.1f},{sy(y):.1f}"
                            for index, (x, y) in enumerate(item.points))
            out.append(f'<path d="{path}" fill="none" stroke="{item.color}" stroke-width="1.6"/>')
        else:
            for x, y in item.points:
                out.append(f'<circle cx="{sx(x):.1f}" cy="{sy(y):.1f}" r="3.2" fill="{item.color}"/>')
        for x, y, text in item.annotations:
            out.append(_text(sx(x), sy(y) - 6, text, size=9, fill=item.color, anchor="middle"))
        return out


def sequence_coverage_svg(
    sequence: str,
    n_terminal_cuts,
    c_terminal_cuts,
    *,
    highlight=(),
    title: str = "",
    subtitle: str = "",
    per_row: int = 40,
) -> str:
    """Fragment coverage drawn on the sequence itself.

    A cut after residue k is drawn as a tick between residue k and k+1: down-left
    for an N-terminal ion, up-right for a C-terminal one. Reading where the ticks
    stop is how you see whether a site is actually bracketed, which no generic
    chart type shows.
    """
    n_cuts, c_cuts, marked = set(n_terminal_cuts), set(c_terminal_cuts), set(highlight)
    rows = [sequence[i:i + per_row] for i in range(0, len(sequence), per_row)] or [""]
    cell, row_height = 19, 62
    width = 74 + per_row * cell
    height = 46 + len(rows) * row_height + (14 if subtitle else 0)
    top_offset = 40 + (14 if subtitle else 0)

    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" font-family="ui-monospace,SFMono-Regular,Consolas,monospace">',
        f'<rect width="{width}" height="{height}" fill="#ffffff"/>',
    ]
    if title:
        out.append(_text(20, 24, title, size=15, weight="600", family="ui-sans-serif,Segoe UI,Helvetica,Arial"))
    if subtitle:
        out.append(_text(20, 40, subtitle, size=11, fill=MUTED, family="ui-sans-serif,Segoe UI,Helvetica,Arial"))

    for row_index, row in enumerate(rows):
        base_y = top_offset + row_index * row_height + 22
        start = row_index * per_row
        out.append(_text(20, base_y + 4, f"{start + 1}", size=10, fill=MUTED))
        for column, residue in enumerate(row):
            position = start + column + 1
            x = 68 + column * cell
            if position in marked:
                out.append(f'<rect x="{x - 8}" y="{base_y - 13}" width="17" height="18" rx="3" fill="#fff3cd"/>')
            out.append(_text(x, base_y, residue, size=13, fill=TEXT, anchor="middle"))
            if position in n_cuts and column < len(row):
                out.append(
                    f'<path d="M{x + 9},{base_y - 15} v-9 h-7" fill="none" stroke="{PALETTE[0]}" stroke-width="1.6"/>'
                )
            if position in c_cuts and column < len(row):
                out.append(
                    f'<path d="M{x + 9},{base_y + 7} v9 h7" fill="none" stroke="{PALETTE[1]}" stroke-width="1.6"/>'
                )
    out.append(_text(width - 210, height - 12, "N-terminal ions", size=10, fill=PALETTE[0],
                     family="ui-sans-serif,Segoe UI,Helvetica,Arial"))
    out.append(_text(width - 100, height - 12, "C-terminal ions", size=10, fill=PALETTE[1],
                     family="ui-sans-serif,Segoe UI,Helvetica,Arial"))
    out.append("</svg>")
    return "\n".join(out)


def _text(x, y, value, *, size=12, fill=TEXT, anchor="start", weight="400", family=None) -> str:
    extra = f' font-family="{family}"' if family else ""
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" fill="{fill}" text-anchor="{anchor}" '
        f'font-weight="{weight}"{extra}>{html.escape(str(value))}</text>'
    )


def _compact(value: float) -> str:
    if value == 0:
        return "0"
    magnitude = abs(value)
    if magnitude >= 1e6 or magnitude < 1e-3:
        return f"{value:.1e}"
    if magnitude >= 1000:
        return f"{value:,.0f}"
    if magnitude >= 10:
        return f"{value:.1f}".rstrip("0").rstrip(".")
    return f"{value:.4g}"


def _floor_log10(value: float) -> int:
    exponent = 0
    if value <= 0:
        return 0
    while value < 1:
        value *= 10
        exponent -= 1
    while value >= 10:
        value /= 10
        exponent += 1
    return exponent


# --------------------------------------------------------------------- output


def add_plot_arguments(target) -> None:
    """The visualization flags every analysis script shares."""
    target.add_argument("--plot", help="write the figure here; defaults to an .svg beside --out")
    target.add_argument("--plot-pdf", help="also write the primary matplotlib figure as vector PDF")
    target.add_argument("--no-plot", action="store_true", help="skip figures entirely")
    target.add_argument(
        "--plot-data", action="store_true",
        help="also return the drawn series in the JSON, so the figure can be rebuilt elsewhere",
    )


def _figure_extension(renderable) -> str:
    """PNG for a Plot drawn with matplotlib, SVG for hand-written markup."""
    return "png" if isinstance(renderable, Plot) else "svg"


def _render_matplotlib(plot: "Plot", path) -> str:
    """Render a Plot with matplotlib. Returns the media type written.

    matplotlib is imported here, never at module load, so a script that runs
    --no-plot -- or one that only draws hand-written SVG -- stays standard-library
    only and starts fast.

    This is the renderer for envelope / theoretical-vs-experimental stick plots,
    where a real plotting library earns its dependency. Error heatmaps and
    sequence-fragment ladders stay hand-drawn (see sequence_coverage_svg); a
    general library lays those domain-specific layouts out badly.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    suffix = path.suffix.lower().lstrip(".")
    fmt = suffix if suffix in ("png", "svg", "pdf") else "png"
    fig, ax = plt.subplots(figsize=(plot.width / 100.0, max(2.4, plot.height / 100.0)), dpi=150)
    ax.set_axisbelow(True)
    ax.grid(axis="y", color=GRID, linewidth=0.6, alpha=0.65)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(AXIS)
    ax.spines["bottom"].set_color(AXIS)
    ax.tick_params(axis="both", labelsize=8, colors=AXIS)
    xs_all = [x for item in plot.series for x, _ in item.points]
    span = (max(xs_all) - min(xs_all)) if len(xs_all) > 1 else 1.0
    for item in plot.series:
        if not item.points:
            continue
        xs = [x for x, _ in item.points]
        ys = [y for _, y in item.points]
        label = item.label or None
        if item.kind == "bar":
            width = (span / max(len(xs), 1)) * 0.7 if span else 0.7
            ax.bar(xs, ys, width=width, color=item.color, alpha=0.85, label=label)
        elif item.kind == "line":
            ax.plot(xs, ys, color=item.color, linewidth=1.6, label=label)
        elif item.kind == "point":
            ax.scatter(xs, ys, color=item.color, s=16, alpha=0.88, linewidths=0, label=label)
        else:  # stem
            ax.vlines(xs, 0.0, ys, color=item.color, linewidth=1.3, label=label)
            ax.scatter(xs, ys, color=item.color, s=11, linewidths=0)
        for ax_x, ax_y, text in item.annotations:
            ax.annotate(
                str(text), (ax_x, ax_y), textcoords="offset points", xytext=(0, 4),
                ha="center", fontsize=7, color=item.color,
            )
    heading = plot.title + (f"\n{plot.subtitle}" if plot.subtitle else "")
    if heading.strip():
        ax.set_title(heading, fontsize=10, loc="left")
    if plot.x_label:
        ax.set_xlabel(plot.x_label, fontsize=9)
    if plot.y_label:
        ax.set_ylabel(plot.y_label, fontsize=9)
    if any(item.label for item in plot.series):
        ax.legend(fontsize=8, frameon=False)
    ax.margins(x=0.02)
    fig.tight_layout()
    fig.savefig(path, format=fmt, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return {"svg": "image/svg+xml", "pdf": "application/pdf"}.get(fmt, "image/png")


def _write_figure(renderable, path):
    """Write one figure. Returns (written_path, media_type, renderer).

    A hand-written SVG string is written as-is (standard library) — this is the
    path for the domain-specific layouts (`sequence_coverage_svg`, error
    heatmaps) that a general chart library renders badly. A `Plot` object is
    always drawn with matplotlib; if matplotlib is unavailable that raises, so
    an analysis script fails loudly instead of emitting a degraded hand-written
    SVG. --plot-data still carries the exact series either way.
    """
    if isinstance(renderable, Plot):
        media_type = _render_matplotlib(renderable, path)
        return path, media_type, "matplotlib"
    path.write_text(renderable, encoding="utf-8")
    return path, "image/svg+xml", "hand-svg"


def emit_figures(result: dict, args, figures) -> None:
    """Write figures and record them on the result.

    `figures` is a sequence of (suffix, kind, title, renderable, series). The first
    entry is the primary figure. `renderable` is either a hand-written SVG string
    (error heatmaps, sequence-fragment ladders) or a Plot object (envelope stick
    plots, rendered with matplotlib as PNG). Nothing is written unless the caller
    asked for an output path, either with --plot or with --out; visualization is on
    by default but never invents a location. --plot-data records the exact drawn
    series regardless of --no-plot.
    """
    from pathlib import Path as _Path

    figures = [item for item in figures if item and item[3]]
    if not figures:
        return
    if getattr(args, "plot_data", False):
        result["figure_data"] = [
            {"kind": kind, "title": title, "series": series} for _, kind, title, _, series in figures
        ]
    if getattr(args, "no_plot", False):
        return

    explicit = getattr(args, "plot", None)
    out = getattr(args, "out", None) or getattr(args, "out_dir", None)
    primary_ext = _figure_extension(figures[0][3])
    if explicit:
        primary = _Path(explicit)
    elif out:
        base = _Path(out)
        if base.is_dir() or getattr(args, "out_dir", None):
            primary = base / f"figure.{primary_ext}"
        else:
            primary = base.with_suffix(f".{primary_ext}")
    else:
        return

    written = []
    for index, (suffix, kind, title, renderable, _series) in enumerate(figures):
        ext = _figure_extension(renderable)
        target = primary if index == 0 else primary.with_suffix(f".{suffix}.{ext}")
        target.parent.mkdir(parents=True, exist_ok=True)
        written_path, media_type, renderer = _write_figure(renderable, target)
        written.append({
            "path": str(written_path), "kind": kind, "title": title,
            "media_type": media_type, "renderer": renderer,
        })
    pdf = getattr(args, "plot_pdf", None)
    if pdf:
        renderable = figures[0][3]
        if not isinstance(renderable, Plot):
            raise ValueError("--plot-pdf requires a matplotlib Plot figure")
        target = _Path(pdf)
        target.parent.mkdir(parents=True, exist_ok=True)
        written_path, media_type, renderer = _write_figure(renderable, target)
        written.append({
            "path": str(written_path), "kind": f"{figures[0][1]}-pdf",
            "title": figures[0][2], "media_type": media_type, "renderer": renderer,
        })
    result["figures"] = written
