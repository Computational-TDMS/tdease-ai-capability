# Volcano tutorial — build one from scratch, and get the encodings right

Teaching companion to `volcano.md` (the spec). Two reference implementations ship with
this skill — read them alongside this walkthrough:
- `examples/volcano-site-stats/` — synthetic, **shape channel** demo (class → glyph);
  its click-selection card and presence panel follow the interaction spec

## Step 1 — audit the data before any pixel

Columns you need: identifier, log2FC, p (or padj), and *optionally* a tertiary measure
(mean abundance), a class, a presence tag. Count: quantified vs single-presence rows —
the excluded count goes into the subtitle, never silently dropped
(e.g. "N = 512 quantified (1387 single-presence excluded) · 96 significant" — the real numbers from the run, never rounded or invented).

## Step 2 — scales: symmetric x, headroomed y

```js
const xMax = Math.max(...data.map(d => Math.abs(d.log2fc))) * 1.05;   // 0 centred
const yMax = Math.ceil(Math.max(...data.map(d => -Math.log10(d.padj)))) + 0.5; // label headroom
const x = d3.scaleLinear().domain([-xMax, xMax]).range([M.l, W - M.r]);
const y = d3.scaleLinear().domain([0, yMax]).range([H - M.b, M.t]);
```

## Step 3 — the audit trail (thresholds + tints + subtitle)

Dashed cutoff lines, corner tints at alpha 0.05, and a subtitle stating N, both cutoffs
and the adjustment method. A volcano without this is decoration.

## Step 4 — the encoding budget (the part everyone gets wrong)

A volcano has **two position channels already spent** (x=FC, y=significance). What
remains, in Cleveland-McGill accuracy order:

| Channel | Accuracy | Use for | Rules |
| --- | --- | --- | --- |
| colour | medium (nominal) | the 4-state verdict (NS/P/F/P+F) | mandatory legend; one meaning only |
| shape | medium-low | a class variable (≤4 glyphs) | grayscale-surviving glyphs; legend mandatory |
| **size (area)** | **low** | a tertiary measure (mean abundance) | `scaleSqrt`, range 1.6–5.5pt, **size legend (Q1/med/Q3) mandatory** |
| alpha | lowest | **nothing — keep constant 0.7–0.8** | never data-mapped on top of size |

**Budget rule: at most ONE soft channel beyond colour.** Size + alpha both data-mapped
compete and neither decodes. In a size-channel volcano, size encodes mean abundance and
alpha stays 0.75 — and the honest caveat: that dataset's abundances span only 22–25
(log2), so the size channel is intentionally subtle. When the tertiary variable's spread
is narrow, say so (or drop the channel) rather than exaggerating the range — a size
legend with three nearly-equal bubbles is the visual confession.

**Size mapping mechanics** (the usual bugs): use `d3.scaleSqrt` — radius linear scaling
quadruples area for 2× data (lies); clamp range to 1.6–5.5pt so labels stay placeable;
`NaN`-safe radius fallback (~2.2pt) for rows missing the measure.

## Step 5 — emphasis without lying

- Mass points: **no stroke** (strokes fake density).
- Significant subset: thin own-colour darker stroke (0.9px) — emphasis by membership,
  not by a new hue (the size-channel reference does exactly this).
- Labels: top-N significant by padj (≤10), greedy 8-offset repel, right-edge flip
  (a clipped label is a bug — we shipped one and fixed it; see QA list).

## Step 6 — interaction that earns its keep

- Tooltip: identifier, presence/class, all measures, formatted (exponential for tiny p).
- **Click → pinned selection**: ring on the point + info card inside the plot; click
  empty space releases (`select()`). Tooltips evaporate; a
  selection card is what reviewers cite in emails.
- Hover hit areas: transparent column rects (±5px), not bare circle hit boxes.

## Step 7 — the excluded data gets a face

Single-presence rows can't be placed (no FC) but they are results. Give them a count
panel inside the plot ("single-presence (not placed) · WT-only ×1067 · M-only ×325")
— the subtitle claims it, the panel shows it.

## Step 8 — title discipline (learned the hard way)

The title states ONLY what the figure supports. One real-data title claimed
"family X separates the two conditions" — knowledge imported from elsewhere in the
project, not from this plot. Corrected to what the canvas itself showed: the broad
direction of change plus the single strongest outlier. An insight-led title backed by
nothing on the canvas is a claim without evidence.

## QA checklist (runs against every volcano)

1. Subtitle: N, excluded, both cutoffs, adjustment method.
2. Four-state colours + legends for every active channel (colour/shape/size).
3. No label collisions; no label clipped at edges; min font ≥5pt at print size.
4. Significant subset emphasized; mass cloud strokeless.
5. Size legend matches `scaleSqrt` truth; alpha constant.
6. Click-select works; tooltip fields complete.
7. Greyscale decode still works (shape channel survives).
8. Provenance footer auto-rendered from meta; `--meta colour_scale` matches reality.
