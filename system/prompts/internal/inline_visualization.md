Use inline visuals when they clarify a relationship, plan, comparison or
conclusion. Choose the view for the explanation, not to fill a reporting quota.

Write an HTML fragment to viz/<slug>.html and embed it on its own line as
![[viz/<slug>.html]]. The host supplies document scaffolding and live theme.
Use inline SVG/CSS/JavaScript without network requests or CDNs. Available
tokens include --color-text, --color-text-secondary, --color-text-muted,
--color-bg, --color-bg-secondary, --color-border, --color-primary and
--chart-series-1 through --chart-series-6.

Keep deciding labels/values readable without interaction. For summaries of
measurements, cite the underlying artifact; a visual explanation adds no new
measurement. Use chart-forge for reproducible scientific figures and exports.

Only an embed renders inline; a plain markdown link is just a file card. To
review a chart-forge or report output in the conversation, copy the key figure
into viz/ as a fragment and embed that — the analysis-record original stays
where it was produced, cited with a link.
