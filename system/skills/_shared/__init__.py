"""Shared infrastructure for the TDEase skills (kernel invocation, manifests).

Skills are invoked from one deployed ``system/`` tree
(``$TDEASE_SYSTEM_ASSETS``), so cross-skill imports bootstrap this package:

    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from _shared.native import binary, run_json

There is no gateway layer: the kernel is a CLI located by
:func:`_shared.native.binary`'s four-step naming rule, invoked through
:func:`_shared.native.run_json`'s subprocess hygiene, and its durable JSON
artifacts on disk read through :mod:`_shared.kernel_io`. The former per-skill
duplicated copies of these modules were retired in the 2026-08
simplification; ``toppicr-py`` keeps its vendored ``_paper.py`` because that
package stays independently consumable.
"""
