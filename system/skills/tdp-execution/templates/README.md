# tdp-execution templates

Byte-copies of `system/templates/<id>/` for guided walkthroughs inside this
skill. `system/templates/` is the single authority — when the copies drift,
the system copy wins and this directory must be refreshed. A sync test
(`tests/test_tdp_execution_skill.py`) fails when the copies diverge, so an
authoritative template change and this copy land in the same commit.
