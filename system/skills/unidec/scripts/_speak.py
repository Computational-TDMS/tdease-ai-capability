"""Speak: a script reports its run instead of going quiet.

Three rules, mirroring the kernel's one-error-shape doctrine:
- fail loud: any exception becomes one JSON envelope on stdout, exit 1
- every run ends in a `#` comment footer on stdout (rows, warnings, timing)
  — TSV consumers skip `#` lines, humans and agents cannot miss it
- suspicious results warn in the footer; they never just print empty output
"""

import json
import sys
import time

_START = time.perf_counter()
_NOTES = []
# Bound at import time: engine calls run under redirect_stdout(sys.stderr),
# and the report must still reach the real stdout.
_STDOUT = sys.stdout


def _emit(text):
    print(text, file=_STDOUT)


def note(message):
    _NOTES.append(str(message))


def warn(message):
    _NOTES.append(f"WARN: {message}")


def fail(code, message, **detail):
    envelope = {"ok": False, "error": {"code": code, "message": message}}
    if detail:
        envelope["error"]["detail"] = detail
    _emit(json.dumps(envelope))
    sys.exit(1)


def finish(rows, extra=None):
    lines = [f"# rows={rows}"]
    lines += [f"# {n}" for n in _NOTES]
    for key, value in (extra or {}).items():
        lines.append(f"# {key}={value}")
    lines.append(f"# seconds={time.perf_counter() - _START:.1f}")
    _emit("\n".join(lines))
