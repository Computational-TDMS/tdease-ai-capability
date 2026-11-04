"""TopPIC unexpected-modification column grammar shared by result readers."""
import math
import re


MOD_RE = re.compile(
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)\s*:\s*\[(\d+)(?:-(\d+))?\]"
)


def parse_mods(value, *, strict=False):
    text = str(value or "").strip()
    if text in {"", "-", "0"}:
        return []
    matches = list(MOD_RE.finditer(text))
    if strict and (not matches or MOD_RE.sub("", text).strip(" ;,\t\r\n")):
        raise ValueError(f"unrecognized unexpected-modification syntax: {text!r}")
    result = []
    for match in matches:
        mass, start, end = float(match[1]), int(match[2]), int(match[3] or match[2])
        if not math.isfinite(mass) or start < 1 or end < start:
            raise ValueError(f"invalid unexpected modification: {match[0]!r}")
        result.append((mass, start, end, start == 1))
    return result
