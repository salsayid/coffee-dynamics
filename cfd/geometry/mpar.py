"""Write ANSA .ansa_mpar files from the saved defaults with a few keys overridden."""
import os
import re

HERE = os.path.dirname(__file__)
DEFAULTS = {
    "surface": os.path.join(HERE, "ansa_defaults", "bm_surface.ansa_mpar"),
    "layers": os.path.join(HERE, "ansa_defaults", "bm_layers.ansa_mpar"),
    "volume": os.path.join(HERE, "ansa_defaults", "bm_volume.ansa_mpar"),
}


def write(kind, out, **over):
    text = open(DEFAULTS[kind]).read().splitlines()
    seen = set()
    for i, line in enumerate(text):
        m = re.match(r"^(\s*)([^#=][^=]*?)(\s*)=(.*)$", line)
        if not m:
            continue
        key = m.group(2).strip()
        if key in over:
            text[i] = f"{m.group(1)}{m.group(2)}{m.group(3)}= {over[key]}"
            seen.add(key)
    missing = set(over) - seen
    if missing:
        raise KeyError(f"keys not in {kind} defaults: {sorted(missing)}")
    open(out, "w").write("\n".join(text) + "\n")
    return out
