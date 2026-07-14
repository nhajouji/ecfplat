"""One-time converter: ecqf_ord_pcbij_ext.json -> ecqf_ord_ext.sqlite.

The 87MB JSON ext store (117k ordinary isogeny classes, 4 <= p < 8192)
cannot be loaded eagerly on the 512MB web instance. This packs it into a
SQLite file keyed (a, p) with per-row zlib-compressed payloads, so the
site can look classes up lazily (see pycode/ext_store.py). Rerun only
when the JSON store changes.
"""

import json
import re
import sqlite3
import zlib
from pathlib import Path

ROOT = Path(__file__).parent.parent
SRC = ROOT / "pycode" / "data" / "ecqf_ord_pcbij_ext.json"
DST = ROOT / "pycode" / "data" / "ecqf_ord_ext.sqlite"


def main():
    print(f"loading {SRC.name} …")
    ext = json.load(open(SRC))
    if DST.exists():
        DST.unlink()
    conn = sqlite3.connect(DST)
    conn.execute("CREATE TABLE bij (a INTEGER, p INTEGER, payload BLOB, "
                 "PRIMARY KEY (a, p)) WITHOUT ROWID")
    rows = []
    for key, val in ext.items():
        a, p = (int(v) for v in re.findall(r"-?\d+", key))
        blob = zlib.compress(json.dumps(val, separators=(",", ":")).encode())
        rows.append((a, p, blob))
    conn.executemany("INSERT INTO bij VALUES (?, ?, ?)", rows)
    conn.commit()
    n, pmax = conn.execute("SELECT COUNT(*), MAX(p) FROM bij").fetchone()
    conn.close()
    print(f"{DST.name}: {n} pairs, p_max = {pmax}, "
          f"{DST.stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
