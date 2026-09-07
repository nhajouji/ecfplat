"""Lazy access to the extended ordinary bijection store (p < 8192).

The full store is 117k isogeny classes. On the web instance (512 MB) it is
consulted one (a, p) pair at a time from ecqf_ord_ext.sqlite (built by
scripts/build_ext_sqlite.py); reads are from an immutable file, so a single
cross-thread connection is safe under Streamlit's threading.

Resolution order:
  1. the sqlite at $ECFPLAT_EXT_DB, if that env var is set;
  2. the sqlite at its default location, pycode/data/ecqf_ord_ext.sqlite;
  3. fallback: the JSON source ecqf_ord_pcbij_ext.json, loaded whole into
     memory once (fine on a laptop, not on the web instance -- keep the
     sqlite in the deployed repo).
Both files hold the same data; the JSON is the factory-side source of truth
and the sqlite is its deploy-friendly repack.
"""

import json
import os
import re
import sqlite3
import zlib
from pathlib import Path

_DATA_DIR = Path(__file__).parent / "data"
_DB = Path(os.environ.get("ECFPLAT_EXT_DB", _DATA_DIR / "ecqf_ord_ext.sqlite"))
_JSON = _DATA_DIR / "ecqf_ord_pcbij_ext.json"
_conn = None
_json_store = None    # {(a, p): {j: (A, B, C)}} once the fallback has loaded


def _c():
    global _conn
    if _conn is None:
        _conn = sqlite3.connect(f"file:{_DB}?mode=ro&immutable=1",
                                uri=True, check_same_thread=False)
    return _conn


def _json_fallback():
    """Load the whole JSON store once and serve lookups from memory."""
    global _json_store
    if _json_store is None:
        with open(_JSON) as f:
            raw = json.load(f)
        _json_store = {}
        for key, val in raw.items():
            a, p = (int(v) for v in re.findall(r"-?\d+", key))
            _json_store[(a, p)] = {int(j): tuple(qf) for j, qf in val.items()}
    return _json_store


def available() -> bool:
    return _DB.exists() or _JSON.exists()


def has_pair(a: int, p: int) -> bool:
    if _DB.exists():
        row = _c().execute("SELECT 1 FROM bij WHERE a = ? AND p = ?",
                           (a, p)).fetchone()
        return row is not None
    if _JSON.exists():
        return (a, p) in _json_fallback()
    return False


def get_pair(a: int, p: int):
    """{j: (A, B, C)} for the class of trace a over F_p, or None."""
    if _DB.exists():
        row = _c().execute("SELECT payload FROM bij WHERE a = ? AND p = ?",
                           (a, p)).fetchone()
        if row is None:
            return None
        val = json.loads(zlib.decompress(row[0]))
        return {int(j): tuple(qf) for j, qf in val.items()}
    if _JSON.exists():
        return _json_fallback().get((a, p))
    return None


def stats():
    """(pair count, max p) — for the coverage note."""
    if _DB.exists():
        return _c().execute("SELECT COUNT(*), MAX(p) FROM bij").fetchone()
    if _JSON.exists():
        store = _json_fallback()
        return (len(store), max(p for _, p in store))
    return (0, 0)
