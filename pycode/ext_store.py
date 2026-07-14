"""Lazy access to the extended ordinary bijection store (p < 8192).

The full store is 117k isogeny classes — far too big to hold in memory on
the web instance — so it lives in ecqf_ord_ext.sqlite (built by
scripts/build_ext_sqlite.py) and is consulted one (a, p) pair at a time.
Reads are from an immutable file, so a single cross-thread connection is
safe under Streamlit's threading."""

import json
import sqlite3
import zlib
from pathlib import Path

_DB = Path(__file__).parent / "data" / "ecqf_ord_ext.sqlite"
_conn = None


def _c():
    global _conn
    if _conn is None:
        _conn = sqlite3.connect(f"file:{_DB}?mode=ro&immutable=1",
                                uri=True, check_same_thread=False)
    return _conn


def available() -> bool:
    return _DB.exists()


def has_pair(a: int, p: int) -> bool:
    if not available():
        return False
    row = _c().execute("SELECT 1 FROM bij WHERE a = ? AND p = ?",
                       (a, p)).fetchone()
    return row is not None


def get_pair(a: int, p: int):
    """{j: (A, B, C)} for the class of trace a over F_p, or None."""
    if not available():
        return None
    row = _c().execute("SELECT payload FROM bij WHERE a = ? AND p = ?",
                       (a, p)).fetchone()
    if row is None:
        return None
    val = json.loads(zlib.decompress(row[0]))
    return {int(j): tuple(qf) for j, qf in val.items()}


def stats():
    """(pair count, max p) — for the coverage note."""
    if not available():
        return (0, 0)
    return _c().execute("SELECT COUNT(*), MAX(p) FROM bij").fetchone()
