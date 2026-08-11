"""Prototype: build the Explorer's class view as a static site.

Proves the load-bearing claim behind a Netlify port — that the (a, p) view
needs no Python at visit time. Everything the page shows is computed here, at
build time, and shipped as JSON; the browser only draws.

The trick is that explorer_viz's applets are already data-driven: each emits
"const DATA = <json>;" inside an IIFE. We generate one applet, split it at
that line into (template, data), keep the template once as shared JS, and
write only the per-class data into a per-prime JSON shard. The drawing code
is therefore byte-identical to what the live site runs — this prototype
inherits the applets rather than reimplementing them.

    python3 netlify-proto/build.py 101 1009

Output lands in netlify-proto/public/, which is what you'd hand to Netlify.
"""
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
sys.path.insert(0, str(REPO / "pycode"))

from ecqf import ECQFIsogenyClass, class_graph_descriptor      # noqa: E402
from ecqf_tools import ec_eq_str_base, abc_to_tau, abc_to_tau_str  # noqa: E402
from nt import primeQ, primefact                               # noqa: E402
from palette import row_colors                                 # noqa: E402
import explorer_viz                                            # noqa: E402

OUT = ROOT / "public"
DATA_DIR = OUT / "data"
GUARD_NODES = 600            # mirrors GUARDS["volcano_nodes"] in the app


# ── splitting a generated applet into (template, data) ───────────────────────

def split_applet(html: str, init_name: str):
    """Turn a generated applet into a function we can call with fetched data.

    Returns (markup+script template, data dict). The applet's IIFE becomes
    window.<init_name>(DATA) so the page can invoke it once the JSON arrives.
    """
    lines = html.split("\n")
    hits = [i for i, ln in enumerate(lines) if ln.startswith("const DATA = ")]
    if len(hits) != 1:
        raise RuntimeError(f"expected exactly one 'const DATA =' line, got {len(hits)}")
    i = hits[0]
    data = json.loads(lines[i][len("const DATA = "):].rstrip(";"))
    del lines[i]
    tmpl = "\n".join(lines)

    # the IIFE wrapper -> a named init function
    open_iife = "<script>\n(() => {"
    if open_iife not in tmpl:
        raise RuntimeError("applet does not open with the expected IIFE")
    tmpl = tmpl.replace(open_iife, f"<script>\nwindow.{init_name} = (DATA) => {{", 1)

    close_iife = "})();\n</script>"
    if tmpl.count(close_iife) != 1:
        raise RuntimeError("applet does not close with the expected IIFE")
    tmpl = tmpl.replace(close_iife, "};\n</script>", 1)
    return tmpl, data


# ── per-class payload (mirrors class_view in pages/1_Explorer.py) ────────────

def ls_for(cls):
    ls = {2, 3}
    if cls.cond > 1:
        ls.update(int(q) for q in primefact(cls.cond))
    return sorted(ls)


def centered(c: int, p: int) -> int:
    c %= p
    return c - p if 2 * c > p else c


def fd_points(cls, p: int):
    """The class's CM points, ordered like qfs_ordered so node indices agree."""
    n = len(cls.qfs_all)
    df = cls.ecqf_df() if cls.js_to_qf is not None else None
    colors = row_colors(n)
    qf_to_color, qf_to_js, qf_to_fg = {}, {}, {}
    if df is not None:
        for i, (_, row) in enumerate(df.iterrows()):
            qf = tuple(row["qf_coefs"])
            qf_to_color.setdefault(qf, colors[i])
            qf_to_js.setdefault(qf, []).append(row["j_inv"])
            qf_to_fg.setdefault(qf, tuple(int(v) for v in row["EC_coefs"]))
    pts = []
    for qf in cls.qfs_ordered:
        x, y = abc_to_tau(qf)
        js = qf_to_js.get(qf)
        label = f"j={js[0]}" + ("…" if js and len(js) > 1 else "") if js else str(qf)
        sub = (f"⟨{', '.join(str(v) for v in qf)}⟩ · End disc {cls.endo_disc_dict[qf]}"
               f" · τ = {abc_to_tau_str(qf)}")
        if js:
            sub += " · j = " + ", ".join(str(j) for j in js)
        if qf in qf_to_fg:
            sub += " · " + ec_eq_str_base(tuple(centered(v, p) for v in qf_to_fg[qf]))
        pts.append({"x": float(round(x, 5)), "y": float(round(y, 5)),
                    "color": qf_to_color.get(qf, "#4da3d8"),
                    "label": label, "sub": sub})
    return pts


def node_js(cls):
    """j-invariants per lattice class, indexed like qfs_ordered.

    A lattice class can carry several curves (quadratic twists), so this is a
    list per node. The first is the canonical address for that node's link.
    """
    out = [[] for _ in cls.qfs_ordered]
    if cls.js_to_qf is None:
        return out
    idx_of = {qf: i for i, qf in enumerate(cls.qfs_ordered)}
    for jsig, qf in cls.js_to_qf.items():
        j = jsig[0] if isinstance(jsig, tuple) else jsig
        out[idx_of[qf]].append(int(j))
    return [sorted(v) for v in out]


def table_rows(cls):
    """One row per curve, tagged with the index of its lattice class.

    Several curves can share a lattice class (quadratic twists), so `idx` is
    not unique across rows — it is the join key onto qfs_ordered, which is what
    the volcano and the CM points index. That is what lets the table take part
    in the shared selection.
    """
    if cls.js_to_qf is None:
        return None
    idx_of = {qf: i for i, qf in enumerate(cls.qfs_ordered)}
    rows = []
    for _, r in cls.ecqf_df().iterrows():
        qf = tuple(r["qf_coefs"])
        rows.append({
            "idx": idx_of[qf],
            "sig": str(r["ec_invs"]),
            "j": int(r["j_inv"]),
            "fg": [int(v) for v in r["EC_coefs"]],
            "qf": [int(v) for v in qf],
            "ed": int(r["endo_disc"]),
            "ec": int(r["endo_cond"]),
            "cc": int(r["endo_cocond"]),
            "frob": str(r["frobmat"]),
            "tau": str(r["tau_s"]),
        })
    return {"rows": rows,
            "discs": sorted({r["ed"] for r in rows})}


FD_POINTS_GUARD = 2000       # mirrors GUARDS["fd_points"] in the app


def curve_payloads(cls, a: int, p: int, templates: dict):
    """One entry per j-invariant — the address the static site uses.

    The torus applet takes (qf, a, p, frobmat, n_points) and derives the points
    in the browser, so each entry costs ~0.1 KB regardless of #E(F_p). Curve
    data therefore rides along in the per-prime shard instead of needing files
    of its own.
    """
    if cls.js_to_qf is None:
        return {}
    idx_of = {qf: i for i, qf in enumerate(cls.qfs_ordered)}
    js_of = node_js(cls)
    N = p + 1 - a
    n_pts = N if N <= FD_POINTS_GUARD else 0
    out = {}
    for jsig, qf in cls.js_to_qf.items():
        j = int(jsig[0] if isinstance(jsig, tuple) else jsig)
        frm = cls.qf_to_frob_mats[qf]
        tmpl, data = split_applet(
            explorer_viz.curve_torus_html(qf, a, p, frm.vec, n_pts), "initTorus")
        templates.setdefault("torus", tmpl)
        model = cls.js_to_models.get(jsig)
        idx = idx_of[qf]
        out[str(j)] = {
            "j": j, "idx": idx, "a": a, "p": p,
            "qf": [int(v) for v in qf],
            "model": [centered(int(v), p) for v in model] if model else None,
            "endoDisc": int(cls.endo_disc_dict[qf]),
            "endoCond": int(cls.endo_cond_dict[qf]),
            "tau": abc_to_tau_str(qf),
            "frob": str(frm.vec),
            # other curves on the same lattice class — quadratic twists
            "twins": [v for v in js_of[idx] if v != j],
            "torus": data,
        }
    return out


def class_payload(a: int, p: int, templates: dict):
    cls = ECQFIsogenyClass(a, p)
    n = len(cls.qfs_all)
    if n > GUARD_NODES:
        return None
    d = a * a - 4 * p

    # Address curves by j-invariant, not by position in qfs_ordered: the
    # position is an implementation detail and would silently repoint every
    # saved link if the ordering ever changed. Nodes with no curve tables get
    # no link rather than a broken one.
    js = node_js(cls)
    hrefs = [f"curve.html?p={p}&j={v[0]}" if v else None for v in js]

    graph_html = explorer_viz.isogeny_graph_html(
        {l: class_graph_descriptor(cls, l) for l in ls_for(cls)},
        node_hrefs=hrefs, height_px=740)
    fd_html = explorer_viz.fd_points_html(fd_points(cls, p),
                                          node_hrefs=hrefs, height_px=560)

    g_tmpl, g_data = split_applet(graph_html, "initGraph")
    f_tmpl, f_data = split_applet(fd_html, "initFD")
    templates.setdefault("graph", g_tmpl)
    templates.setdefault("fd", f_tmpl)

    chi = (f"x^2 + {p}" if a == 0
           else f"x^2 + {-a}x + {p}" if a < 0 else f"x^2 - {a}x + {p}")
    return {
        "a": a, "p": p, "d": d, "chi": chi, "N": p + 1 - a,
        "fieldDisc": cls.field_disc, "cond": cls.cond, "n": n,
        "supersingular": a == 0,
        "hasCurves": cls.js_to_qf is not None,
        "graph": g_data, "fd": f_data, "table": table_rows(cls),
        "curves": curve_payloads(cls, a, p, templates),
    }


# ── build ────────────────────────────────────────────────────────────────────

def build(primes):
    if OUT.exists():
        shutil.rmtree(OUT)
    DATA_DIR.mkdir(parents=True)
    templates, index = {}, {}

    for p in primes:
        if not primeQ(p):
            raise SystemExit(f"{p} is not prime")
        shard, a = {}, -int((4 * p) ** 0.5) - 1
        while a * a >= 4 * p:
            a += 1
        while a * a < 4 * p:
            if a == 0 or a % p != 0:
                payload = class_payload(a, p, templates)
                if payload is not None:
                    shard[str(a)] = payload
            a += 1
        path = DATA_DIR / f"{p}.json"
        path.write_text(json.dumps(shard, separators=(",", ":")))
        index[str(p)] = sorted(int(k) for k in shard)
        print(f"  p={p:<6} {len(shard):>4} classes  {path.stat().st_size/1024:8.1f} KB")

    (DATA_DIR / "index.json").write_text(json.dumps(index, separators=(",", ":")))
    (OUT / "class.html").write_text(
        (ROOT / "class.tmpl.html").read_text()
        .replace("__GRAPH_APPLET__", templates["graph"])
        .replace("__FD_APPLET__", templates["fd"]))
    (OUT / "curve.html").write_text(
        (ROOT / "curve.tmpl.html").read_text()
        .replace("__TORUS_APPLET__", templates.get("torus", "")))
    (OUT / "index.html").write_text((ROOT / "index.tmpl.html").read_text())
    # pretty URLs are a host concern; these are what the design doc specifies
    (OUT / "_redirects").write_text(
        "/explorer/c/:p/:a       /class.html?p=:p&a=:a       200\n"
        "/explorer/c/:p/:a/j:j   /curve.html?p=:p&j=:j       200\n"
        "/Explorer               /class.html                 200\n")
    print(f"\nbuilt -> {OUT}")


if __name__ == "__main__":
    args = [int(x) for x in sys.argv[1:]] or [101, 1009]
    build(args)
