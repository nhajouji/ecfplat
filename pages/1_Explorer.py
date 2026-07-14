"""The Explorer — one drill-down page for the whole (a, p) / discriminant world.

Navigation is by query params, so every view is a shareable URL:
    (nothing)          the two entry doors
    ?d=-231            discriminant view (characteristic 0)
    ?a=6&p=101         class view (Frobenius view over F_p)
    ?a=6&p=101&node=4  curve view (one lattice class / one curve)
Canvas widgets link between views with target=_parent anchors.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "pycode"))

import streamlit as st
import streamlit.components.v1 as components

from nt import primeQ, primefact, discfac, primesBetween
from ecqf import (QFIsogenyClass, ECQFIsogenyClass, class_graph_descriptor,
                  disc_to_aps, GUARDS, P_MAX)
from ecqf_tools import (ec_eq_str_base, frob_to_mw_gens, abc_to_tau,
                        abc_to_tau_str, ec_look_up, ap_in_pc_data)
from palette import row_colors
import explorer_viz

STORE_P_MAX = 8191          # curve tables cover 99.86% of classes to here
                            # (the rest render structure-only, gracefully)

# ── helpers ───────────────────────────────────────────────────────────────────

@st.cache_resource(max_entries=24, show_spinner=False)
def _load_ec_class(a: int, p: int) -> ECQFIsogenyClass:
    return ECQFIsogenyClass(a, p)


@st.cache_resource(max_entries=24, show_spinner=False)
def _load_qf_class(d: int) -> QFIsogenyClass:
    return QFIsogenyClass(d)


def _qp_int(key):
    try:
        return int(st.query_params[key])
    except (KeyError, ValueError):
        return None


def _ls_for(cls) -> list:
    """Isogeny degrees offered by the graph widget: always 2 and 3, plus every
    prime dividing the conductor (those are the volcano directions)."""
    ls = {2, 3}
    if cls.cond > 1:
        ls.update(int(q) for q in primefact(cls.cond))
    return sorted(ls)


def _descs_for(cls, link_qs: str = None):
    return {l: class_graph_descriptor(cls, l) for l in _ls_for(cls)}


def _centered(c: int, p: int) -> int:
    c %= p
    return c - p if 2 * c > p else c


def _valid_disc(d) -> bool:
    return d is not None and d < 0 and d % 4 in (0, 1)


def _long_to_short(a1, a2, a3, a4, a6, p):
    """(a1..a6) mod p -> the short model (A, B): complete square and cube.
    Valid for p >= 5 (needs 48, 864 invertible)."""
    b2 = (a1 * a1 + 4 * a2) % p
    b4 = (2 * a4 + a1 * a3) % p
    b6 = (a3 * a3 + 4 * a6) % p
    c4 = (b2 * b2 - 24 * b4) % p
    c6 = (-b2 ** 3 + 36 * b2 * b4 - 216 * b6) % p
    A = (-c4 * pow(48, p - 2, p)) % p
    B = (-c6 * pow(864, p - 2, p)) % p
    return A, B


def _resolve_curve(f, g, p):
    """Analyze y² = x³ + fx + g over F_p; route to the right view.

    Returns ('curve', a, node) when the class is in the tables (the full
    curve page), ('eq', data) when only the computed invariants exist,
    or ('singular', None)."""
    if (4 * pow(f, 3, p) + 27 * pow(g, 2, p)) % p == 0:
        return ("singular", None)
    data = ec_look_up((f % p, g % p), p)
    a = data["frob_tr"]
    if ap_in_pc_data((a, p)):
        cls = _load_ec_class(a, p)
        key = data["j"] if a != 0 else (data["j"], data["s"])
        qf = (cls.js_to_qf or {}).get(key)
        if qf is not None:
            return ("curve", a, cls.qfs_ordered.index(qf))
    return ("eq", data)


def _crumbs(*parts):
    """Breadcrumb line: ('label', href-or-None), current view last with None."""
    bits = [f"[{label}]({href})" if href else f"**{label}**"
            for label, href in parts]
    st.markdown(" › ".join(bits))


# Friendly table names: Nadir's column keys are coding conveniences.
COL_RENAME = {
    "j_inv": "j-invariant",
    "ec_invs": "signature",
    "EC_coefs": "curve (f, g)",
    "qf_coefs": "lattice (a, b, c)",
    "endo_disc": "End-ring disc",
    "endo_cond": "End-ring conductor",
    "endo_cocond": "index [End(E) : ℤ[φ]]",
    "frobmat": "Frobenius matrix",
    "tau_s": "τ",
    "tau_xy": "τ = (x, y)",
}

LEGEND = r"""
**Reading the table.** Each row pairs an isomorphism class of **CM lattices**
with an isomorphism class of **elliptic curves over $\mathbb{F}_p$**; the two
sides are matched by an equivalence of categories.

* **signature** — determines the $\mathbb{F}_p$-isomorphism class. For ordinary
  classes it is just the $j$-invariant; for supersingular classes it is a pair
  $(j, \pm 1)$ recording a quadratic-twist sign.
* **j-invariant** — the $j$-invariant of the curve.
* **curve (f, g)** — coefficients of the short Weierstrass model
  $y^2 = x^3 + fx + g$.
* **lattice (a, b, c)** — the reduced binary quadratic form
  $ax^2 + bxy + cy^2$ that pins down the lattice class.
* **End-ring disc / conductor** — discriminant and conductor of the
  endomorphism ring of that row.
* **index [End(E) : ℤ[φ]]** — how far the Frobenius order sits inside the full
  endomorphism ring.
* **Frobenius matrix** — the action of the Frobenius lift on the lattice in the
  ordered basis $1, \tau$.
* **τ** — the CM point of the lattice class in the upper half-plane.
"""


# ── the four views ────────────────────────────────────────────────────────────

def entry_view():
    st.header("Explorer")
    st.markdown(
        "Everything on this site funnels through one picture: isogeny classes "
        "of elliptic curves over $\\mathbb{F}_p$ seen through their CM "
        "lattices. Pick a door — a **prime** (the Frobenius view) or a "
        "**discriminant** (the characteristic-0 view). Every view you land on "
        "is a shareable URL."
    )
    col_p, col_d = st.columns(2, gap="large")

    with col_p:
        st.subheader("Start from a prime")
        p_raw = st.number_input("prime p", min_value=5, max_value=P_MAX - 1,
                                value=101, step=2,
                                help=f"Curve tables cover 99.86% of classes for "
                                     f"p < {P_MAX}; the few gaps render "
                                     "structure-only.")
        p = int(p_raw)
        if not primeQ(p):
            below = primesBetween(2, p)
            p = int(below[-1]) if below else 5
            st.caption(f"{int(p_raw)} isn't prime — showing the previous prime, **p = {p}**.")
        st.markdown(f"Each admissible trace $a$ (with $a^2 < 4p$) is one "
                    f"isogeny class over $\\mathbb{{F}}_{{{p}}}$. Click a stem.")
        components.html(explorer_viz.hasse_picker_html(p, STORE_P_MAX),
                        height=400, scrolling=False)

        st.markdown("**…or start from an equation** (over the same $p$)")
        form = st.radio("Weierstrass form",
                        ["y² = x³ + fx + g",
                         "y² + a₁xy + a₃y = x³ + a₂x² + a₄x + a₆"],
                        label_visibility="collapsed", key="eq_form")
        if form.startswith("y² ="):
            cf, cg = st.columns(2)
            f_in = int(cf.number_input("f", value=1, step=1, key="eq_f"))
            g_in = int(cg.number_input("g", value=3, step=1, key="eq_g"))
        else:
            cols = st.columns(5)
            a1, a2, a3, a4, a6 = (
                int(c.number_input(lbl, value=v, step=1, key=f"eq_{lbl}"))
                for c, lbl, v in zip(cols, ("a₁", "a₂", "a₃", "a₄", "a₆"),
                                     (0, -1, 1, -10, -20)))
        if st.button("Analyze the curve →", width="stretch"):
            if form.startswith("y² ="):
                A, B = f_in % p, g_in % p
            else:
                A, B = _long_to_short(a1, a2, a3, a4, a6, p)
            hit = _resolve_curve(A, B, p)
            if hit[0] == "singular":
                st.error("That curve is singular over "
                         f"$\\mathbb{{F}}_{{{p}}}$ (Δ ≡ 0 mod p) — not an "
                         "elliptic curve. Nudge a coefficient.")
            else:
                st.query_params.clear()
                if hit[0] == "curve":
                    _, a_tr, node = hit
                    st.query_params.update(
                        {"a": str(a_tr), "p": str(p), "node": str(node)})
                else:
                    st.query_params.update(
                        {"p": str(p), "f": str(A), "g": str(B)})
                st.rerun()

    with col_d:
        st.subheader("Start from a discriminant")
        d_raw = st.number_input("discriminant d", max_value=-3, value=-368, step=1,
                                help="A negative discriminant: d ≡ 0 or 1 (mod 4).")
        d = int(d_raw)
        if not _valid_disc(d):
            st.caption(f"{d} ≡ {d % 4} (mod 4) isn't a discriminant — "
                       f"nearest valid is **{d - (d % 4 - 1) if (d % 4) in (2, 3) else d}**.")
        if st.button("Open the discriminant view →", width="stretch"):
            while not _valid_disc(d):
                d += 1
            st.query_params.clear()
            st.query_params["d"] = str(d)
            st.rerun()
        st.markdown(
            "The discriminant view shows the pure lattice side: the class "
            "group, the ℓ-isogeny structure of all orders between "
            "$\\mathbb{Z}[\\sqrt{d}\\,]$-level and the maximal order — no "
            "prime chosen, no distinguished Frobenius — and then lists every "
            "$(a, p)$ where that structure is realized by curves."
        )


def disc_view(d: int):
    _crumbs(("⌂ Explorer", "?"), (f"discriminant {d}", None))
    if not _valid_disc(d):
        st.error(f"{d} is not a negative discriminant (need d < 0, d ≡ 0 or 1 mod 4).")
        return
    cls = _load_qf_class(d)
    n = len(cls.qfs_all)
    if n > GUARDS["volcano_nodes"]:
        st.warning(f"This class has {n} lattice classes — beyond the drawing "
                   f"guard ({GUARDS['volcano_nodes']}). Try a smaller |d|.")
        return
    d0, c = cls.field_disc, cls.cond
    st.header(f"Discriminant {d}")
    st.markdown(f"$d = {d} = {d0} \\cdot {c}^2$ &nbsp;·&nbsp; field disc ${d0}$, "
                f"conductor ${c}$ &nbsp;·&nbsp; **{n}** lattice classes across "
                f"all conductor levels.")

    st.markdown("##### The ℓ-isogeny structure (characteristic 0)")
    components.html(explorer_viz.isogeny_graph_html(_descs_for(cls)),
                    height=740, scrolling=False)

    st.markdown("##### Where it lives over finite fields")
    st.markdown(
        "Pairs $(a, p)$ with $a^2 - 4p = d\\,m^2$: for $m = 1$ the Frobenius "
        "order has exactly this discriminant (primitive); larger $m$ puts $d$ "
        "higher up the volcano (imprimitive)."
    )
    pairs = disc_to_aps(d)
    if not pairs:
        st.info(f"No prime p < {P_MAX} realizes this discriminant.")
        return
    shown = pairs[:80]
    rows = ["| a | p | m | curve tables | open |", "|---|---|---|---|---|"]
    for e in shown:
        a_, p_, m_ = e["a"], e["p"], e["m"]
        instore = "✓" if e["in_store"] else "structure only"
        rows.append(f"| {a_} | {p_} | {m_} | {instore} | "
                    f"[({a_}, {p_}) →](?a={a_}&p={p_}) |")
    st.markdown("\n".join(rows))
    if len(pairs) > len(shown):
        st.caption(f"…and {len(pairs) - len(shown)} more pairs up to p < {P_MAX}.")


def class_view(a: int, p: int):
    d = a * a - 4 * p
    _crumbs(("⌂ Explorer", "?"), (f"disc {d}", f"?d={d}"),
            (f"(a, p) = ({a}, {p})", None))
    cls = _load_ec_class(a, p)
    n = len(cls.qfs_all)
    if n > GUARDS["volcano_nodes"]:
        st.warning(f"{n} classes — beyond the drawing guard.")
        return
    has_curves = cls.js_to_qf is not None
    N = p + 1 - a
    st.header(f"Isogeny class (a, p) = ({a}, {p})")
    chi = (f"x^2 + {p}" if a == 0
           else f"x^2 + {-a}x + {p}" if a < 0
           else f"x^2 - {a}x + {p}")
    facts = (f"$\\chi(x) = {chi}$ &nbsp;·&nbsp; "
             f"$\\#E(\\mathbb{{F}}_p) = {N}$ &nbsp;·&nbsp; disc ${d}$ "
             f"(field disc ${cls.field_disc}$, cond ${cls.cond}$) &nbsp;·&nbsp; "
             f"**{n}** curves")
    if a == 0:
        facts += " &nbsp;·&nbsp; supersingular"
    st.markdown(facts)
    if not has_curves:
        st.info("This class is one of the few (0.14%) whose curve tables are "
                "still being computed, so this view shows the lattice side "
                "only — the graph structure is exact.")

    link_base = f"?a={a}&p={p}"
    st.markdown("##### The ℓ-isogeny graph")
    components.html(
        explorer_viz.isogeny_graph_html(_descs_for(cls), link_base=link_base),
        height=740, scrolling=False)

    # FD picture — points ordered like qfs_ordered so &node= indices agree
    st.markdown("##### The lattice classes as CM points")
    df = cls.ecqf_df() if has_curves else None
    colors = row_colors(n)
    qf_to_color, qf_to_js, qf_to_fg = {}, {}, {}
    if df is not None:
        for i, (_, row) in enumerate(df.iterrows()):
            qf = tuple(row["qf_coefs"])
            qf_to_color.setdefault(qf, colors[i])
            qf_to_js.setdefault(qf, []).append(row["j_inv"])
            qf_to_fg.setdefault(qf, tuple(int(v) for v in row["EC_coefs"]))
    pts = []
    for i, qf in enumerate(cls.qfs_ordered):
        x, y = abc_to_tau(qf)
        js = qf_to_js.get(qf)
        label = f"j={js[0]}" + ("…" if js and len(js) > 1 else "") if js else str(qf)
        sub = (f"⟨{', '.join(str(v) for v in qf)}⟩ · End disc {cls.endo_disc_dict[qf]}"
               f" · τ = {abc_to_tau_str(qf)}")
        if js:
            sub += " · j = " + ", ".join(str(j) for j in js)
        if qf in qf_to_fg:
            fg = tuple(_centered(v, p) for v in qf_to_fg[qf])
            sub += f" · {ec_eq_str_base(fg)}"
        pts.append({"x": float(round(x, 5)), "y": float(round(y, 5)),
                    "color": qf_to_color.get(qf, "#4da3d8"),
                    "label": label, "sub": sub})
    components.html(explorer_viz.fd_points_html(pts, link_base=link_base),
                    height=560, scrolling=False)

    # the table, as a drawer
    if df is not None:
        with st.expander(f"data table — all {n} curves"):
            discs = sorted(set(df["endo_disc"]))
            c1, c2 = st.columns([1, 2])
            with c1:
                ring = st.selectbox("filter by endomorphism ring",
                                    ["all"] + [str(v) for v in discs],
                                    key="tbl_ring")
            with c2:
                cols = st.multiselect("columns", list(COL_RENAME),
                                      default=["j_inv", "EC_coefs", "qf_coefs",
                                               "endo_disc", "endo_cond", "tau_s"],
                                      format_func=lambda k: COL_RENAME[k],
                                      key="tbl_cols")
            view = df if ring == "all" else df[df["endo_disc"] == int(ring)]
            if cols:
                st.dataframe(view[cols].rename(columns=COL_RENAME),
                             width="stretch", hide_index=True)
            st.markdown(LEGEND)


def curve_view(a: int, p: int, node: int):
    d = a * a - 4 * p
    cls = _load_ec_class(a, p)
    if not (0 <= node < len(cls.qfs_ordered)):
        st.error(f"node {node} is out of range for this class.")
        return
    qf = cls.qfs_ordered[node]
    has_curves = cls.js_to_qf is not None

    js, fg = [], None
    if has_curves:
        for jsig, q in cls.js_to_qf.items():
            if q == qf:
                j = jsig[0] if isinstance(jsig, tuple) else jsig
                js.append(j)
                if fg is None:
                    fg = tuple(int(v) for v in cls.js_to_models[jsig])
    title = f"j = {js[0]}" if js else f"lattice ⟨{', '.join(str(v) for v in qf)}⟩"
    _crumbs(("⌂ Explorer", "?"), (f"disc {d}", f"?d={d}"),
            (f"(a, p) = ({a}, {p})", f"?a={a}&p={p}"), (title, None))
    st.header(f"Curve view — {title} over 𝔽_{p}")

    frm = cls.qf_to_frob_mats[qf]
    N = p + 1 - a
    gendic = frob_to_mw_gens(frm, 1)
    gens = sorted(gendic, key=lambda g: gendic[g], reverse=True)
    iso = " ⊕ ".join(f"ℤ/{gendic[g]}" for g in gens) if gens else "trivial"

    lines = []
    if fg is not None:
        fgc = tuple(_centered(v, p) for v in fg)
        lines.append(f"**Equation** &nbsp; {ec_eq_str_base(fgc)}")
        lines.append(f"**j-invariant{'s' if len(js) > 1 else ''}** &nbsp; "
                     + ", ".join(str(j) for j in js))
    lines.append(f"**Lattice class** &nbsp; ⟨{', '.join(str(v) for v in qf)}⟩ "
                 f"&nbsp;·&nbsp; τ = {abc_to_tau_str(qf)}")
    lines.append(f"**Endomorphism ring** &nbsp; disc {cls.endo_disc_dict[qf]}, "
                 f"conductor {cls.endo_cond_dict[qf]}")
    lines.append(f"**Frobenius on ⟨1, τ⟩** &nbsp; matrix {frm.vec}")
    lines.append(f"**Group structure** &nbsp; E(𝔽_{p}) ≅ {iso} &nbsp; "
                 f"(order {N})")
    if gens:
        lines.append("**Generators** &nbsp; " + ", ".join(
            f"{g} of order {gendic[g]}" for g in gens)
            + " &nbsp;·&nbsp; coordinates w.r.t. the basis 1, τ on ℂ/Λ")
    st.markdown("  \n".join(lines))

    n_pts = N if N <= GUARDS["fd_points"] else 0
    components.html(explorer_viz.curve_torus_html(qf, a, p, frm.vec, n_pts),
                    height=480, scrolling=False)

    # walk the graph: this node's neighbours, per degree
    st.markdown("##### Walk the isogeny graph")
    descs = _descs_for(cls)
    any_nbrs = False
    for l, desc in descs.items():
        nbrs = []
        for e in desc["edges"]:
            other = e["t"] if e["s"] == node else (e["s"] if e["t"] == node else None)
            if other is None:
                continue
            nd = desc["nodes"][other]
            up_down = ("horizontal" if e["kind"] == "h" else
                       ("up ↑" if nd["endo_cond"] < cls.endo_cond_dict[qf] else "down ↓"))
            lbl = (f"j = {nd['curves'][0]['j']}" if nd["curves"]
                   else f"⟨{', '.join(str(v) for v in nd['qf'])}⟩")
            nbrs.append(f"[{lbl} ({up_down})](?a={a}&p={p}&node={other})")
        if nbrs:
            any_nbrs = True
            st.markdown(f"**ℓ = {l}:** &nbsp; " + " · ".join(nbrs))
    if not any_nbrs:
        st.caption("No ℓ-isogenies for the offered degrees (they are inert here).")


def eq_view(p: int, f: int, g: int):
    """A curve given by equation, outside the precomputed tables (or a
    shared ?p&f&g URL — which redirects to the full curve view when the
    class turns out to be covered)."""
    hit = _resolve_curve(f, g, p)
    if hit[0] == "curve":
        _, a_tr, node = hit
        st.query_params.clear()
        st.query_params.update({"a": str(a_tr), "p": str(p), "node": str(node)})
        st.rerun()
    fgc = tuple(_centered(v, p) for v in (f, g))
    _crumbs(("⌂ Explorer", "?"), (f"{ec_eq_str_base(fgc)} over 𝔽_{p}", None))
    if hit[0] == "singular":
        st.error(f"y² = x³ + {fgc[0]}x + {fgc[1]} is singular over "
                 f"$\\mathbb{{F}}_{{{p}}}$ (Δ ≡ 0 mod p) — not an elliptic "
                 "curve.")
        return
    data = hit[1]
    a = data["frob_tr"]
    N = p + 1 - a
    d = data["frob_disc"]
    st.header(f"Curve view — {ec_eq_str_base(fgc)} over 𝔽_{p}")
    lines = [
        f"**j-invariant** &nbsp; {data['j']}",
        f"**Cardinality** &nbsp; #E(𝔽_{p}) = {N}",
        f"**Trace of Frobenius** &nbsp; a = {a}",
        f"**Frobenius discriminant** &nbsp; {d} "
        f"(field disc {data['endo_fun_disc']}, cond {data['frob_cond']})",
    ]
    if a == 0:
        lines.append(f"**Supersingular** &nbsp; signature "
                     f"({data['j']}, {data['s']:+d}) — this pins down the "
                     f"𝔽_{p}-isomorphism class")
    st.markdown("  \n".join(lines))
    st.info("This class is one of the few (0.14%) whose curve tables are "
            "still being computed — the invariants above hold regardless, "
            f"and the class's exact lattice structure is at the "
            f"[isogeny-class view](?a={a}&p={p}).")
    if N <= 2000:
        with st.expander(f"the {N} points of E(𝔽_{p})"):
            sq = {}
            for y in range(p):
                sq.setdefault(y * y % p, []).append(y)
            pts = ["𝒪 (at infinity)"]
            for x in range(p):
                rhs = (x ** 3 + f * x + g) % p
                for y in sorted(set(sq.get(rhs, []))):
                    pts.append(f"({x}, {y})")
            st.markdown(" · ".join(pts))


# ── router ────────────────────────────────────────────────────────────────────

_a, _p, _d, _node, _f, _g = (_qp_int(k) for k in ("a", "p", "d", "node", "f", "g"))

if _a is not None and _p is not None and primeQ(_p) and _a * _a < 4 * _p \
        and (_a == 0 or _a % _p != 0) and _p < P_MAX:
    if _node is not None:
        curve_view(_a, _p, _node)
    else:
        class_view(_a, _p)
elif _p is not None and _f is not None and _g is not None \
        and primeQ(_p) and 5 <= _p < P_MAX:
    eq_view(_p, _f % _p, _g % _p)
elif _d is not None:
    disc_view(_d)
else:
    if _a is not None or (_p is not None and _f is None):
        st.error("That (a, p) isn't admissible — need p prime, a² < 4p, p ∤ a "
                 f"(or a = 0), p < {P_MAX}.")
    entry_view()
