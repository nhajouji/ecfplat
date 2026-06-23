# Handoff: interactive fundamental-domain picture for the Isogeny Class page

## Goal
On the **Isogeny Class** page ([pages/1_Isogeny_Class.py](../pages/1_Isogeny_Class.py)),
add an interactive picture beside the lattice-class dataframe. Each row of the
dataframe is an isomorphism class of CM lattices = a reduced binary quadratic
form of discriminant `d = a² − 4p`. Draw one **point per row** in the standard
fundamental domain of `SL(2,Z)\H`, linked to the table:

- **Click a point →** highlight the corresponding dataframe row (and vice versa).
- **Open the curve:** for v1, single-click selects + the existing "View this
  curve in EC Search →" button navigates. (True double-click→navigate is a
  later enhancement — see "Future".)

## Toolkit decision (already made)
- **Plotly** via `st.plotly_chart(fig, on_select="rerun", selection_mode="points")`
  is the standard for interactive applets going forward. Native Streamlit event
  support, no new dependency, and it scales to interactive 3D (`Scatter3d`,
  `Surface`) and frame-based animation for future pages.
- A **custom React component** is the eventual escape hatch for shader-quality /
  real-time / 3D-with-click-events graphics. Not needed for this feature.
- Altair was rejected (2D only, weak animation).

## The math is already in the repo — reuse it
- `abc_to_tau((A,B,C)) -> np.array([Re τ, Im τ])` — [pycode/ecqf_tools.py:52](../pycode/ecqf_tools.py)
  maps a quadratic form to its CM point τ in the fundamental domain.
- Static matplotlib reference to port: `uhmodgam_fd_plot`, `qfs_uhfd_plot`,
  `qfs_uhfd_plot_w_highlights` — [pycode/graphic_tools.py:213](../pycode/graphic_tools.py).
  These already draw the FD (gray rectangle Re∈[−½,½] minus the unit disc) and
  the τ points, with a red highlight.
- Per-row data: `df = isoclass.ecqf_df()` has columns `ec_invs`, `EC_coefs`
  (= (f, g)), `qf_coefs` (= the form (A,B,C)). Canonical order is
  `isoclass.qfs_ordered`. `ecqf_df` is defined at
  [pycode/ecqf_tools.py:622](../pycode/ecqf_tools.py).

So the static picture is essentially `[abc_to_tau(qf) for qf in <df qf order>]`.
The new work is interactivity + bidirectional linking with the table.

## Plan

### 1. Build the Plotly figure
- Points: `tau = abc_to_tau(row["qf_coefs"])` for each df row, **in df-row order**
  (so the Plotly point index equals the df row index — critical, see Gotchas).
- Fundamental-domain boundary as Plotly shapes: the two vertical segments
  `Re τ = ±½` and the arc `|τ| = 1` between `e^{iπ/3}` and `e^{i2π/3}`. Optional
  faint fill for the FD interior.
- Hover text per point: row index, form `(A,B,C)`, and `ec_invs` (let the reader
  see which curve a point is).
- Reasonable y-cap so tall cusps don't squash the picture; equal aspect.

### 2. Click point → highlight row (and row → highlight point)
- Drive everything off a single `st.session_state.selected_row`.
- `event = st.plotly_chart(fig, on_select="rerun", selection_mode="points")`;
  read the selected point index from `event.selection` → set `selected_row`.
- The dataframe already supports `st.dataframe(..., on_select="rerun",
  selection_mode="single-row")` ([pages/1_Isogeny_Class.py:70](../pages/1_Isogeny_Class.py));
  read its selection too. Whichever was touched this run updates `selected_row`.
- Re-render: redraw the Plotly figure with the selected τ marked (red/larger),
  and highlight the table row (pandas `Styler.apply` row highlight, rendered via
  `st.dataframe`/`st.write` of the styled frame — note: `st.dataframe`'s own
  selection + a Styler highlight can coexist; verify the styling shows).

### 3. Navigation (v1)
- Reuse the existing mechanism verbatim:
  `st.session_state.ec_prefill = {"f": f, "g": g, "p": p}` then
  `st.switch_page("pages/2_EC_Search.py")` — see
  [pages/1_Isogeny_Class.py:89](../pages/1_Isogeny_Class.py). Keep the button;
  it now also works when the selection came from clicking a point.

### 4. Layout
- Put the picture next to the table (e.g. `st.columns([1,1])` inside the
  existing "Isogeny class table" tab), or a sub-tab. Keep the existing table +
  button behavior intact.

## Gotchas
- **Index alignment is everything.** Plot points in the exact order of the df
  rows so `point index == df row index`. `ecqf_df()` row order vs
  `qfs_ordered` — confirm they match (they should), or build the τ list from the
  df's `qf_coefs` column directly to be safe.
- **Orientation/conjugation convention.** The curve↔lattice bijection has a
  global complex-conjugation freedom; make sure the τ shown for a row maps to the
  same curve the row's `EC_coefs` point to (sanity-check one class end to end:
  click point → row → "View curve" → correct (f,g,p)).
- **Crowding for large h(d).** Points cluster near the corners/cusp; rely on
  hover labels and allow zoom. Cap Im τ for the view.
- **Notation:** the form's leading coefficient `A` vs the Frobenius trace `a`.
- **Plotly selection payload shape:** `st.plotly_chart` returns selection as
  `event.selection["points"]` (list of dicts with `point_index`/`point_number`);
  confirm the exact key on the installed Plotly/Streamlit and map to the row.
  (Pinned env: Streamlit 1.58, Python 3.13 — see requirements.txt.)

## Future (not v1)
- **Double-click → navigate:** needs either `streamlit-plotly-events` (lightly
  maintained; check compat) or a small custom React component that emits
  `plotly_doubleclick` and calls back. The custom-component route doubles as the
  foundation for shader/3D graphics later.
- **3D / animation** applets (tori, Frobenius orbits in motion): Plotly first;
  custom component when it needs to be real-time or shader-quality.

## Definition of done (v1)
Picture beside the table; clicking a point highlights its row and clicking a row
highlights its point; the "View this curve" button opens the right EC page for
the selected class; verified on a couple of `(a,p)` with h(d) > 1 (e.g. one small
and one with a dozen+ classes) in the running app.
