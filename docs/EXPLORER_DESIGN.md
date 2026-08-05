# The Explorer, redesigned for a static site

Design for the Explorer after the move off Streamlit. This is a **redesign, not
a port** — the goal is what the section should be when no framework is shaping
it, keeping the feature set but not the layout.

Companion to `docs/DEPLOY.md` (current hosting) and `netlify-proto/` (the
working proof that the class view needs no server).

## Why redesign rather than clone

Measured before deciding, so the reasoning is on the record:

- The full precompute over all **117,155** classes is **~190 MB gzipped**
  (mean 1.6 KB/class), sharded one file per prime: **1,026 files, ~181 KB each**.
  Nothing here needs a server at visit time.
- `explorer_viz`'s applets are already data-driven (`const DATA = <json>` inside
  an IIFE), so the drawing code ports **byte-identical**. See
  `netlify-proto/build.py:split_applet`.

Three constraints in the current design are Streamlit artifacts, not choices:

| Artifact | Consequence | After |
|---|---|---|
| Each applet is a sandboxed iframe | `navParent` exists purely to defeat the sandbox; applets cannot talk to each other | Applets are page elements; they share state |
| Iframes can't self-size | Every applet has a hardcoded `height=` | Applets size to content |
| Query-param routing on a JS app | No class has a real address; none are indexable | Real URLs, one page per prime prerendered |

## Decisions

Settled with NH, 2026-08-05:

- **Landing page is just the search box.** LMFDB idiom — you came to look
  something up. No featured example, no overview panel.
- **Downloads are per-object.** A curve or a class or a prime, on its own page.
  Anyone wanting the whole dataset clones the repo; there is no bulk endpoint
  and therefore no schema to keep stable for outsiders.
- **Prerender prime pages only** (1,026 files, class facts inline so the content
  is indexed). Class and curve pages render client-side from the shards.

## Object model and URLs

A discriminant and a prime are the two axes; a class sits where they cross; a
curve sits inside a class. The current breadcrumbs already say this
(`Explorer › disc -368 › (a,p) = (6,101)`) — the URLs should too.

```
/explorer/                 search
/explorer/p/101/           F_101 — every class over it          [prerendered]
/explorer/d/-368/          the discriminant + every (a,p) realizing it
/explorer/c/101/6/         the isogeny class (a, p) = (6, 101)
/explorer/c/101/6/j30/     one curve, keyed by j-invariant
```

**Curves are keyed by j, not by node index.** Today a curve is `&node=4`, a
position in `cls.qfs_ordered`. That ordering is an implementation detail: if it
ever changes, every saved or cited link silently resolves to a *different*
curve. `j` is intrinsic to the object, so the address survives any regeneration
of the data. This is the one change that must not be skipped — the site is
meant to be cited.

Where a class has several curves sharing a j (quadratic twists), the page shows
them together; the twist is not a separate address.

## Pages

### `/explorer/` — search

A single field, and nothing else. It parses:

| Input | Resolves to |
|---|---|
| `101` | `/explorer/p/101/` |
| `-368` | `/explorer/d/-368/` |
| `(6,101)`, `6, 101`, `a=6 p=101` | `/explorer/c/101/6/` |
| `j = 30 mod 101` | the curve with that j over F_101 |
| `y^2 = x^3 + 3x + 2 mod 101` | resolve the curve, then its class |
| `y^2 + xy = x^3 - 10x - 20 mod 101` | long Weierstrass form, same |

Ambiguity rule: a bare positive integer is a prime, a bare negative integer is a
discriminant. Non-prime positive input offers the nearest prime rather than
erroring. Input that names a class outside the covered range says so plainly and
states the bound (`p < 8192`, `|d| <= 32768`) rather than failing silently.

The equation forms deserve top billing here. They are the least prominent thing
on the current entry page and are probably what a researcher most often has in
hand: they have a curve, and want to know its class, its CM order, its lattice.

### `/explorer/p/101/` — the prime slice  *(prerendered)*

The Hasse picker, plus a table of every class over F_p with facts inline
(trace, #E, disc, field disc, conductor, curve count). The inline table is what
makes prerendering worth doing: this is the page a search engine reads, and it
carries real content rather than an empty shell.

Download: JSON and CSV for the whole slice.

### `/explorer/c/101/6/` — the isogeny class

Header facts as today. Then the three views of the same lattice classes —
**volcano graph, CM points, table — sharing one selection.** Hover or select in
any of them and the other two follow.

This is the payoff of leaving the iframes behind, and the strongest reason the
port improves the site rather than merely relocating it: the volcano and the CM
picture are two pictures of the same twelve objects, and right now they cannot
be made to say so.

Requires a small addition to `explorer_viz`: each applet emits its selection and
accepts one back. Guard the emit with optional chaining (`window.onSelect?.(i)`)
so the Streamlit build is unaffected while both live side by side.

Download: JSON for the class.

### `/explorer/c/101/6/j30/` — one curve

The curve view as it stands (torus picture, Frobenius matrix, the lattice and
its endomorphism ring), addressed stably.

Download: JSON for the curve.

### `/explorer/d/-368/` — the discriminant

The characteristic-0 lattice side, and every `(a, p)` realizing it.

**Keep the `|d| <= 4·P_MAX` bound** introduced in `dad7a8c`. On the static site
it stops being a safety guard — there is no process to stall — and becomes a
build-scope statement: 16,384 valid discriminants, all precomputable.

## Build and deploy

The heavy generation runs on NH's laptop, never in Netlify's build (which would
time out):

1. `build.py` precomputes every class → per-prime JSON shards, plus the
   prerendered prime pages.
2. Roughly **8 hours single-threaded, ~1 hour across 10 cores.** Run
   deliberately, when the data changes — not per deploy.
3. Publish by uploading the built folder, so ~190 MB never enters git.

Generated output stays gitignored, as in `netlify-proto/`.

## Dropped, and why

- **The two-door entry view** — replaced by search. Two doors was `st.columns`.
- **`?p&f&g` as its own view** — folds into search; it was a separate route only
  because there was nowhere else to put an equation.
- **The table's column multiselect** — a Streamlit affordance. Ship sensible
  columns with a "show everything" toggle.
- **`&node=` addressing** — superseded by `j` (see above).

## Open

- What the curve page shows when a class has no curve tables (the 0.14%).
- Whether `/explorer/d/` pages are worth prerendering too (16,384 of them —
  cheap next to 117k, and the discriminant is the more citable object).
- Whether the search box should accept an LMFDB label and link out.
