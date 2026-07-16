import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "pycode"))

import streamlit as st
import streamlit.components.v1 as components

import blog_viz

IMG = Path(__file__).parent / "blog_img"

# ──────────────────────────────────────────────────────────────────────────────
# Post registry.  Add a dict here (newest first) and write its render function;
# the index and the ?post=<slug> routing pick it up automatically.
# ──────────────────────────────────────────────────────────────────────────────


def _img(name, caption):
    st.image(str(IMG / name), caption=caption, width="stretch")


# ── Post 1 — sum of two squares ───────────────────────────────────────────────
def render_sum_of_squares():
    st.markdown("# Minkowski and Frobenius")
    st.caption("How a 350-year-old theorem of Fermat — plus a little geometry of "
               "numbers — says something concrete about the Bitcoin curve. · July 2026")

    st.markdown(r"""
Someone once asked, on a Bitcoin forum, what the curve **secp256k1** *looks like*.
It's a fair question with an awkward answer: the curve lives over a field with
about $10^{77}$ elements, so there is no honest picture to draw. But we can ask a
better question — **what can we learn about it using down-to-earth, visualizable
methods?**

secp256k1 is $y^2 = x^3 + 7$, which has $j$-invariant $0$, and its prime $p$ is
$\equiv 1 \pmod 3$, so the curve is *ordinary*. To pin down its arithmetic (the
number of points, the trace of Frobenius) it turns out we need to write

$$p = x^2 + 3y^2.$$

That has a friendlier cousin that everyone can picture, and it's where I want to
start: for a prime $p \equiv 1 \pmod 4$, **Fermat** tells us

$$p = a^2 + b^2.$$

Fermat's theorem promises the two squares *exist*. The interesting question is
sharper: given $p$, how do you **find** $a$ and $b$? And it turns out the whole
thing is secretly a problem about a lattice.
""")

    st.markdown("### The answer is a shortest vector")
    st.markdown(r"""
Fix an integer $r$ with $r^2 \equiv -1 \pmod p$ (this exists exactly when
$p \equiv 1 \pmod 4$). Now look at the lattice

$$\Lambda = \{(a,b)\in\mathbb{Z}^2 : a \equiv r\,b \pmod p\}.$$

Two things are true and easy to check: $\Lambda$ has one point per $p$ cells of
the integer grid (covolume $p$), and **every** point of $\Lambda$ has
$a^2 + b^2 \equiv 0 \pmod p$. So if we can find a *short enough* nonzero point —
one with $a^2 + b^2 < 2p$ — then its norm is a positive multiple of $p$ below
$2p$, which leaves only one option: $a^2 + b^2 = p$.

Minkowski's theorem guarantees such a short point exists. But I don't want to
just know it exists — I want to **find it**. The problem has become purely
geometric: *find the shortest nonzero vector of $\Lambda$.*
""")
    _img("gauss_minkowski.png",
         "The lattice Λ (index p in ℤ[i]) with the disc a²+b² < 2p. Its shortest "
         "vector is the sum-of-two-squares answer; here 113 = (−8)² + 7².")

    st.markdown("### Finding it: a descent you can watch")
    st.markdown(r"""
Here is the trick. Read the pair $(a,b)$ as a **Gaussian integer** $a + bi$, and
run the **Euclidean algorithm in $\mathbb{Z}[i]$**. Every step is a single,
completely visual move: *find the nearest multiple of the current divisor* — the
nearest point of a rotated, scaled copy of the square lattice — and take the
remainder, which is strictly shorter. The "how many times $p$" left over (the
cofactor $N/p$) strictly decreases, so by plain well-ordering the process must
stop — and it stops exactly at the shortest vector, whose norm is $p$.

Pick a prime below and step through it yourself. Watch the red vector spiral in
onto the circle $|z| = \sqrt p$, where it lands on the answer.
""")
    components.html(blog_viz.sos_descent_html(), height=820, scrolling=True)

    st.markdown("### Why it never gets stuck — and the payoff for secp256k1")
    st.markdown(r"""
The descent can't stall, and the reason is a small piece of magic: $\mathbb{Z}[i]$
is **norm-Euclidean**, so every division genuinely produces a shorter remainder.
More broadly, *this trick works for exactly those curves whose complex
multiplication is by a norm-Euclidean order.*

For the Bitcoin curve, $j = 0$, the relevant order isn't the Gaussian integers
but the **Eisenstein integers** $\mathbb{Z}[\omega]$ — the hexagonal lattice —
which is also norm-Euclidean. The same descent finds $p = x^2 + 3y^2$, and from
that single representation the **six candidate traces of Frobenius** of
$y^2 = x^3 + 7$ fall out (rotate by the six units of $\mathbb{Z}[\omega]$). On the
actual 256-bit secp256k1 prime this runs in about a millisecond.
""")
    _img("eisenstein_minkowski.png",
         "The j = 0 story: the hexagonal Eisenstein lattice. Its six shortest "
         "vectors give p = x² + 3y² and the six candidate Frobenius traces.")

    st.markdown("### When is it slowest? (a Fibonacci surprise)")
    st.markdown(r"""
There's a lovely coda. For ordinary integers, the Euclidean algorithm is slowest
on consecutive **Fibonacci numbers** — the golden ratio. Our version inherits
that: the slowest primes for the real algorithm are the **Fibonacci primes**
$p = F_{2k+1}$, and thanks to the identity $F_{2k+1} = F_k^2 + F_{k+1}^2$ the
answer it grinds out is a pair of consecutive Fibonacci squares. (Try $28657$ in
the widget: it returns $89^2 + 144^2$.)

But the *complex* algorithm uses nearest-integer division, and its worst case is
a different animal: the **silver ratio** $1+\sqrt2$ instead of the golden ratio.
Swapping the real Euclid for the complex one swaps golden for silver — a small,
concrete instance of the classical fact that the nearest-integer continued
fraction beats the ordinary one.
""")
    _img("worstcase.png",
         "Left: the plain algorithm is slowest at Fibonacci primes, riding the "
         "golden-ratio ceiling. Right: the two algorithms' worst cases point "
         "along the golden vs silver directions.")

    st.markdown("### The bottom line")
    st.markdown(r"""
The naive ways to find $a, b$ — searching pairs, or counting points — cost about
$p$ operations, which for a 256-bit prime is hopeless. The lattice descent costs
a number of steps polynomial in the *number of digits* of $p$: milliseconds, not
millennia. (Number theorists will recognize this as **Cornacchia's algorithm**,
rederived here as pure geometry of numbers.) So we can't draw secp256k1 — but we
can hold its arithmetic in our hands, and every step of getting there is a
picture.
""")


# ── Post 2 — the representation-theory story ──────────────────────────────────
def render_three_languages():
    st.markdown("# One square, three languages")
    st.caption("The same class group, read as geometry, as representation theory, "
               "and as particle physics — and a ± ambiguity that turns out to be "
               "the same ghost in all three. · July 2026")

    st.markdown(r"""
Some years ago I learned to draw a picture that, it turned out, lives in two
completely different corners of mathematics at once. Take a finite group $G$ and
a chosen representation $\rho$. Draw one dot per irreducible representation of
$G$, and an arrow from $v$ to $w$ whenever $w$ shows up inside $v \otimes \rho$.
That graph is called a **fusion graph** (or McKay graph). For the finite groups
$G \subset \mathrm{SL}_2(\mathbb{C})$, with $\rho$ the defining $2$-dimensional
representation, these graphs are exactly the **affine ADE diagrams** — the very
same diagrams that classify how the surface singularity $\mathbb{C}^2/G$ gets
resolved. Two different constructions, one list of pictures. That is the McKay
correspondence.

The case I want is the simplest one: $G$ **cyclic**. Then the fusion graph is
just a cycle — the affine $A_{n-1}$ diagram — and that cycle is the Cayley graph
$\mathrm{Cay}(\mathbb{Z}/n, \{\pm 1\})$. If you've spent any time around
elliptic curves, a cycle like that will look familiar: it is also the shape of a
**horizontal isogeny cycle**. This post is about why that resemblance is not a
coincidence — and what it teaches us about a stubborn $\pm$ sign.
""")
    st.caption("Aside: those same ADE diagrams also govern the singular fibres of "
               "elliptic surfaces, through the work of Kodaira and Tate. That was "
               "most of my dissertation, and it is a story for another day.")

    st.markdown("### The setup: fusion graphs of a class group")
    st.markdown(r"""
Let $A$ be a finite abelian group. For us it will be the **class group** of an
imaginary quadratic order — the same group that acts, simply transitively, on
each layer of an isogeny volcano. Two facts about a finite *abelian* group make
everything clean: it has exactly $|A|$ irreducible representations, and every one
of them is $1$-dimensional. So the irreps are just the **characters**,

$$\mathrm{Irr}(A) = \mathrm{Hom}(A, \mathbb{C}^\times),$$

and once we fix a set of generators we can label the vertices of any fusion graph
by the elements of $A$ itself. Everything now rides on the single choice of the
distinguished representation $\rho$.

And here is the first surprise: **almost every choice of $\rho$ is useless.** A
$1$-dimensional $\rho = \psi$ just permutes the characters, $\chi \mapsto
\chi\psi$, so its fusion graph is a disjoint union of bare cycles carrying no
information. To get a graph worth reading, we need to choose $\rho$ well.
""")

    st.markdown("### The reveal: real means undirected means ±")
    st.markdown(r"""
Here is the choice, and it is the spine of the whole post. Ask: *when is the
fusion graph undirected?* The answer is exactly:

> the fusion graph is undirected $\iff$ $\rho$ is **self-dual** $\iff$ (over an
> abelian group, where there is no quaternionic type) $\rho$ is a **real**
> representation.

A real representation of $A$ is one of two things. Either a character that already
lands in $\{\pm 1\}$ — a real $1$-dimensional irrep, with $\chi^2 = 1$ — or, for a
genuinely complex character $\chi$, the $2$-dimensional real representation
$\chi \oplus \bar\chi$. Feed that second kind in as $\rho$ and watch what happens:

$$\chi' \otimes (\chi \oplus \bar\chi) = \chi'\chi \;\oplus\; \chi'\bar\chi,$$

so every vertex connects to $\chi' \cdot \chi^{\pm 1}$. The fusion graph is
**exactly the Cayley graph** $\mathrm{Cay}(A, \{\chi^{\pm 1}\})$. This is the
point I care about most: the $\pm$ is not something we lost by "only looking at
the graph." A real $2$-dimensional representation simply *is* the unordered pair
$\{\chi, \bar\chi\}$ — it never recorded which of the two it came from. The
ambiguity is built into the word "real."
""")
    st.info(r"""
**A physicist's translation** (which I trust further than I understand it). A
representation is a *particle* — this is Wigner's classification. Tensoring
representations is *fusing* particles, the physics of adding angular momenta. The
dual representation is the *antiparticle*, and complex conjugation is *charge
conjugation* $C$. So a self-dual — real — representation is a particle that is
**its own antiparticle**, like the photon or a Majorana fermion, and the reason
the graph is undirected is that a $C$-symmetric detector cannot tell a particle
from its antiparticle. (The honest home for "a finite abelian group whose objects
fuse by the group law" is the theory of **abelian anyons**; that is a model of a
class of systems, not one experiment on a bench.)
""")

    st.markdown("### The dictionary")
    st.markdown(r"""
Every column below is the *same* finite abelian group $A$ — the class group —
read three ways. Complex conjugation, inversion, and charge conjugation are one
map wearing three hats.
""")
    st.markdown(r"""
| | Class group $A$ | Representation theory | Particle physics |
|---|---|---|---|
| **objects** | CM points / $j$-invariants | irreps $\mathrm{Hom}(A,\mathbb{C}^\times)$ | particles |
| **conjugation** | inversion $x \mapsto -x$ (Frobenius) | dual $\chi \mapsto \bar\chi$ | charge conjugation $C$ (antiparticle) |
| **self-dual / fixed** | $2$-torsion $A[2]$ (real $j$'s) | real $1$-dim irreps ($\chi^2=1$) | own antiparticle (Majorana / $\pi^0$) |
| **conjugate pair** | $\{x, -x\}$ | real $2$-dim irrep $\chi \oplus \bar\chi$ | particle / antiparticle pair |
| **the $\pm$ ghost** | inversion is invisible to undirected Cayley graphs | self-dual $\Rightarrow$ undirected fusion graph | a $C$-symmetric detector can't tell matter from antimatter |
""")

    st.markdown("### One square")
    st.markdown(r"""
The smallest example worth drawing: discriminant $-39$, whose class group is
$\mathbb{Z}/4$, and whose horizontal $\ell = 2$ isogeny graph is a $4$-cycle.
Below it is drawn as a diamond, coloured in all three languages at once, with the
**real axis horizontal**. The two real $j$-invariants land on that axis (they are
the $2$-torsion, self-dual, their own antiparticles); the complex-conjugate pair
$744 \pm 18197\,i$ sits on the other diagonal. Press the mirror: reflecting across
the real axis swaps the conjugate pair, fixes the two real vertices, and leaves
the graph **completely unchanged**. That invisible reflection — complex
conjugation, Frobenius, and charge conjugation, all the same involution — is the
whole difficulty in one gesture.
""")
    components.html(blog_viz.three_languages_square_html(), height=1000, scrolling=True)

    st.markdown("### A bigger square: when one generator isn't enough")
    st.markdown(r"""
On the $4$-cycle a single character generated everything, so "spanning set" was an
overstatement — one edge type sufficed. The moment $A$ stops being cyclic, the
word *set* earns its keep. Take discriminant $-231$, with class group
$\mathbb{Z}/2 \times \mathbb{Z}/6$: twelve curves, and **no single prime's
isogenies reach them all**. You need a set of generators, and the picture is a
**hexagonal prism**.
""")
    components.html(blog_viz.three_languages_prism_html(), height=900, scrolling=True)
    st.markdown(r"""
The $\ell = 2$ isogenies give a class of **order $6$** — the two hexagons, a real
$2$-dimensional representation $\chi \oplus \bar\chi$. The $\ell = 3$ isogenies
give a class of **order $2$** — the rungs, a real $1$-dimensional representation.
Both flavours of real distinguished representation, side by side in a single
spanning set, and neither one reaches the whole group alone. The mirror is now a
**plane**: it swaps matter with antimatter (left with right), holds the four
self-dual vertices still, and — once again — leaves every edge exactly where it
was. Choosing which side to call matter is a real choice the graph cannot make
for you.
""")

    st.markdown("### The landing")
    st.markdown(r"""
That last choice is the whole game. A spanning set of these undirected fusion
graphs pins down $A$ only **up to inversion** $x \mapsto -x$ — the one
automorphism that is invisible to every symmetric Cayley graph at once. A
**rigid** spanning set is precisely one carrying enough extra structure to remove
that final ambiguity. And that is exactly the computation this site runs to match
CM lattices with elliptic curves over $\mathbb{F}_p$ at the floor of every
isogeny volcano — now restated with not an elliptic curve in sight.

The honest caveat, stated plainly: this **repackages** the rigid-spanning-set
problem as a curve-free toy; it does not reprove it. But I find it genuinely
clarifying that the $\pm$ you fight in an isogeny graph, the $\pm$ in "which
character did this real representation come from," and the matter/antimatter $\pm$
are the *same ghost* — inversion, wearing three hats.
""")


# ──────────────────────────────────────────────────────────────────────────────
POSTS = [
    {
        "slug": "sum-of-two-squares",
        "title": "Minkowski and Frobenius",
        "date": "July 2026",
        "teaser": "Fermat says a prime p ≡ 1 (mod 4) is a sum of two squares. "
                  "Here's how to *find* the squares — as the shortest vector of a "
                  "lattice — with a descent you can watch, and what it says about "
                  "the Bitcoin curve.",
        "render": render_sum_of_squares,
    },
    {
        "slug": "three-languages",
        "title": "One square, three languages",
        "date": "July 2026",
        "teaser": "The same class group, read as geometry, as representation "
                  "theory, and as particle physics — and a ± ambiguity that turns "
                  "out to be the same ghost in all three.",
        "render": render_three_languages,
        "draft": True,
    },
]


def render_index():
    st.markdown("# Blog")
    st.markdown(
        "Short, self-contained stories from the making of this site — each one a "
        "familiar computation seen from an unfamiliar angle."
    )
    st.divider()
    for post in POSTS:
        tag = " &nbsp;·&nbsp; *draft*" if post.get("draft") else ""
        st.markdown(
            f"### [{post['title']}](/Blog?post={post['slug']})",
        )
        st.markdown(
            f"<span style='opacity:.7'>{post['date']}{tag}</span>",
            unsafe_allow_html=True,
        )
        st.markdown(post["teaser"])
        st.markdown(f"[Read &rarr;](/Blog?post={post['slug']})", unsafe_allow_html=True)
        st.divider()


# ── routing ──────────────────────────────────────────────────────────────────
slug = st.query_params.get("post")
post = next((p for p in POSTS if p["slug"] == slug), None)

if post is None:
    render_index()
else:
    st.markdown("[&larr; All posts](/Blog)", unsafe_allow_html=True)
    post["render"]()
    st.divider()
    st.markdown("[&larr; All posts](/Blog)", unsafe_allow_html=True)
