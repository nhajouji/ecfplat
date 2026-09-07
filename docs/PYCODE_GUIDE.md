# The pycode guide

What every module in `pycode/` is for, what its key entry points are, and who
uses it — plus the ledger of everything removed in the 2026-08 cleanup and
where each removed thing's successor lives.  Companion documents:
`notebooks/userguide.ipynb` (the same API, runnable), `docs/X0L_LEGACY.md`
(design notes for the removed X₀(ℓ) toolkit), `docs/EXPLORER_DESIGN.md`
(the static-site plan).

Layer map, bottom to top — imports only ever point downward:

```
nt ─ alg_classes ─ misctools                      elementary number theory / algebra
   qfs ─ identities ─ clgp                        binary quadratic forms & class groups
      modularpolynomials ─ ecfp ─ velu            curves over F_p, modular polynomials
         graph_tools ─ ecqf_bij                   the bijection engine
            hilbert_crt ─ modpoly_crt             the CRT polynomial pipelines
            rigid_cache ─ ss_bij_cache ─ ext_store ─ bij_factory ─ phi_factory   caches & factories
               ecqf ─ ecqf_tools                  the website-facing classes
                  *_viz ─ palette ─ shortlinks ─ ss_graph(_page) ─ pell_duality  applets & research
```

## Core math

| module | purpose | key API | used by |
|---|---|---|---|
| `nt.py` | Elementary number theory: gcd/CRT, trial-division factoring, quadratic characters, primes, square roots mod p, Frobenius extension degrees, sums-of-squares descent (the blog-post-1 material). | `primefact`, `discfac`, `primeQ`, `primesBetween`, `quad_rec` (Legendre symbol), `crt_list`/`crt_pair`, `sqrt_mod_prime`, `frob_ext_degrees`, `hall_multiplier` (live in the MW chain), `sos`/`esos`/`x2_3y2` | everything |
| `alg_classes.py` | The algebra layer, three strata: abstract groups/rings (`AbGrp`, `Ring`, `Field`, matrices), concrete fields (`GF_p`, `GF_pn`, `QuadExt`), and polynomial rings (`PolyRing`/`Poly`, the ring-pattern generation; plus the *legacy* `Polynomial`/`PolyFp`, still live only via `modularpolynomials` — scheduled for removal after that reroute). Kept intact per Nadir: parts serve side projects. | `Mat_n_Z`, `MatrixElement`, `GF_p`, `GF_pn_auto`, `ZnProduct`, `AbGrElt`, `abgrp_basis`, `poly_ring`, `poly_crt`, `ZZ`/`CC` | qfs, clgp, velu, hilbert_crt, ecqf_tools |
| `misctools.py` | Four small tuple/dict helpers. | `sort_tuple`, `merge_dicts` (note: mutates its first argument), `compdiv_dics` | graph_tools, ecqf_bij |
| `qfs.py` | Binary quadratic forms and the SL₂(ℤ) action: reduction, class-group enumeration, ℓ-isogeny neighbours (horizontal/ascending), Γ₀(ℓ) cosets and the Fricke involution. | `qf_mod_gamma`, `_qf_reduce`, `get_qfs_strict`/`get_qfs_all`, `qfs_ordered_by_cond`, `qf_isogs_hor`, `qf_parents`, `qf_isog_cycle`, `class_group_id`/`class_group_inv`, `qf_x0_endos`, `gamma_0_orb`, `fricke_inv`, `qf_ev` (kept for prime-representation searches), `qf_2_mat`/`mat_2_qf` (mult-by-aτ matrix; NB `qf_to_mat` is the *Gram* matrix — different map) | ecqf, ecqf_bij, clgp, hilbert_crt, modpoly_crt, rigid_cache, viz |
| `identities.py` | Class-number identities and counts: twisted Euler φ, class-number sums over conductors, the supersingular trace/Eichler counts, discriminant closures. Canonical home of `twisted_phi` (the `nt` copy was deleted). | `twisted_phi(_sum)`, `clgr_size_gen`, `clgr_sum`, `supsingtrace`, `disc_closure`, `x0l_fp_card` | pell_duality, modpoly_crt |
| `clgp.py` | Class groups as honest groups: Gauss composition (`_compose_raw`, Cohen 5.4.7), `ClassGroup`/`ClGrpElt` on the `AbGrp` machinery, generator searches scored by prime, the fast `rigid_lset` path. Part of the package API. | `class_group(d)`, `ClassGroup.order/basis/reduce`, `class_from_prime`, `rigid_lset` | scripts/validate_clgp.py; integration with ecqf_bij planned |

## Curves over F_p and modular polynomials

| module | purpose | key API | used by |
|---|---|---|---|
| `ecfp.py` | Curves over F_p by F_p methods only: j ↔ (f,g) model conversions, the signature invariant (curve vs quadratic twist), vectorized trace-of-Frobenius tables, Atkin isogeny codomains. **Planned home of an `EllipticCurveFp` class and, later, the Vélu engine.** | `j_to_fg`/`fg_to_j`, `signature`, `js_to_fg`, `trace_frob`, `_trace_table` (cached, GPU-pluggable via `TRACE_TABLE_IMPL`), `trfr_to_js`, `trfr_to_models`, `fp_isog_codomains` | ecqf_bij, velu, hilbert_crt, ss_bij_cache |
| `modularpolynomials.py` | The modular-polynomial store: Atkin Φₗ* for the 15 Atkin primes, classical Φₗ from the CRT pipeline, Hilbert polynomial dictionaries, the q-expansion coefficients. Loads its JSON stores at import (lazy-loading is a planned packaging fix). | `atkin_polys_dict`, `atk_at_j`, `eval_atk`, `classical_modpoly`, `register_modpoly`, `hilb_polys_dict`, `heeg_js`, `modular_prime_pool`, `small_bij_check` | ecfp, ecqf_bij, rigid_cache, hilbert_crt, phi_factory |
| `velu.py` | The Vélu isogeny engine over F_p and its extensions: codomain computation from a kernel point, eigenline walks, 2-isogeny signatures, the supersingular coset walk. Works with the `alg_classes` field towers. | `velu_l_isog_codomain_fast`, `velu_l_isog_codomain_twist`, `velu_walk_cycle`, `velu_nbr_data_ord`, `velu_nbr_data_ss_walk`, `two_isogeny_sigs`, `embed_fp`, `random_point`, `ec_mul` | ecqf_bij, ss_bij_cache |
| `graph_tools.py` | Abstract graph utilities for the bijection: neighbour-data tree searches (the `_zn` chain), cycle walks, the cube/backtracking assembly for noncyclic class groups, `compute_bijection_zn`. | `nbrdata_tree_search_zn` (implemented via `_XC1` → `_XCyc` → base — the whole chain is live), `cycle_from_neighbor_data`, `oriented_cycle_from_nbdata`, `nbrdata_to_isomat`, `compute_bijection_zn`, `verify_bijection` (future pytest seed) | ecqf, ecqf_bij, rigid_cache |
| `ecqf_bij.py` | **The core algorithm**: the bijection between curves over F_p (j-side / signature-side) and quadratic forms (lattice side), via rigid ℓ-set searches, descriptor dispatch (pin/powkey/sib/lift), and the Vélu or Atkin edges. | `ecqf_full_bijection_ord`, `ecqf_full_bijection_ss`, `disc_rigid_lset_search`, `qf_isog_data`, `ecfp_nbr_data_ord(_X1)`, `ss_signatures`, `get_ancestor_data_ord`, `ssprimes` | rigid_cache, bij_factory, ss_bij_cache, hilbert_crt, ldata (via rigid_cache) |

## Pipelines, caches, stores

| module | purpose | key API / CLI | status |
|---|---|---|---|
| `hilbert_crt.py` | Hilbert class polynomials by CRT over ordinary primes: harvest CM j-invariants, identify endomorphism discs (ancestor data or the elimination trick), assemble H_d exactly. Validated 81/81 against the known table. | `hilbert_via_crt`, `hilbert_poly_search(d, N)`, `hilbert_library`, `find_aps`, `certification_report` | library + notebook driver |
| `modpoly_crt.py` | Classical Φ_p by nonstandard CRT: diagonal from Hilbert polys, special-value rows from the 13 class-number-1 j's, interpolation, Kronecker-congruence check. | `phi_p_via_crt`, `phi_monomials`, `phi_diagonal`, `special_value_rows`, `solve_phi_from_pairs` | library |
| `phi_factory.py` | Batch driver for `modpoly_crt` over p ≤ 71 with checkpointing (`data/phi_factory_state.json`). Φ₂…Φ₆₇ shipped; Φ₇₁ pending a laptop run. | `run_factory`, `factory_targets` | pipeline script |
| `bij_factory.py` | Batch driver extending the bijection store beyond p = 1024 (output: the 83 MB ext JSON), plus the Hilbert harvest factory. | `compute_and_save`, `load_ext_bijections`, `hilbert_factory` | pipeline script |
| `rigid_cache.py` | Per-discriminant cache of rigid ℓ-set searches + bijections (`rigid_lset_cache.json`, 32 MB) with the runtime accessor the factories use; also regenerates the lightweight `qf_ldata.json` (absorbed from the deleted `ldata_cache.py`). | `ecqf_ord_bij_cached`, `get_disc_entry`, `populate`; CLI `python rigid_cache.py [--ldata] --min -4096` | library + CLI |
| `ss_bij_cache.py` | Regenerates and validates the supersingular bijection table (`ecqf_ss_pcbij_velu_4_1024.json`, primes to 8192). `validate_entry` is the 4-leg validation battery (sets, genuine-SS, Φ-edges, root convention) — future pytest material. | `populate`, `validate_entry`, `bijection_entry`; CLI `python ss_bij_cache.py --max 2048` | pipeline + validation |
| `ext_store.py` | Lazy access to the 117k-class extended ordinary store. Resolution: `$ECFPLAT_EXT_DB` → default sqlite → whole-JSON fallback (laptop only; the deployed site needs the sqlite for memory). | `get_pair(a, p)`, `has_pair`, `stats` | ecqf (website hot path) |

## Website-facing classes

| module | purpose | key API | used by |
|---|---|---|---|
| `ecqf.py` | The assembly point: `QFIsogenyClass(d)` (characteristic-0 view) and `ECQFIsogenyClass(a, p)` (Frobenius view, with j's/models/MW when the stores cover (a,p)); the graph descriptor the canvas volcano renders; the disc ↔ (a,p) bridge. | `QFIsogenyClass`, `ECQFIsogenyClass` (`.ecqf_df()`, `.qf_to_mw_gens_dict(k)`, `.qf_to_mwgroups_alltups(k)`, `.isog_cycle`), `class_graph_descriptor`, `disc_to_aps`, `coverage`, `P_MAX`/`GUARDS` | Explorer page, netlify-proto |
| `ecqf_tools.py` | The per-form machinery under `ecqf`: the eagerly-loaded 1K stores, model selection over F_p, and the **Mordell–Weil chain** — Frobenius matrix on the lattice basis (1, τ), kernel generators via Hall multipliers, point enumeration. | `ecqf_ord_1K_pc`/`ecqf_ss_1K_pc`, `ap_in_pc_data`, `ecfp_js_to_model`, `qf_ap_FrMat`, `qf_to_ERGM_1T` (mult-by-aτ; ≡ `qfs.qf_2_mat` — consolidation planned), `frob_to_mw_gens`, `qf_mat_ker_gens`, `pts_from_gendic`, `ec_look_up`, `abc_to_tau(_str)` | ecqf, Explorer page |

## Visualization & site

| module | purpose | used by |
|---|---|---|
| `basics_viz.py` | 17 canvas applets for Background ch. 1–2 (group law, F_p pictures, CM tori, Vélu, Frobenius flow). | Background |
| `modular_viz.py` | 8 applets for Background ch. 3 (moduli, X₀(ℓ), modular polynomials); owns the shared `_HEAD` CSS. | Background |
| `explorer_viz.py` | The Explorer's 5 linked canvas views (volcano graph, Hasse picker, FD points, curve torus, SS graph). Pure data-in/HTML-out; also consumed by `netlify-proto/build.py` (which splits on its exact output strings — don't reformat casually). | Explorer, static proto |
| `slide_viz.py`, `hero_viz.py`, `hopf_viz.py`, `blog_viz.py` | Talk-derived applets, the homepage hero, the Three.js Hopf torus, blog-post applet loader. | Homepage/Background/Blog |
| `ss_graph.py` / `ss_graph_page.py` | The F_p supersingular ℓ-isogeny graph (built on `pell_duality`'s F_q toolkit) and its Streamlit page (blog post 3). `ss_graph_page` is app glue that will move to `pages/` at packaging time. | Blog |
| `plotly_tools.py` | **Notebook-only**: interactive plotly figures (FD pictures, volcano layouts) for pure-Python environments. Not a deploy dependency (plotly isn't in requirements.txt) — install plotly locally to use it. | notebooks |
| `palette.py`, `shortlinks.py`, `trace_gpu.py` | Shared colors; `/Gallery?g=` short links; optional torch backend that plugs into `ecfp.TRACE_TABLE_IMPL`. | Explorer/Gallery/notebooks |
| `pell_duality.py` | **Active research module** (kept whole, deliberately): the p–ℓ duality experiments — F_q arithmetic, division polynomials, kernels/Vélu at supersingular vertices, the 2p+2 identity, genus formulas, pairing catalogs. `ss_graph` imports its F_q/Vélu toolkit. Stage results cache to `experiments/pell_duality_results.json` (gitignored). | ss_graph, ongoing experiments |

## Data files (`pycode/data/`)

| file | producer | consumer |
|---|---|---|
| `atkinpolys.json`, `hilbpolys.json`, `jcoefs.json`, `jq_coeffs.json` | sourced/precomputed | `modularpolynomials` (import-time) |
| `classical_modpolys.json` (13 MB) | `phi_factory` | `modularpolynomials.classical_modpoly` |
| `hilbpolys_crt.json`, `hilbert_roots_ext.json` | `hilbert_crt` / `bij_factory.hilbert_factory` | `hilbert_crt.hilbert_library` |
| `ecqf_ord_pcbij_4_1024.json` (1.6 MB) | original pipeline | `ecqf_tools.ecqf_ord_1K_pc` (eager) |
| `ecqf_ss_pcbij_velu_4_1024.json` (1.6 MB) | `ss_bij_cache` | `ecqf_tools.ecqf_ss_1K_pc` (eager) |
| `ecqf_ord_pcbij_ext.json` (83 MB) | `bij_factory` | factory-side source; `ext_store` JSON fallback |
| `ecqf_ord_ext.sqlite` (39 MB) | `scripts/build_ext_sqlite.py` | `ext_store` (the deployed site's lazy store) |
| `rigid_lset_cache.json` (32 MB), `qf_ldata.json` | `rigid_cache` (`--ldata` for the second) | factories; `ecqf_bij` (import-time) |
| `phi_factory_state.json` (8.5 MB) | `phi_factory` checkpoints | only the Φ₇₁ entry still matters |
| `hopf_disc8_k2.json` | Steve-side export | `hopf_viz` |

## Deletion ledger — 2026-08 cleanup (branch `cleanup`, tag `pre-cleanup`)

Everything below was verified zero-caller (fresh grep per deletion) before
removal; `git show pre-cleanup:pycode/<file>` recovers any of it.  Old name →
where that job is done now:

| removed | replacement / note |
|---|---|
| `graphic_tools.py` (matplotlib) | `plotly_tools.py` kept as the notebook option; site applets are canvas |
| `ldata_cache.py` | `python rigid_cache.py --ldata` (`rigid_cache.populate_ldata`) |
| `ecfp.sig_to_qfs_dict`, `get_j_to_qfs_dict` | `ecqf_bij.ecqf_full_bijection_ord` / `_ss` (these were also broken — NameError) |
| `ecfp.frobmat` | `ecqf_tools.qf_ap_FrMat(qf, (a, p))` |
| `ecfp.kernel_gen_cyc`, `divide_cyclic_gen`, `mw_gens` | `ecqf_tools.qf_mat_ker_cyc` / `divide_cyclic_gen` / `frob_to_mw_gens` (the working copies; `pts_from_gendic` fixed in the same pass) |
| `ecfp.get_endo_disc_cands`, `j_to_disc`, `endo_db_check`, `supp*` | endomorphism discs now come from ancestor data: `hilbert_crt.class_endo_discs` |
| `ecfp.trfr_to_leaves`/`trfr_to_aec_data`/`check_leaf`, `all_ssl_cycles_from_jp`, `j0_to_j1s_in_sslcycs`, `cubic_qrs`, `trace_and_2tor`, `get_precomputed_ssdict` | leaf/cycle logic lives in `ecqf_bij`'s neighbour-data builders |
| `nt.twisted_phi(_sum)` | `identities.twisted_phi(_sum)` (were byte-identical copies) |
| `nt.class_no_formula`, `jacobi_symbol` | `len(qfs.get_qfs_strict(d))` or `clgp.ClassGroup.order` (~1000× faster) |
| `nt.int_sqrt`, `get_rou_mod`, `mod_sfd`, `pf_to_int`, `quad_gcd`, `ap_to_lm` | `math.isqrt`; the rest had no successors needed |
| `qfs.qf_to_fun_dom`, `minv`, `find_rrep_g0`, `qf_to_gamma_0_fd`, `qf_mod_gamma_0`, `x0_endos_all`, `iso_taus_x0_l`, `isos_x0_l_all` | **see `docs/X0L_LEGACY.md`** — kept: `qf_x0_endos`, `gamma_0_orb`, `fricke_inv` |
| `qfs.qf_to_dc`, `prod_tup`, `qf_isogs_asc`, `qfs_isogs_int`, `cycs_from_ancestors` | `discfac(qf_disc(qf))`, `math.prod`, `qf_parents`; last two had a bug / no callers |
| `ecqf.adjacency_matrix` | `graph_tools.nbrdata_to_isomat(cls.get_neighbor_data_all(l), cls.qfs_ordered)` |
| `ecqf.isog_cycle_partition`, `ecqf_mw_df`, `qf_to_mwgr_arr_single`, `.l_dict`/`.ord_dict`/`.qfs_leaves` | `isog_cycle` per form; `ecqf_df` + `qf_to_mw_gens_dict`; per-form smallest split ℓ: `ecqf_tools.qf_l_order(qf)` on demand |
| `ecqf_tools.export_points`, `ap_FrbMats_1T`, `qf_cond` | `qf_ap_FrMat` in a loop; `discfac(qf_disc(qf))[1]` |
| `ecqf_bij.zn_ecqf_bij`, `qf_reps_pm`, `_desc_order` | `ecqf_full_bijection_ord` |
| `velu.point_order`, `velu_nbr_data_ss` | `velu_nbr_data_ss_walk` |
| `modularpolynomials.modpoly_nbrs`, `modpoly_primes`, `modpoly_from_terms`, `rat_eval_mod`, `count_roots_fp_bf`, `atk_at_j_fpfac`, `_spow` | `modpoly_roots_among`, `modular_prime_pool` |
| `hilbert_crt.crt_prime_candidates`, `extend_disc`, `ATKIN_SET` | `find_aps`, `hilbert_poly_search` |
| `misctools.invert_dict` | `{v: k for k, v in d.items()}` |
| `ss_bij_cache.write_list_form` + `data/ssfp_pc_bij_velu.json` | the dict-form file is the single SS store (list form fed only dead code) |
| `data/ecqf_ord_pcbij_4to256.json`, `testing_discriminants.json`, `ssfp_isog_cycles_10b.json` | superseded / orphaned |

### If an old notebook breaks

The June-2026 notebooks in `experiments/` (`neighbordata`, `QF isographs`,
`testing2`) reference `zn_ecqf_bij`, `qf_reps_pm`, `adjacency_matrix`,
`invert_dict`, `frobmat`, `mw_gens` — map them through the table above, or run
them against the `pre-cleanup` tag.  Everything the **pell-duality, SOS, and
factory** work uses is untouched.  `notebooks/userguide.ipynb` has been
rewritten against the current API.
