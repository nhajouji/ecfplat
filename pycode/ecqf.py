"""ecqf — the assembly point for the ecfplat website.

This module owns the two main classes (QFIsogenyClass, ECQFIsogenyClass) and
the Explorer-facing helpers built on top of them. Pages should import from
here; the underlying machinery (quadratic forms, Frobenius matrices, the
precomputed j <-> qf stores, lookups) lives in ecqf_tools and friends.

Entry points, matching the Explorer's two doors:
  * discriminant-first:  QFIsogenyClass(d)      — characteristic-0 view
  * (a, p)-first:        ECQFIsogenyClass(a, p) — Frobenius view over F_p
  * the bridge:          disc_to_aps(d)         — pairs with a^2 - 4p = d m^2
"""

import numpy as np
import pandas as pd

from nt import discfac, quad_rec, primeQ
from qfs import (qfs_ordered_by_cond, qf_disc, qf_isogs_hor, qf_parents,
                 class_group_id)
from graph_tools import nbrdata_to_isomat, cycle_from_neighbor_data
from ecqf_tools import (
    # lookups and stores
    ec_look_up, ap_in_pc_data, get_aps_pc, get_ssps_pc,
    ecqf_ord_1K_pc, ecqf_ss_1K_pc,
    # per-form machinery used by the classes
    qf_l_order, qf_ap_FrMat, ecfp_js_to_model,
    frob_to_mw_gens, pts_from_gendic, mw_arr_from_gens,
    abc_to_tau, abc_to_tau_str, ec_eq_str_base, ec_eq_str,
)


# ══ The main classes ══════════════════════════════════════════════════════════

class QFIsogenyClass:
    """The isogeny class of a negative discriminant d — characteristic 0.

    No prime, no distinguished endomorphism: just the forms of every
    conductor dividing cond(d), their l-neighbour structure (horizontal
    cycles + vertical ascents = the volcano), and the class-group data."""

    def __init__(self, d: int):
        if d >= 0 or (d % 4 > 1):
            raise ValueError(f'{d} must be a negative discriminant')
        self.qfs_all = qfs_ordered_by_cond(d)
        self.qfs_ordered = tuple(self.qfs_all)
        self.neighbor_data_horz = {}
        self.vert_isog_data = {}
        self.neighbor_data = {}
        self.disc = d
        d0, c = discfac(d)
        self.field_disc = d0
        self.cond = c
        self.qfs_leaves = [qf for qf in self.qfs_all if qf_disc(qf) == d]
        # discriminant / conductor of each form's endomorphism ring
        self.endo_disc_dict = {qf: qf_disc(qf) for qf in self.qfs_all}
        self.endo_cond_dict = {qf: discfac(qf_disc(qf))[1] for qf in self.qfs_all}
        self.l_dict = {}
        self.ord_dict = {}
        for qf in self.qfs_all:
            l, n = qf_l_order(qf)
            self.l_dict[qf] = l
            self.ord_dict[qf] = n

    def get_isog_neighbors_horz(self, l: int):
        if l in self.neighbor_data_horz:
            return self.neighbor_data_horz[l]
        qfs = self.qfs_ordered
        data = {qf: qf_isogs_hor(qf, l) for qf in qfs}
        self.neighbor_data_horz[l] = data
        return data

    def get_isog_neighbors_asc(self, l: int):
        c = self.cond
        if c % l != 0:
            return {}
        if l in self.vert_isog_data:
            return self.vert_isog_data[l]
        asc_data = {}
        for qf in self.qfs_ordered:
            d_qf = qf_disc(qf)
            c_qf = discfac(d_qf)[1]
            if c_qf % l == 0:
                asc_data[qf] = qf_parents(qf, l)[0]
        self.vert_isog_data[l] = asc_data
        return asc_data

    def get_neighbor_data_all(self, l):
        if l in self.neighbor_data:
            return self.neighbor_data[l]
        # Copy the cached horizontal lists before appending the vertical
        # neighbours, otherwise the in-place append would pollute the cached
        # horizontal-neighbour data (which isog_cycle and the graph layout
        # rely on being purely horizontal).
        horiz = self.get_isog_neighbors_horz(l)
        neighbors_data_l = {qf: list(horiz[qf]) for qf in horiz}
        if self.cond % l == 0:
            vert_data = self.get_isog_neighbors_asc(l)
            for qf0 in vert_data:
                qf1 = vert_data[qf0]
                neighbors_data_l[qf0].append(qf1)
                neighbors_data_l[qf1].append(qf0)
        self.neighbor_data[l] = neighbors_data_l
        return neighbors_data_l

    def adjacency_matrix(self, l):
        data = self.get_neighbor_data_all(l)
        return nbrdata_to_isomat(nbrdata=data, verts_ordered=self.qfs_ordered)

    def isog_cycle(self, qf0: tuple, l: int):
        if qf0 not in self.qfs_all:
            raise ValueError(f'{qf0} is not in isogeny class')
        lnbr_data = self.get_isog_neighbors_horz(l)
        return cycle_from_neighbor_data(qf0, lnbr_data)

    def isog_cycle_partition(self, l):
        cond_dict = self.endo_cond_dict
        qfs_by_cond = {c: [] for c in cond_dict.values()}
        for qf, c in cond_dict.items():
            qfs_by_cond[c].append(qf)
        cycles_by_cond = {c: [] for c in qfs_by_cond if c % l != 0}
        for c in cycles_by_cond:
            qfs_c = qfs_by_cond[c]
            while len(qfs_c) > 0:
                cyc_new = self.isog_cycle(qfs_c[0], l)
                assert len(cyc_new) > 0
                cycles_by_cond[c].append(cyc_new)
                qfs_c_new = [qf for qf in qfs_c if qf not in cyc_new]
                assert len(qfs_c_new) + len(cyc_new) == len(qfs_c)
                qfs_c = qfs_c_new
        return qfs_by_cond


class ECQFIsogenyClass(QFIsogenyClass):
    """The isogeny class of trace a over F_p — the Frobenius view.

    Extends the characteristic-0 class of d = a^2 - 4p with the Frobenius
    matrices, and (when the precomputed stores cover (a, p)) the actual
    j-invariants, Weierstrass models, and Mordell-Weil data."""

    def __init__(self, a: int, p: int):
        super().__init__(a * a - 4 * p)
        self.a_sign = 1
        if a < 0:
            self.a_sign = -1
        self.disc = a * a - 4 * p
        nonsquare = p - 1
        self.ap = (a, p)
        self.char = p
        self.trace = a
        while nonsquare > 1 and quad_rec(nonsquare, p) >= 0:
            nonsquare -= 1
        self.nonsquare = nonsquare
        self.qf_to_frob_mats = {qf: qf_ap_FrMat(qf, (a, p), s=self.a_sign)
                                for qf in self.qfs_all}
        self.js_to_qf = None
        self.js = None
        self.jsigs = None
        self.js_to_models = {}
        if a != 0:
            app = (abs(a), p)
            if app in get_aps_pc():
                self.js_to_qf = ecqf_ord_1K_pc[app]
                self.jsigs = [j for j in self.js_to_qf]
                self.js = self.jsigs
                self.js_to_models = {js: ecfp_js_to_model(js, (a, p), self.nonsquare)
                                     for js in self.js_to_qf}
        if a == 0:
            if p in get_ssps_pc():
                self.js_to_qf = ecqf_ss_1K_pc[p]
                self.jsigs = [j for j in self.js_to_qf]
                self.js = list(set(js[0] for js in self.js_to_qf))
                self.js_to_models = {js: ecfp_js_to_model(js, (a, p), self.nonsquare)
                                     for js in self.js_to_qf}

    def qf_to_mw_gens_dict(self, k: int):
        qf_to_frms = self.qf_to_frob_mats
        return {qf: frob_to_mw_gens(qf_to_frms[qf], k) for qf in qf_to_frms}

    def qf_to_mwgroups_alltups(self, k: int = 1):
        qf_gens = self.qf_to_mw_gens_dict(k)
        return {qf: pts_from_gendic(qf_gens[qf]) for qf in qf_gens}

    def qf_to_mwgr_arr_single(self, k: int = 1, qf: tuple = None):
        if qf is None:
            qf = (self.qfs_all)[0]
        if qf not in self.qf_to_frob_mats:
            raise ValueError('Form not in dictionary')
        frm = self.qf_to_frob_mats[qf]
        mwgens = frob_to_mw_gens(frm, k)
        if len(mwgens) == 0:
            return [np.array([0, 0])]
        return mw_arr_from_gens(qf, mwgens)

    def ecqf_df(self):
        if self.js_to_qf is None:
            raise ValueError('No data available')
        jss = self.jsigs
        if type(jss[0]) == tuple:
            jlist = [js[0] for js in jss]
        else:
            jlist = jss
        fgs = [(self.js_to_models)[js] for js in jss]
        qfs = [(self.js_to_qf)[js] for js in jss]
        frobmats = [self.qf_to_frob_mats[qf].vec for qf in qfs]
        qfds = [qf_disc(qf) for qf in qfs]
        qf_cs = [discfac(d)[1] for d in qfds]
        qf_ccs = [self.cond // c for c in qf_cs]
        tau_xys_arr = [abc_to_tau(qf) for qf in qfs]
        tau_xys = [tuple([np.round(x, 3) for x in xy]) for xy in tau_xys_arr]
        tau_strs = [abc_to_tau_str(qf) for qf in qfs]
        return pd.DataFrame({'ec_invs': jss, 'j_inv': jlist, 'EC_coefs': fgs,
                             'qf_coefs': qfs, 'endo_disc': qfds,
                             'endo_cond': qf_cs, 'endo_cocond': qf_ccs,
                             'frobmat': frobmats, 'tau_s': tau_strs,
                             'tau_xy': tau_xys})

    def ecqf_mw_df(self, k: int):
        if self.js_to_qf is None:
            raise ValueError('No data available')
        jss = self.jsigs
        if type(jss[0]) == tuple:
            jlist = [js[0] for js in jss]
        else:
            jlist = jss
        fgs = [(self.js_to_models)[js] for js in jss]
        qfs = [(self.js_to_qf)[js] for js in jss]
        frmats = [self.qf_to_frob_mats[qf] for qf in qfs]
        qfds = [qf_disc(qf) for qf in qfs]
        qf_cs = [discfac(d)[1] for d in qfds]
        qf_ccs = [self.cond // c for c in qf_cs]
        tau_xys = [abc_to_tau(qf) for qf in qfs]
        tau_strs = [abc_to_tau_str(qf) for qf in qfs]
        frmat_tups = [frm.vec for frm in frmats]
        mwgsets = []
        mwntups = []
        for frm in frmats:
            gendic = frob_to_mw_gens(frm, k)
            genvs = [g for g in gendic]
            genvs.sort(key=lambda g: gendic[g], reverse=True)
            mwgsets.append(genvs)
            mwntups.append(tuple([gendic[g] for g in genvs]))
        return pd.DataFrame({'ec_invs': jss, 'j_inv': jlist, 'EC_coefs': fgs,
                             'qf_coefs': qfs, 'endo_disc': qfds,
                             'endo_cond': qf_cs, 'endo_cocond': qf_ccs,
                             'frobmat': frmat_tups, 'tau_s': tau_strs,
                             'tau_xys': tau_xys, 'MW_gens': mwgsets,
                             'MW_iso_type': mwntups})


# ══ Explorer assembly ═════════════════════════════════════════════════════════

# Hard ceilings for the website. Data coverage ends at 8192; the picture
# guards keep heavy renders from being attempted far beyond their design
# size — above a guard the Explorer degrades to a lighter view, it never
# crashes.
P_MAX = 8192
GUARDS = {
    'volcano_nodes': 600,    # full volcano drawing above this -> summary view
    'fd_points': 2000,       # fundamental-domain scatter cap
    'curve_grid_p': 257,     # F_p point-grid pictures only for p <= this
}


def disc_to_aps(d: int, p_max: int = P_MAX, include_imprimitive: bool = True):
    """All pairs (a, p) with a^2 - 4p = d * m^2, p prime <= p_max, a > 0.

    The bridge from the discriminant-first view to the (a, p) view. m = 1
    gives the classes whose Frobenius discriminant is d itself; m > 1 gives
    the pairs where d sits above the curve's discriminant (only when
    include_imprimitive). Returns dicts {a, p, m, in_store} sorted by p."""
    if d >= 0 or (d % 4 > 1):
        raise ValueError(f'{d} must be a negative discriminant')
    out = []
    m = 1
    while d * m * m + 4 * p_max >= 4:          # smallest a^2 = d m^2 + 4p
        a = 1 if (d * m) % 2 else 2            # parity: a^2 ≡ d m^2 (mod 4)
        while a * a - d * m * m <= 4 * p_max:
            q4 = a * a - d * m * m
            if q4 % 4 == 0:
                p = q4 // 4
                if p > 2 and primeQ(p) and a % p != 0:   # ordinary pairs
                    out.append({'a': a, 'p': p, 'm': m,
                                'in_store': ap_in_pc_data((a, p))})
            a += 2
        if not include_imprimitive:
            break
        m += 1
    out.sort(key=lambda e: (e['p'], e['m']))
    return out


def coverage():
    """What the precomputed stores actually contain (for bounds + gap notes)."""
    aps = get_aps_pc()
    ssps = get_ssps_pc()
    return {
        'ordinary_pairs': len(aps),
        'ordinary_p_max': max(p for (_, p) in aps),
        'ss_primes': len(ssps),
        'ss_p_max': max(ssps),
    }


def class_graph_descriptor(cls: QFIsogenyClass, l: int) -> dict:
    """JSON-able descriptor of the l-isogeny graph of an isogeny class.

    The contract consumed by the canvas volcano component: nodes carry the
    form, its endo-ring data and tau; edges are typed 'h' (horizontal =
    crater/cycle) or 'v' (vertical = tree). For an ECQFIsogenyClass with
    store coverage, nodes also carry their j-invariants and models."""
    qfs = list(cls.qfs_ordered)
    idx = {qf: i for i, qf in enumerate(qfs)}
    jmap = {}
    if isinstance(cls, ECQFIsogenyClass) and cls.js_to_qf:
        for js, qf in cls.js_to_qf.items():
            j = js[0] if isinstance(js, tuple) else js
            jmap.setdefault(qf, []).append(
                {'j': j, 'model': list(cls.js_to_models.get(js, ()))})
    nodes = []
    for qf in qfs:
        tau = abc_to_tau(qf)
        nodes.append({
            'qf': list(qf),
            'endo_disc': cls.endo_disc_dict[qf],
            'endo_cond': cls.endo_cond_dict[qf],
            'tau': [float(np.round(tau[0], 5)), float(np.round(tau[1], 5))],
            'curves': jmap.get(qf, []),
        })
    edges = []
    seen = set()
    horiz = cls.get_isog_neighbors_horz(l)
    for qf, nbrs in horiz.items():
        for nb in nbrs:
            key = (min(idx[qf], idx[nb]), max(idx[qf], idx[nb]))
            if key not in seen and qf != nb:
                seen.add(key)
                edges.append({'s': key[0], 't': key[1], 'kind': 'h'})
    for qf0, qf1 in cls.get_isog_neighbors_asc(l).items():
        key = (min(idx[qf0], idx[qf1]), max(idx[qf0], idx[qf1]))
        if key not in seen:
            seen.add(key)
            edges.append({'s': key[0], 't': key[1], 'kind': 'v'})
    return {'disc': cls.disc, 'cond': cls.cond, 'l': l,
            'nodes': nodes, 'edges': edges}
