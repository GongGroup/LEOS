"""In-situ Allred-Rochow ownership with a shape-based NELECT cap (table-free).

r3 repair of r2-T-B (parent sha256 d58e27ae...). The per-WF polarity mechanism is
the parent's (Gauss screening + in-situ chi + force fallback); what changed:

1. Occupations (change 2). A WF holds its #WF occupation when the block covers
   every WF (n_pair electrons otherwise); the inner-electron count and the valence
   second moment are weighted by it, so an empty or partly empty WF neither
   screens nor defines a radius.

2. Dipole-free shells (change 3). An own WF (one-centre, nucleus within its rms
   sphere) whose centroid offset is below PARITY = 1/8 of its rms radius is a
   parity-pure d/s-like shell. It screens a contested WF at distance d by the
   exact fraction of its density (isotropic Gaussian, per-axis variance
   spread/3) inside the sphere of radius d, and it never defines the valence
   radius. Only the directional (hybrid) own WFs are lone pairs: wholly inside
   (offset + REACH_SIGMA * rms < d) they screen fully, otherwise they give
   rho^2 = <offset^2 + spread> and chi = Z_eff / rho^2 decides between two ends
   that both carry lone pairs; otherwise the force Z_eff / d^2 decides
   (Z_eff = ZVAL - inner electrons - 5/16).

3. One WF extent (change 4). REACH_SIGMA = 1 rms radius is declared and used in
   the three places where the parent had an implicit 1: the anchoring test, the
   inner-shell test and the atom's cloud (held test).

4. Homonuclear shares (change 6). Atoms of the winning element whose pairwise
   bond regions contain the WF share it in proportion to tau / d^2, tau = depth
   of the WF inside the deepest bond region the member forms with another
   member (1 at the midpoint, 0 on the region's boundary). Two-member groups
   keep the parent's 1/d^2 split; the far vertex of an equilateral triangle
   (dm = da exactly, the knife edge) enters with a share that goes to zero
   continuously.

5. NELECT surplus (change 1). With occupations: the parent's least-occupied
   levels; a partly emptied occupation-degenerate level is emptied equally per
   WF. Without occupations: WFs held by no atom are emptied first; then the
   donors are the atoms that hold a dipole-free own WF (centred d/s-like states,
   where a projected manifold puts its empty states); their charges are
   equalised (minimum sum q^2), each giving at most the electrons of its own
   WFs, which inside the atom are emptied dipole-free first, then most diffuse
   first. The loss is averaged within discrete classes (element + multiset of
   the kinds of WFs held, spin channels unordered). No site electrostatics
   decides anything. The answer is flagged 'undetermined' (diagnose()) when no
   dipole-free WF can take the surplus (the rest then comes from the most
   diffuse WFs), when every atom has the same full-supply charge (FeO type) or
   ends at one charge, and when the parent's electrostatic removal, run as a
   strict witness (no all-atoms fallback, all-neutral outcome rejected), is
   missing or gives other integers.

6. Integers (change 5). Largest remainders with the electron total conserved;
   a tied group straddling the cut is split one electron per bonded fragment
   in turn, canonical order (wrapped fractional coordinates) inside and
   between fragments, and its mean is reported by diagnose(); a tie spanning
   several classes is flagged. assign() asserts
   sum(OS) = sum ZVAL - min(NELECT, n_pair * number of WFs).

No electronegativity or oxidation-state table is read.
"""
from __future__ import annotations

import math
import sys

from pyleos.parsers import inv3

PAIR_SCREEN = 5.0 / 16.0   # electrons: partner screening inside one pair
FLOAT_EPS = 1e-6           # numerical equality (occupations printed with 6 decimals)
REACH_SIGMA = 1.0          # extent of a WF, in rms radii
PARITY = 1.0 / 8.0         # dipole-free own WF: centroid offset below PARITY * rms
# derived, not free: Ewald cutoffs converged to machine epsilon; Gaussian helper
EWALD_X = math.sqrt(-math.log(sys.float_info.epsilon))
SQRT2 = math.sqrt(2.0)

META = {
    "name": "In-situ Allred-Rochow ownership with dipole-free shells, depth-weighted "
            "homonuclear shares and a shape-based NELECT cap",
    "mechanism": (
        "One-centre WFs stay on their atom. For a contested WF at distance d from an "
        "atom, the atom's own WFs screen its ZVAL: dipole-free own WFs (centroid offset "
        "< PARITY*rms: d/s-like shells) by the Gaussian fraction of their density inside "
        "the sphere of radius d, hybrid own WFs fully when offset + REACH_SIGMA*rms < d; "
        "the remaining hybrid own WFs are lone pairs whose occupation-weighted second "
        "moment is the valence radius rho^2. Between ends that both carry lone pairs the "
        "larger in-situ Allred-Rochow chi = Z_eff/rho^2 takes the pair, otherwise the "
        "larger Z_eff/d^2 (Z_eff = ZVAL - inner electrons - 5/16). Same-element ends "
        "share in proportion to depth/d^2 (depth of the WF in their bond region). NELECT "
        "surplus: least-occupied levels when occupations exist; otherwise unheld WFs, "
        "then charge equalisation over the atoms holding dipole-free WFs (capped by "
        "their own WFs, emptied dipole-free then most diffuse first), averaged over "
        "discrete classes; ambiguous cells are flagged 'undetermined' by diagnose() "
        "(no dipole-free WF to empty, FeO-type equal charges, or disagreement with the "
        "parent's electrostatic removal run as a witness). Largest remainders with the "
        "electron total conserved. No electronegativity or oxidation-state table."),
    "params": {
        "PAIR_SCREEN": [PAIR_SCREEN,
                        "screening of one electron by its partner in the same orbital: "
                        "variational He-like result Z - 5/16; inherited from r1-T6/r2-T-B, "
                        "not refitted (dev plateau 0.1-0.6, see stage2/sweep_s2.log)"],
        "FLOAT_EPS": [FLOAT_EPS,
                      "float hygiene: equality of occupations (printed with 6 decimals), "
                      "remainders, spreads and site-potential signs; inherited"],
        "REACH_SIGMA": [REACH_SIGMA,
                        "extent of a WF in rms radii, used for anchoring (nucleus within "
                        "one extent of the centre), the inner-shell test (offset + extent "
                        "< d) and the atom's cloud; one rms radius is the natural size of "
                        "a Gaussian-like WF (implicit in the parent, declared now); dev "
                        "plateau in stage2/sweep_s2.log"],
        "PARITY": [PARITY,
                   "dipole-free test offset < PARITY*rms: a d (or s) function has no "
                   "same-centre dipole by parity and is displaced only by inter-site "
                   "mixing, while an s-p hybrid on one centre is displaced by 0.3-0.5 rms "
                   "(hydrogenic sp3 0.51); 1/8 rms ~ 2 percent p admixture (hydrogenic "
                   "estimate a = 0.92 sqrt(eps)). Dev plateau 0.08-0.17 "
                   "(stage2/sweep_s2.log); set between the dev d shells (<= 0.078) and "
                   "lone pairs (>= 0.178) from unlabelled geometry, not from labels"],
    },
    "derived": {
        "EWALD_X": [EWALD_X, "sqrt(-ln(machine epsilon)) from sys.float_info: Ewald "
                    "cutoffs of the witness converged to double precision"],
    },
    "table_free": True,
}


# ---------------------------------------------------------------- geometry
def _min_image(fc, fa, latt):
    """Cartesian vector from point fc to atom fa (fractional) and its |.|^2."""
    df = [fa[k] - fc[k] for k in range(3)]
    df = [x - round(x) for x in df]
    best, bv = math.inf, None
    for n0 in (-1, 0, 1):
        for n1 in (-1, 0, 1):
            for n2 in (-1, 0, 1):
                f = (df[0] + n0, df[1] + n1, df[2] + n2)
                v = [f[0] * latt[0][j] + f[1] * latt[1][j] + f[2] * latt[2][j]
                     for j in range(3)]
                d2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2]
                if d2 < best:
                    best, bv = d2, v
    return bv, best


def _n2(v):
    return v[0] * v[0] + v[1] * v[1] + v[2] * v[2]


def _mid_d2(va, vb):
    return _n2([(va[k] + vb[k]) / 2 for k in range(3)])


def _in_bond_region(va, vb):
    """WF centre (origin) is nearer the a-b midpoint than either nucleus."""
    dm = _mid_d2(va, vb)
    return dm < _n2(va) and dm < _n2(vb)


def _depth(va, vb):
    """Depth of the WF (origin) inside the a-b bond region: 1 at the midpoint,
    0 on the boundary (midpoint as far as the nearer nucleus)."""
    dn = math.sqrt(min(_n2(va), _n2(vb)))
    if dn <= 0.0:
        return 1.0
    return max(0.0, 1.0 - math.sqrt(_mid_d2(va, vb)) / dn)


def _gauss_inside(o, s, R):
    """Fraction of an isotropic Gaussian (per-axis sigma s) centred at distance
    o from the nucleus that lies inside the sphere of radius R."""
    if R <= 0.0:
        return 0.0
    if s <= 0.0:
        return 1.0 if o < R else 0.0
    if o <= FLOAT_EPS * s:
        x = R / s
        return math.erf(x / SQRT2) - math.sqrt(2.0 / math.pi) * x * math.exp(-0.5 * x * x)
    a, b = (R - o) / s, (R + o) / s
    v = 0.5 * (math.erf(a / SQRT2) + math.erf(b / SQRT2)) \
        - s / (o * math.sqrt(2.0 * math.pi)) * (math.exp(-0.5 * a * a) - math.exp(-0.5 * b * b))
    return min(1.0, max(0.0, v))


def _cross(u, v):
    return [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]


def _dot(u, v):
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


# ---------------------------------------------------------------- electrostatic witness
def site_potentials(latt, frac, q):
    """Ewald potential (e/A) at every nucleus from the point charges q of all
    other atoms and images (uniform background for a net charge)."""
    n = len(q)
    a, b, c = latt
    vol = abs(_dot(a, _cross(b, c)))
    rl = [[2.0 * math.pi / vol * x for x in _cross(b, c)],
          [2.0 * math.pi / vol * x for x in _cross(c, a)],
          [2.0 * math.pi / vol * x for x in _cross(a, b)]]
    eta = math.sqrt(math.pi) * (n / (vol * vol)) ** (1.0 / 6.0)
    rcut = EWALD_X / eta
    gcut = 2.0 * eta * EWALD_X

    def reach(vecs, cut):
        out = []
        for i in range(3):
            nrm = _cross(vecs[(i + 1) % 3], vecs[(i + 2) % 3])
            h = abs(_dot(vecs[i], nrm)) / math.sqrt(_dot(nrm, nrm))
            out.append(int(math.ceil(cut / h)) + 1)
        return out

    cart = [[f[0] * a[j] + f[1] * b[j] + f[2] * c[j] for j in range(3)] for f in frac]
    nr, ng = reach(latt, rcut), reach(rl, gcut)
    shifts = [[i * a[k] + j * b[k] + m * c[k] for k in range(3)]
              for i in range(-nr[0], nr[0] + 1)
              for j in range(-nr[1], nr[1] + 1)
              for m in range(-nr[2], nr[2] + 1)]
    V = [0.0] * n
    rc2 = rcut * rcut
    for i in range(n):
        s = 0.0
        for j in range(n):
            if q[j] == 0.0:
                continue
            d0 = [cart[i][k] - cart[j][k] for k in range(3)]
            for t in shifts:
                x, y, z = d0[0] - t[0], d0[1] - t[1], d0[2] - t[2]
                r2 = x * x + y * y + z * z
                if 0.0 < r2 < rc2:
                    r = math.sqrt(r2)
                    s += q[j] * math.erfc(eta * r) / r
        V[i] = s
    for i0 in range(-ng[0], ng[0] + 1):
        for j0 in range(-ng[1], ng[1] + 1):
            for m0 in range(-ng[2], ng[2] + 1):
                if i0 == 0 and j0 == 0 and m0 == 0:
                    continue
                g = [i0 * rl[0][k] + j0 * rl[1][k] + m0 * rl[2][k] for k in range(3)]
                g2 = _dot(g, g)
                if g2 > gcut * gcut:
                    continue
                pref = 4.0 * math.pi / vol * math.exp(-g2 / (4.0 * eta * eta)) / g2
                ph = [_dot(g, x) for x in cart]
                sr = sum(qj * math.cos(p) for qj, p in zip(q, ph))
                si = sum(qj * math.sin(p) for qj, p in zip(q, ph))
                for i in range(n):
                    V[i] += pref * (math.cos(ph[i]) * sr + math.sin(ph[i]) * si)
    qt = sum(q)
    for i in range(n):
        V[i] -= 2.0 * eta / math.sqrt(math.pi) * q[i] + math.pi * qt / (vol * eta * eta)
    return V


def _equalise(q, cap, members, amount):
    """Remove `amount` electrons from `members`, raising their charges to a
    common level (minimum sum of q^2), atom a giving at most cap[a]; exact
    piecewise-linear solution."""
    if amount <= 0.0 or not members:
        return {}

    def removed(level):
        return sum(min(max(level - q[a], 0.0), cap[a]) for a in members)

    knots = sorted(set([q[a] for a in members] + [q[a] + cap[a] for a in members]))
    lev = knots[-1]
    for lo, hi in zip(knots, knots[1:]):
        r_hi = removed(hi)
        if r_hi >= amount:
            r_lo = removed(lo)
            lev = lo + (hi - lo) * (amount - r_lo) / (r_hi - r_lo) if r_hi > r_lo else hi
            break
    return {a: min(max(lev - q[a], 0.0), cap[a]) for a in members}


def _witness_removal(st, q, cap, amount):
    """The parent's electrostatic removal, strict: donors are the atoms at
    negative site potential once the removal is made (fixed point started
    from all atoms), charges equalised. None when no site is negative, when
    the donors cannot take the surplus, or when every atom ends at one charge
    (the all-neutral fixed point carries no information)."""
    cand = tuple(a for a in range(len(q)) if cap[a] > FLOAT_EPS)
    donors, seen, x = cand, [], {}
    qx = list(q)
    for _ in range(len(cand) + 2):
        take = min(amount, sum(cap[a] for a in donors))
        if amount - take > FLOAT_EPS:
            return None
        x = _equalise(q, cap, list(donors), take)
        qx = [q[a] + x.get(a, 0.0) for a in range(len(q))]
        V = site_potentials(st.lattice, st.frac, qx)
        new = tuple(a for a in cand if V[a] < -FLOAT_EPS)
        if not new:
            return None
        if new == donors or new in seen:
            break
        seen.append(donors)
        donors = new
    if max(qx) - min(qx) <= FLOAT_EPS:
        return None
    return x


# ---------------------------------------------------------------- shells
class _Shells:
    """Own (one-centre, anchored) WFs of one atom: rows (offset, rms radius,
    electrons, dipole-free)."""

    def __init__(self, rows):
        self.pp = [(o, r / math.sqrt(3.0), n) for o, r, n, p in rows if p]
        self.hy = [(o, r, n) for o, r, n, p in rows if not p]

    def seen_from(self, d):
        """(inner electrons, valence radius^2 or None) for a WF at distance d."""
        inner = 0.0
        for o, s, n in self.pp:
            inner += n * _gauss_inside(o, s, d)
        w = m = 0.0
        for o, r, n in self.hy:
            if o + REACH_SIGMA * r < d:
                inner += n
            else:
                w += n
                m += n * (o * o + r * r)
        return inner, (m / w if w > FLOAT_EPS else None)


# ---------------------------------------------------------------- ownership
def _analyse(case):
    st = case.structure
    latt = st.lattice
    els = st.elements
    na = len(els)
    zval = case.eoscar.zval
    binv = inv3(latt)
    n_pair = 2 if case.ispin == 1 else 1

    occ = {}
    for wd in (case.eoscar.wf_detail or []):
        occ[(wd.spin, wd.index)] = wd.occ
    n_wf = sum(len(ch.centres) for ch in case.channels)
    has_occ = n_wf > 0 and all((ci + 1, wi + 1) in occ
                               for ci, ch in enumerate(case.channels)
                               for wi in range(len(ch.centres)))

    geo = []
    own = [[] for _ in range(na)]
    for ci, ch in enumerate(case.channels):
        for wi, centre in enumerate(ch.centres):
            fc = [sum(centre[i] * binv[i][j] for i in range(3)) for j in range(3)]
            vecs, d2 = [], []
            for a in range(na):
                v, dd = _min_image(fc, st.frac[a], latt)
                vecs.append(v)
                d2.append(dd)
            near = min(range(na), key=lambda a: (d2[a], a))
            contested = any(a != near and _in_bond_region(vecs[near], vecs[a])
                            for a in range(na))
            spread = max(ch.spreads[wi], 0.0)
            rad = math.sqrt(spread)
            if has_occ:
                e = max(0.0, min(float(n_pair), occ[(ci + 1, wi + 1)]))
            else:
                e = float(n_pair)
            off = math.sqrt(d2[near])
            is_own = (not contested) and off <= REACH_SIGMA * rad
            pp = is_own and off < PARITY * rad
            if is_own:
                own[near].append((off, rad, e, pp))
            geo.append({"vecs": vecs, "d2": d2, "near": near, "contested": contested,
                        "rad": rad, "spread": spread, "e": e, "own": is_own, "pp": pp,
                        "slot": e / n_pair, "ch": ci})
    shells = [_Shells(rows) for rows in own]

    notes = []
    for g in geo:
        vecs, d2, near = g["vecs"], g["d2"], g["near"]
        if not g["contested"]:
            g["hosts"] = {near: 1.0}
            continue
        zeff, rho2, weak = [], [], []
        for a in range(na):
            inner, r2 = shells[a].seen_from(math.sqrt(d2[a]))
            z = zval[els[a]] - inner - PAIR_SCREEN
            zeff.append(z)
            rho2.append(r2)
            weak.append(d2[a] / z if z > 0.0 else math.inf)   # 1 / force
        order = sorted(range(na), key=lambda a: (weak[a], d2[a], a))
        top = order[0]
        partners = [b for b in order if els[b] != els[top]
                    and _in_bond_region(vecs[top], vecs[b])]
        if rho2[top] is not None:
            best = top
            for b in partners:
                if rho2[b] is not None and zeff[b] / rho2[b] > zeff[best] / rho2[best]:
                    best = b
            if best != top:
                top = best
                partners = [b for b in order if els[b] != els[top]
                            and _in_bond_region(vecs[top], vecs[b])]
        group = [top]
        for a in order:
            if a == top or els[a] != els[top]:
                continue
            if all(_in_bond_region(vecs[a], vecs[b]) for b in group):
                group.append(a)
            else:
                break
        if len(group) >= 2:
            tau = {a: max(_depth(vecs[a], vecs[b]) for b in group if b != a) for a in group}
            force = {a: tau[a] / d2[a] for a in group}
            tot = sum(force.values())
            g["hosts"] = {a: f / tot for a, f in force.items()}
        else:
            g["hosts"] = {top: 1.0}
        if partners:
            b = partners[0]
            both = rho2[top] is not None and rho2[b] is not None
            notes.append({
                "winner": els[top], "loser": els[b], "rule": "chi" if both else "force",
                "force_margin": (None if weak[top] == math.inf else
                                 math.inf if weak[top] == 0.0 else weak[b] / weak[top]),
                "chi_margin": ((zeff[top] / rho2[top]) / (zeff[b] / rho2[b])
                               if both and zeff[b] > 0.0 else None),
                "zeff": (zeff[top], zeff[b]),
                "d": (math.sqrt(d2[top]), math.sqrt(d2[b])),
            })

    # held test: an atom's cloud reaches as far as its anchored WFs do
    cloud = [0.0] * na
    for g in geo:
        for a in range(na):
            r = math.sqrt(g["d2"][a])
            if r <= REACH_SIGMA * g["rad"]:
                cloud[a] = max(cloud[a], r + REACH_SIGMA * g["rad"])
    for g in geo:
        g["held"] = any(g["d2"][a] <= cloud[a] * cloud[a] for a in range(na))
    return {"geo": geo, "notes": notes, "has_occ": has_occ, "n_pair": n_pair,
            "na": na, "els": els, "zval": zval, "st": st}


def pair_margins(case):
    """Diagnostics per (winner element, loser element): contested WFs, how many
    decided by chi or by force, and the smallest margins (winner/loser)."""
    out = {}
    for n in _analyse(case)["notes"]:
        rec = out.setdefault((n["winner"], n["loser"]),
                             {"n": 0, "by_chi": 0, "by_force": 0,
                              "min_chi_margin": None, "min_force_margin": None})
        rec["n"] += 1
        if n["rule"] == "chi":
            rec["by_chi"] += 1
            m = n["chi_margin"]
            if m is not None and (rec["min_chi_margin"] is None or m < rec["min_chi_margin"]):
                rec["min_chi_margin"] = m
        else:
            rec["by_force"] += 1
            m = n["force_margin"]
            if m is not None and (rec["min_force_margin"] is None
                                  or m < rec["min_force_margin"]):
                rec["min_force_margin"] = m
    return out


# ---------------------------------------------------------------- surplus
def _classes(an):
    """Discrete classes: element + multiset (spin channels unordered) of the
    kinds of WFs an atom holds: dipole-free own, hybrid own, k-centre share,
    heteronuclear WF won, other."""
    na, els = an["na"], an["els"]
    tok = [dict() for _ in range(na)]
    for g in an["geo"]:
        h = g["hosts"]
        for a in h:
            if len(h) > 1:
                t = ("h", len(h))
            elif g["own"]:
                t = ("p",) if g["pp"] else ("l",)
            elif g["contested"]:
                t = ("w",)
            else:
                t = ("u",)
            tok[a].setdefault(g["ch"], []).append(t)
    groups = {}
    for a in range(na):
        sig = (els[a], tuple(sorted(tuple(sorted(v)) for v in tok[a].values())))
        groups.setdefault(sig, []).append(a)
    return sorted(groups.values())


def _pool(lost, classes):
    """Atoms of one class lose the class mean."""
    out = list(lost)
    for cls in classes:
        if len(cls) > 1:
            m = sum(lost[a] for a in cls) / len(cls)
            for a in cls:
                out[a] = m
    return out


def _empty_wfs(wfs, amount, lost):
    """Remove `amount` electrons from the listed WFs, equally per electron
    held; the loss goes to each WF's hosts. Returns what could not be taken."""
    if amount <= 0.0 or not wfs:
        return amount
    have = sum(g["left"] for g in wfs)
    if have <= FLOAT_EPS:
        return amount
    take = min(amount, have)
    for g in wfs:
        f = take * g["left"] / have
        for a, s in g["hosts"].items():
            lost[a] += f * s
        g["left"] -= f
    return amount - take


def _tiers(wfs, key):
    """WFs sorted by key (a tuple whose last entry is a float); consecutive
    WFs with the same leading entries and the last entry equal within
    FLOAT_EPS form one tier, emptied together."""
    out = []
    for g in sorted(wfs, key=key):
        if out:
            k0, k1 = key(out[-1][0]), key(g)
            if k0[:-1] == k1[:-1] and abs(k0[-1] - k1[-1]) <= FLOAT_EPS:
                out[-1].append(g)
                continue
        out.append([g])
    return out


def _shape_removal(an, ne, amount, diag):
    """No occupations: per-atom loss from the shape rule (see module doc, 5)."""
    geo, na, els, zval = an["geo"], an["na"], an["els"], an["zval"]
    lost = [0.0] * na
    q = [zval[els[a]] - ne[a] for a in range(na)]
    donors = sorted({g["near"] for g in geo if g["own"] and g["pp"] and g["left"] > FLOAT_EPS})
    diag["donors"] = donors
    if donors:
        cap = [0.0] * na
        for g in geo:
            if g["own"]:
                cap[g["near"]] += g["left"]
        x = _equalise(q, cap, donors, amount)
        for a in donors:
            v = x.get(a, 0.0)
            if v <= FLOAT_EPS:
                continue
            mine = [g for g in geo if g["own"] and g["near"] == a and g["left"] > FLOAT_EPS]
            for tier in _tiers(mine, lambda g: (not g["pp"], -g["spread"])):
                v = _empty_wfs(tier, v, lost)
                if v <= FLOAT_EPS:
                    break
        amount -= sum(lost)
    if amount > FLOAT_EPS:
        diag["flags"].append("undetermined: surplus beyond the dipole-free WFs; "
                             "taken from the most diffuse WFs")
        diag["undetermined"] = True
        rest = [g for g in geo if g["left"] > FLOAT_EPS]
        for tier in _tiers(rest, lambda g: (-g["spread"],)):
            amount = _empty_wfs(tier, amount, lost)
            if amount <= FLOAT_EPS:
                break
    return lost


def _populations(case, diag):
    an = _analyse(case)
    geo, na, els, zval, n_pair = an["geo"], an["na"], an["els"], an["zval"], an["n_pair"]
    ne = [0.0] * na
    for g in geo:
        for a, f in g["hosts"].items():
            ne[a] += n_pair * f
        g["left"] = float(n_pair)     # electrons still in the WF during the cap
    diag["mode"] = "closed"
    if case.eoscar.nelect_wann is None:
        return ne, an
    excess = n_pair * len(geo) - int(round(case.eoscar.nelect_wann))
    diag["excess"] = excess
    if excess <= 0:
        return ne, an
    lost = [0.0] * na
    if an["has_occ"]:
        # occupation levels, least occupied first; unheld WFs first in a level
        diag["mode"] = "occupations"
        levels = []
        for g in sorted(geo, key=lambda g: (g["slot"], g["held"])):
            if levels and abs(levels[-1][-1]["slot"] - g["slot"]) <= FLOAT_EPS \
                    and levels[-1][-1]["held"] == g["held"]:
                levels[-1].append(g)
            else:
                levels.append([g])
        for lev in levels:
            if excess <= FLOAT_EPS:
                break
            size = n_pair * len(lev)
            if excess >= size - FLOAT_EPS:
                _empty_wfs(lev, size, lost)
                excess -= size
                continue
            if len(lev) > 1:
                diag["flags"].append("occupation-degenerate level shared equally")
            _empty_wfs(lev, excess, lost)
            excess = 0
        for a in range(na):
            ne[a] -= lost[a]
        return ne, an
    # ------------------------------------------------ no occupations
    diag["mode"] = "no occupations"
    unheld = [g for g in geo if not g["held"]]
    if unheld:
        before = excess
        excess = _empty_wfs(unheld, excess, lost)
        diag["unheld_removed"] = before - excess
    ne0 = [ne[a] - lost[a] for a in range(na)]
    classes = _classes(an)
    if excess > FLOAT_EPS:
        q0 = [zval[els[a]] - ne0[a] for a in range(na)]
        if len(classes) > 1 and max(q0) - min(q0) <= FLOAT_EPS * na:
            diag["flags"].append("undetermined: every atom at the same full-supply charge")
            diag["undetermined"] = True
        x = _pool(_shape_removal(an, ne0, excess, diag), classes)
        if len(classes) > 1:
            cap = [max(0.0, ne0[a]) for a in range(na)]
            w = _witness_removal(an["st"], q0, cap, excess)
            if w is None:
                diag["witness"] = None
            else:
                w = _pool([w.get(a, 0.0) for a in range(na)], classes)
                diag["witness"] = [ne0[a] - w[a] for a in range(na)]
        for a in range(na):
            lost[a] += x[a]
        qf = [zval[els[a]] - ne[a] + lost[a] for a in range(na)]
        if len(classes) > 1 and max(qf) - min(qf) <= FLOAT_EPS * na:
            diag["flags"].append("undetermined: every atom ends at the same charge")
            diag["undetermined"] = True
    diag["lost"] = lost
    for a in range(na):
        ne[a] -= lost[a]
    return ne, an


def _fragments(an):
    """Bonded fragment of every atom: atoms linked by WFs with two or more
    hosts (union-find)."""
    parent = list(range(an["na"]))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for g in an["geo"]:
        h = sorted(a for a, f in g["hosts"].items() if f > FLOAT_EPS)
        for a in h[1:]:
            ra, rb = find(h[0]), find(a)
            if ra != rb:
                parent[max(ra, rb)] = min(ra, rb)
    return [find(a) for a in range(an["na"])]


def _apportion(ne, frac, diag, classes=None, frag=None):
    """Integer populations summing to round(sum(ne)) by largest remainders.
    Atoms with equal remainders form a tie group; a group that straddles the
    cut is split one electron per bonded fragment in turn (so that each O2
    unit of a superoxide gets one), canonical order (wrapped fractional
    coordinates) inside and between fragments, and reported with its mean."""
    total = int(round(sum(ne)))
    base = [int(math.floor(x + FLOAT_EPS)) for x in ne]
    rest = total - sum(base)
    rem = [x - b for x, b in zip(ne, base)]
    order = sorted(range(len(ne)), key=lambda a: -rem[a])
    groups = []
    for a in order:
        if groups and abs(rem[groups[-1][0]] - rem[a]) < FLOAT_EPS:
            groups[-1].append(a)
        else:
            groups.append([a])

    def canon(a):
        return tuple(round(f % 1.0, 9) % 1.0 for f in frac[a])

    def pick_spread(g, k):
        by = {}
        for a in sorted(g, key=canon):
            by.setdefault(frag[a] if frag else a, []).append(a)
        queues = sorted(by.values(), key=lambda m: canon(m[0]))
        pick = []
        while len(pick) < k:
            for qu in queues:
                if qu and len(pick) < k:
                    pick.append(qu.pop(0))
        return pick

    cl = {a: i for i, c in enumerate(classes or []) for a in c}
    out = list(base)
    if rest < 0:
        for g in reversed(groups):
            if rest >= 0:
                break
            k = min(len(g), -rest)
            for a in pick_spread(g, k):
                out[a] -= 1
            rest += k
    for g in groups:
        if rest <= 0:
            break
        if len(g) <= rest:
            for a in g:
                out[a] += 1
            rest -= len(g)
        else:
            pick = pick_spread(g, rest)
            for a in pick:
                out[a] += 1
            diag.setdefault("ties", []).append({
                "atoms": sorted(g), "mean_electrons": sum(ne[a] for a in g) / len(g),
                "promoted": sorted(pick), "classes": len({cl.get(a, -1 - a) for a in g})})
            rest = 0
    return out


def _run(case):
    diag = {"flags": [], "undetermined": False}
    ne, an = _populations(case, diag)
    classes = _classes(an)
    frag = _fragments(an)
    n = _apportion(ne, an["st"].frac, diag, classes, frag)
    zval, els = an["zval"], an["els"]
    os_ = [int(round(zval[e])) - k for e, k in zip(els, n)]
    if "witness" in diag:
        w = diag["witness"]
        if w is None:
            diag["flags"].append("undetermined: no electrostatic witness "
                                 "(no negative site, or the all-neutral fixed point)")
            diag["undetermined"] = True
        else:
            wn = _apportion(w, an["st"].frac, {}, None, frag)
            diag["witness_os"] = [int(round(zval[e])) - k for e, k in zip(els, wn)]
            if diag["witness_os"] != os_:
                diag["flags"].append("undetermined: the electrostatic witness gives "
                                     "other integers")
                diag["undetermined"] = True
    for t in diag.get("ties", []):
        if t["classes"] > 1:
            diag["flags"].append("undetermined: a tie at the rounding cut spans "
                                 "several classes")
            diag["undetermined"] = True
    supply = an["n_pair"] * len(an["geo"])
    placed = supply if case.eoscar.nelect_wann is None else \
        min(supply, int(round(case.eoscar.nelect_wann)))
    target = int(round(sum(zval[e] for e in els))) - placed
    assert sum(os_) == target, (sum(os_), target)
    diag["sum_os"], diag["cell_charge"] = sum(os_), target
    diag["populations"] = ne
    return os_, diag


def electrons(case):
    """Fractional electron population per atom (POSCAR order), after the cap."""
    return _run(case)[1]["populations"]


def assign(case):
    return _run(case)[0]


def diagnose(case):
    """Diagnostics: mode (closed / occupations / no occupations), excess,
    'undetermined' (bool) and the flags that set it, donors of the shape rule,
    the witness integers, tied groups with their mean, sum(OS) and the cell
    charge, fractional populations."""
    return _run(case)[1]
