"""Wannier-centre-offset ownership (table-free). Template for all variants.

Ownership rule (the only thing that differs from the scaffold): a contested
(bond) Wannier function (WF) between its nearest atom A and the atoms B of
another element that share a bond region with A goes to the candidate with the
smallest effective distance from the WF centre. The effective distance encodes
where the neutral reference point of the bond lies:

  V0   midpoint            eff = d                      (Abu-Farsakh & Qteish)
  V1a  radius ratio        eff = d / r, r = rms radius of the atom's own WFs
                           (lone pairs / core; r^2 = <offset^2 + spread>);
                           if any candidate has no own WF: midpoint
  V1b  radius ratio, zero  as V1a, but an atom without own WFs has r = 0 and
       size                loses to every candidate that has own WFs
  V1d  lone-pair radius    as V1b, but r is taken from the directional own
       ratio, zero size    WFs (lone pairs) only; dipole-free shells (d, s) give
                           no valence radius
  V1c  hybrid centroid     eff = d - a, a = mean centre offset of the atom's
                           own directional (lone-pair) WFs, 0 if it has none:
                           the reference is the midpoint between the two hybrid
                           centroids, t_ref = 1/2 + (a_A - a_B) / (2 d_AB)
  V3   V1c + core force     as V1c; a bond in which no candidate has a lone
                           pair (a = 0 for all) has no in-situ size, and its
                           reference is the point where the pulls of the two
                           bare cores balance, Z_A/d_A^2 = Z_B/d_B^2, i.e.
                           eff = d / sqrt(Z), Z = ZVAL - electrons of the
                           atom's own (dipole-free) WFs. Hybrid with the force
                           idea of r3-T-B, without its screening constant
  V3s  V3, pair screening  as V3 with Z reduced by PAIR_SCREEN = 5/16 (the
                           constant of r3-T-B, inherited, not tuned here)
  V2   pair-pooled V0      per case, the normalised offsets (eff_B - eff_A) /
                           (d_A + d_B) of all bond WFs of an element pair are
                           averaged; the pair gets one owner element, cycles
                           between pairs are broken at the weakest pair, and
                           every bond WF of the pair is assigned that way
  V2c  pair-pooled V1c

Bond type: among the atoms that share a bond region with A, the one whose
midpoint with A is closest to the WF centre decides; if it is of A's element
the WF is a homonuclear bond (shared, no other element competes).

Scaffold (re-implemented from pyleos/methods/frozen/r3-T-B-cap-repair.py, see
the report for the exact list): minimum-image geometry, bond-region test
(centre nearer the A-B midpoint than either nucleus), own-WF test (not
contested, nucleus within REACH_SIGMA rms radii), dipole-free test (PARITY),
homonuclear shares depth/d^2, NELECT surplus rule (occupation levels, else
unheld WFs then charge equalisation over atoms holding dipole-free WFs, class
pooling), largest-remainder integers. The Ewald witness of r3-T-B (diagnostic
only) is not copied. No electronegativity table and no oxidation-state bounds.
"""
from __future__ import annotations

import math

from pyleos.parsers import inv3

VARIANT = "V1c"

FLOAT_EPS = 1e-6           # numerical equality (inherited from r3-T-B)
REACH_SIGMA = 1.0          # extent of a WF in rms radii (inherited from r3-T-B)
PARITY = 1.0 / 8.0         # dipole-free own WF: offset < PARITY * rms (inherited)
PAIR_SCREEN = 5.0 / 16.0   # used by V3s only (inherited from r3-T-B)

META = {
    "name": "WFC-offset ownership, variant " + VARIANT,
    "mechanism": "bond WF goes to the end with the smaller effective distance; "
                 "variant fixes the neutral reference point (see module doc)",
    "params": {
        "FLOAT_EPS": [FLOAT_EPS, "float hygiene, inherited from r3-T-B, not tuned"],
        "REACH_SIGMA": [REACH_SIGMA, "own-WF / held tests, inherited from r3-T-B, not tuned"],
        "PARITY": [PARITY, "dipole-free test, inherited from r3-T-B, not tuned"],
        "PAIR_SCREEN": [PAIR_SCREEN, "variant V3s only: screening of one electron by its "
                        "partner (He-like 5/16), inherited from r3-T-B, not tuned; "
                        "dev plateau 0.05-0.60 (sensitivity_dev.py)"],
    },
    "table_free": True,
}


# ---------------------------------------------------------------- geometry
def _min_image(fc, fa, latt):
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
    dm = _mid_d2(va, vb)
    return dm < _n2(va) and dm < _n2(vb)


def _depth(va, vb):
    dn = math.sqrt(min(_n2(va), _n2(vb)))
    if dn <= 0.0:
        return 1.0
    return max(0.0, 1.0 - math.sqrt(_mid_d2(va, vb)) / dn)


def _equalise(q, cap, members, amount):
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


# ---------------------------------------------------------------- in-situ sizes
def _atom_sizes(own, na):
    """Per atom, from its own WFs (offset, rms radius, electrons, dipole-free):
    r = rms radius of all own WFs about the nucleus (None without own WFs),
    a = electron-weighted mean offset of the directional own WFs (0.0 if none),
    rh = rms radius of the directional own WFs only (None if none)."""
    r, a, rh = [None] * na, [0.0] * na, [None] * na
    for i in range(na):
        w = m = 0.0
        wh = mh = m2 = 0.0
        for o, rad, e, pp in own[i]:
            w += e
            m += e * (o * o + rad * rad)
            if not pp:
                wh += e
                mh += e * o
                m2 += e * (o * o + rad * rad)
        if w > FLOAT_EPS:
            r[i] = math.sqrt(m / w)
        if wh > FLOAT_EPS:
            a[i] = mh / wh
            rh[i] = math.sqrt(m2 / wh)
    return r, a, rh


def _eff(variant, d, r, a, any_undefined):
    """Effective distance of a candidate at distance d from the WF centre."""
    if variant in ("V0", "V2"):
        return d
    if variant in ("V1c", "V2c"):
        return d - a
    if variant == "V1a":
        return d if any_undefined else d / r
    if variant in ("V1b", "V1d"):
        if any_undefined == "all":
            return d
        return math.inf if r is None else d / r
    raise ValueError(variant)


# ---------------------------------------------------------------- ownership
def _analyse(case, variant=None, owner_hook=None):
    variant = variant or VARIANT
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
                        "slot": e / n_pair, "ch": ci, "wi": wi})
    rr, aa, rh = _atom_sizes(own, na)
    rsel = rh if variant == "V1d" else rr

    # candidates of every contested WF: nearest atom + atoms of another element
    # that share a bond region with it (deepest bond region first)
    for g in geo:
        if not g["contested"]:
            continue
        near, vecs = g["near"], g["vecs"]
        region = sorted([b for b in range(na) if b != near
                         and _in_bond_region(vecs[near], vecs[b])],
                        key=lambda b: (_mid_d2(vecs[near], vecs[b]), b))
        # bond type from the deepest bond region: with a same-element atom the
        # WF is a homonuclear bond and no other element competes for it
        g["homo"] = els[region[0]] == els[near]
        g["partners"] = [] if g["homo"] else [b for b in region if els[b] != els[near]]
        cands = [near] + g["partners"]
        undef = [rsel[c] is None for c in cands]
        flag = "all" if all(undef) else (True if any(undef) else False)
        if variant in ("V1b", "V1d"):
            flag = "all" if all(undef) else False
        if variant in ("V3", "V3s"):
            if any(aa[c] > 0.0 for c in cands):
                g["eff"] = {c: math.sqrt(g["d2"][c]) - aa[c] for c in cands}
            else:
                z = {c: zval[els[c]] - sum(e for _o, _r, e, _p in own[c])
                     - (PAIR_SCREEN if variant == "V3s" else 0.0) for c in cands}
                if all(v <= FLOAT_EPS for v in z.values()):
                    g["eff"] = {c: math.sqrt(g["d2"][c]) for c in cands}
                else:
                    g["eff"] = {c: (math.sqrt(g["d2"][c] / z[c]) if z[c] > FLOAT_EPS
                                    else math.inf) for c in cands}
            continue
        g["eff"] = {c: _eff(variant, math.sqrt(g["d2"][c]), rsel[c], aa[c], flag)
                    for c in cands}

    # per-case pair pooling (V2, V2c)
    pair_owner, cycles = {}, 0
    if variant in ("V2", "V2c"):
        acc = {}
        for g in geo:
            if not g["contested"] or not g["partners"]:
                continue
            near = g["near"]
            seen = set()
            for b in g["partners"]:
                if els[b] in seen:
                    continue
                seen.add(els[b])
                x = (g["eff"][b] - g["eff"][near]) / (math.sqrt(g["d2"][b])
                                                      + math.sqrt(g["d2"][near]))
                key = tuple(sorted((els[near], els[b])))
                s = acc.setdefault(key, [0.0, 0.0])
                s[0] += g["e"] * (x if els[near] == key[0] else -x)
                s[1] += g["e"]
        mean = {k: (v[0] / v[1] if v[1] > FLOAT_EPS else 0.0) for k, v in acc.items()}
        for k, m in mean.items():
            if abs(m) > FLOAT_EPS:
                pair_owner[k] = k[0] if m > 0 else k[1]
        # transitivity: reverse the weakest pair of every directed cycle
        while True:
            beats = {}
            for k, w in pair_owner.items():
                beats.setdefault(w, set()).add(k[1] if w == k[0] else k[0])
            cyc = _find_cycle(beats)
            if cyc is None:
                break
            cycles += 1
            edges = [tuple(sorted((cyc[i], cyc[(i + 1) % len(cyc)])))
                     for i in range(len(cyc))]
            weak = min(edges, key=lambda k: (abs(mean[k]), k))
            pair_owner[weak] = weak[1] if pair_owner[weak] == weak[0] else weak[0]
            mean[weak] = math.copysign(math.inf, mean[weak])   # never flipped twice

    for g in geo:
        vecs, d2, near = g["vecs"], g["d2"], g["near"]
        if not g["contested"]:
            g["hosts"] = {near: 1.0}
            continue
        cands = [near] + g["partners"]
        top = min(cands, key=lambda c: (g["eff"][c], d2[c], c))
        if owner_hook is not None:
            top = owner_hook(g, cands, els, top)
        elif pair_owner and g["partners"]:
            cel = sorted({els[c] for c in cands})
            win = [e for e in cel
                   if all(pair_owner.get(tuple(sorted((e, f)))) == e
                          for f in cel if f != e)]
            if len(win) == 1:
                top = min([c for c in cands if els[c] == win[0]],
                          key=lambda c: (g["eff"][c], d2[c], c))
        g["top"] = top
        order = sorted(range(na), key=lambda a: (d2[a], a))
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
            if tot > 0.0:
                g["hosts"] = {a: f / tot for a, f in force.items()}
            else:
                g["hosts"] = {top: 1.0}
        else:
            g["hosts"] = {top: 1.0}

    cloud = [0.0] * na
    for g in geo:
        for a in range(na):
            r = math.sqrt(g["d2"][a])
            if r <= REACH_SIGMA * g["rad"]:
                cloud[a] = max(cloud[a], r + REACH_SIGMA * g["rad"])
    for g in geo:
        g["held"] = any(g["d2"][a] <= cloud[a] * cloud[a] for a in range(na))
    return {"geo": geo, "has_occ": has_occ, "n_pair": n_pair, "na": na, "els": els,
            "zval": zval, "st": st, "r": rr, "a": aa, "rh": rh, "pair_owner": pair_owner,
            "cycles": cycles}


def _find_cycle(beats):
    """A directed cycle (list of nodes) in the graph node -> set of nodes, or None."""
    state, stack = {}, []

    def visit(u):
        state[u] = 1
        stack.append(u)
        for v in sorted(beats.get(u, ())):
            if state.get(v, 0) == 1:
                return stack[stack.index(v):]
            if state.get(v, 0) == 0:
                c = visit(v)
                if c:
                    return c
        stack.pop()
        state[u] = 2
        return None

    for u in sorted(beats):
        if state.get(u, 0) == 0:
            c = visit(u)
            if c:
                return c
    return None


# ---------------------------------------------------------------- surplus
def _classes(an):
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
    out = list(lost)
    for cls in classes:
        if len(cls) > 1:
            m = sum(lost[a] for a in cls) / len(cls)
            for a in cls:
                out[a] = m
    return out


def _empty_wfs(wfs, amount, lost):
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
    out = []
    for g in sorted(wfs, key=key):
        if out:
            k0, k1 = key(out[-1][0]), key(g)
            if k0[:-1] == k1[:-1] and abs(k0[-1] - k1[-1]) <= FLOAT_EPS:
                out[-1].append(g)
                continue
        out.append([g])
    return out


def _shape_removal(an, ne, amount):
    geo, na, els, zval = an["geo"], an["na"], an["els"], an["zval"]
    lost = [0.0] * na
    q = [zval[els[a]] - ne[a] for a in range(na)]
    donors = sorted({g["near"] for g in geo if g["own"] and g["pp"] and g["left"] > FLOAT_EPS})
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
        rest = [g for g in geo if g["left"] > FLOAT_EPS]
        for tier in _tiers(rest, lambda g: (-g["spread"],)):
            amount = _empty_wfs(tier, amount, lost)
            if amount <= FLOAT_EPS:
                break
    return lost


def _populations(case, variant=None, owner_hook=None):
    an = _analyse(case, variant, owner_hook)
    geo, na, els, zval, n_pair = an["geo"], an["na"], an["els"], an["zval"], an["n_pair"]
    ne = [0.0] * na
    for g in geo:
        for a, f in g["hosts"].items():
            ne[a] += n_pair * f
        g["left"] = float(n_pair)
    an["mode"] = "closed"
    if case.eoscar.nelect_wann is None:
        return ne, an
    excess = n_pair * len(geo) - int(round(case.eoscar.nelect_wann))
    an["excess"] = excess
    if excess <= 0:
        return ne, an
    lost = [0.0] * na
    if an["has_occ"]:
        an["mode"] = "occupations"
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
            _empty_wfs(lev, excess, lost)
            excess = 0
        for a in range(na):
            ne[a] -= lost[a]
        return ne, an
    an["mode"] = "no occupations"
    unheld = [g for g in geo if not g["held"]]
    if unheld:
        excess = _empty_wfs(unheld, excess, lost)
    ne0 = [ne[a] - lost[a] for a in range(na)]
    classes = _classes(an)
    if excess > FLOAT_EPS:
        x = _pool(_shape_removal(an, ne0, excess), classes)
        for a in range(na):
            lost[a] += x[a]
    for a in range(na):
        ne[a] -= lost[a]
    return ne, an


def _fragments(an):
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


def _apportion(ne, frac, frag=None):
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
            for a in pick_spread(g, rest):
                out[a] += 1
            rest = 0
    return out


def _run(case, variant=None, owner_hook=None):
    ne, an = _populations(case, variant, owner_hook)
    frag = _fragments(an)
    n = _apportion(ne, an["st"].frac, frag)
    zval, els = an["zval"], an["els"]
    os_ = [int(round(zval[e])) - k for e, k in zip(els, n)]
    supply = an["n_pair"] * len(an["geo"])
    placed = supply if case.eoscar.nelect_wann is None else \
        min(supply, int(round(case.eoscar.nelect_wann)))
    target = int(round(sum(zval[e] for e in els))) - placed
    assert sum(os_) == target, (sum(os_), target)
    an["populations"] = ne
    return os_, an


def assign(case):
    return _run(case)[0]
