"""Bond-graph IUPAC counting under the NELECT cap (with_table), round 2.

Every Wannier function (WF) is attached to a site of the structure's own
bond graph, and all WFs - shared ones included - then go through ONE
capacity-constrained min-cost flow with the shipped cost and bounds, so that
every WF can be left empty by the same cost criterion and the electrons add
up to NELECT.

 1. Bond graph, no distance cutoff. Two atom images are bonded when their
    diametral sphere holds no other atom (Gabriel graph of the structure; an
    atom exactly on the sphere does not break the bond).
 2. Attachment. A WF belongs to the nearest site. Sites: an atom; the
    midpoint of a bonded pair; the centroid of three or more mutually bonded
    atoms of one element (face, interstice). A pair or larger site can win
    only when its centre is nearer to the WF than the nearest atom - on a
    bond axis the boundary lies at 1/4 of the bond length from each nucleus
    (the 1/4 attachment boundary). Sites at exactly equal distance merge.
 3. Ownership. A WF on an atom or on a heteronuclear pair is unshared; it
    counts as owned by the atom of its site to which it is cheapest (shipped
    cost): the atom itself, or normally the more electronegative end of the
    pair. This replaces the r1 "hetero_minus" bookkeeping: a heteronuclear
    bond WF is counted in one valence shell only, that of its owner.
 4. Shares. A WF on a homonuclear site is split equally between the site's
    atoms (IUPAC equal split; each atom image counts once).
    Valence-shell rule: per spin channel an atom holds at most n_max/2 WFs,
    owned and shared WFs counted alike (n_max = ZVAL - OS_MIN of the shipped
    capacity bounds: the octet for main-group atoms, the shipped wide bound
    for transition metals). An excess is taken from the atom's homonuclear
    MULTIPLE bonds only (pairs joined by two or more WFs of one spin): that
    fraction of its shares in them goes to the partners. Single bonds are
    never converted, so hypervalent sigma chains keep the equal split.
 5. Crossing pairs. A same-element pair of another element whose midpoint is
    also nearer to the WF than any atom and which is not bonded only because
    atoms of the winning site sit inside its sphere (the two diagonals of an
    M2X2 rhombus), or two bonded pairs at exactly equal distance, are both
    offered to the flow; inside the pair that receives the WF it is split
    equally.
 6. Occupations (only when the EOSCAR lists one per WF). The depletion of a
    WF (full minus occupation) is booked on its owner(s). An atom's
    depletion sum is handed on to atoms that are linked to it by a WF, are
    less electronegative and carry more depletion than it does, in
    proportion to their depletion (ionic approximation: the hole of a shared
    orbital belongs to the electropositive partner). The sums are rounded to
    whole holes by largest remainder and every atom's holes are taken out of
    its own most depleted WFs before the flow. Two atoms are linked when
    they are the two atoms nearest to some WF or share a shared WF's site.
 7. Flow. Unshared WFs may go to any atom, shared WFs only to the atoms of
    their site and only up to each atom's share, at the mean cost of that
    split; cost d^2 + beta (1 - chi/chi_max), capacities from the shipped
    bounds, total capped at NELECT. What a shared WF ships is divided by its
    shares. Electrons below the cap that a site cannot hold (its atoms are
    full) are shipped like those of an unshared WF.
 8. Cap symmetry (no occupations). Atoms with one WF signature (element; per
    spin, as a multiset over the spins, the number of owned WFs and the
    shared sites) are pooled when the class hosts a WF that the cap left
    (partly) empty: which of several equivalent WFs stays empty is not
    determined by the input.
 9. Group symmetry. In a group of atoms connected by shared WFs, when some
    atom's count is fractional, atoms that the WF bond graph cannot tell
    apart (same element, same elements among the atoms linked to them) get
    equal counts.
10. Integers. Fractions exist only inside such groups; each group's total is
    apportioned by largest remainder, atoms of equal count being kept equal
    whenever whole blocks of them can absorb the remainder; otherwise equal
    remainders are served fewer electrons first, then the atom nearer to its
    shared WFs, then position.
"""
from __future__ import annotations

import math
from fractions import Fraction

from pyleos.capacity import eneg
from pyleos.parsers import inv3
from pyleos.transport import BIG_NEG, MinCostFlow, _build_costs

BETA = 10.0
TOL = 1e-6      # the one numerical tolerance (relative; cost quantum in A^2)

META = {
    "name": "bond-graph IUPAC counting under the NELECT cap",
    "mechanism": (
        "Gabriel (empty-diametral-sphere) bond graph of the structure, no "
        "distance cutoff. Each WF attaches to the nearest site: an atom, a "
        "bonded pair's midpoint, or the centroid of three or more mutually "
        "bonded atoms of one element (face, interstice). A shared site wins "
        "only when its centre is nearer than the nearest atom, i.e. beyond "
        "1/4 of the bond length from a nucleus (declared rule: the 1/4 "
        "attachment boundary); sites at exactly equal distance merge. "
        "Declared rule (ownership, replaces r1 hetero_minus): a WF on an "
        "atom or a heteronuclear pair is owned by the atom of its site to "
        "which the shipped cost is lowest, and is counted in that atom's "
        "valence shell only. WFs on homonuclear sites are split equally "
        "between the site's atoms (IUPAC). Declared rule (valence shell): "
        "per spin an atom holds at most n_max/2 WFs, owned and shared alike "
        "(n_max from the shipped capacity bounds: the octet for main-group "
        "atoms, the shipped wide bound for transition metals); an excess is "
        "taken from its homonuclear multiple bonds only (two or more WFs of "
        "one spin on the pair) and goes to the partners; single bonds are "
        "never converted. Declared rule (crossing pairs): a pair of another "
        "element that is unbonded only because the winning pair sits inside "
        "its sphere is offered to the flow as an alternative, split equally "
        "if chosen. All WFs go through one min-cost flow with the shipped "
        "cost and bounds under the NELECT cap; shared WFs are restricted to "
        "their site and shares and can be left empty like any other WF; "
        "electrons a full site cannot hold are shipped as unshared. With one "
        "occupation per WF the holes are placed before the flow: per-atom "
        "depletion sums, handed on to WF-linked, less electronegative atoms "
        "that carry more depletion, rounded by largest remainder and taken "
        "from the atom's most depleted WFs (two atoms are WF-linked when "
        "they are the two atoms nearest to some WF or share a shared WF's "
        "site). Without occupations, "
        "atoms with one WF signature (spin channels as a multiset) share "
        "equally when their class hosts a WF the cap left empty. Declared "
        "rule (group symmetry): in a group of atoms connected by shared WFs "
        "that has fractional counts, atoms with the same element and the "
        "same elements among their WF-linked atoms get equal counts. "
        "Fractions are rounded inside each group by largest remainder; atoms "
        "of equal count stay equal whenever whole blocks of them can absorb "
        "the remainder, otherwise ties go to fewer electrons first, then "
        "nearer to its shared WFs, then position."),
    "params": {
        "beta": [BETA, "shipped transport electronegativity weight (plateau "
                       "[8,12]); unchanged from the baseline"],
        "TOL": [TOL, "numerical tolerance, not a physical scale: relative "
                     "slack of every geometric comparison (empty-sphere "
                     "test, nearest-site ties, candidate ball, second-"
                     "nearest atom) and the quantum (A^2) to which flow "
                     "costs and tie-break keys are rounded so that the "
                     "solver works on exact integers; only exact "
                     "degeneracies are affected"],
    },
    "table_free": False,
}


# ---------------------------------------------------------------- geometry
def _d2(p, q):
    return (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2


class _Geometry:
    def __init__(self, st):
        latt = st.lattice
        self.latt = latt
        self.binv = inv3(latt)
        self.frac = st.frac
        self.cart = [[f[0] * latt[0][j] + f[1] * latt[1][j] + f[2] * latt[2][j]
                      for j in range(3)] for f in st.frac]
        # spacing of the lattice planes along each fractional coordinate
        self.height = [1.0 / math.sqrt(sum(self.binv[i][k] ** 2 for i in range(3)))
                       for k in range(3)]
        self._blk = {}

    def images(self, xc, rmax):
        """(d2, atom, cell, cart) of every atom image within rmax of xc.
        An image n cells away from the nearest-fraction image is at least
        (|n| - 1/2) plane spacings off, which bounds the cells to visit."""
        binv, latt = self.binv, self.latt
        xf = [xc[0] * binv[0][j] + xc[1] * binv[1][j] + xc[2] * binv[2][j]
              for j in range(3)]
        m0, m1, m2 = (int(math.floor(rmax / h + Fraction(1, 2)))
                      for h in self.height)
        r2 = rmax * rmax
        out = []
        for a, fa in enumerate(self.frac):
            b0 = round(xf[0] - fa[0])
            b1 = round(xf[1] - fa[1])
            b2 = round(xf[2] - fa[2])
            ca = self.cart[a]
            for n0 in range(b0 - m0, b0 + m0 + 1):
                for n1 in range(b1 - m1, b1 + m1 + 1):
                    for n2 in range(b2 - m2, b2 + m2 + 1):
                        x = ca[0] + n0 * latt[0][0] + n1 * latt[1][0] + n2 * latt[2][0]
                        y = ca[1] + n0 * latt[0][1] + n1 * latt[1][1] + n2 * latt[2][1]
                        z = ca[2] + n0 * latt[0][2] + n1 * latt[1][2] + n2 * latt[2][2]
                        d2 = (x - xc[0]) ** 2 + (y - xc[1]) ** 2 + (z - xc[2]) ** 2
                        if d2 <= r2:
                            out.append((d2, a, (n0, n1, n2), (x, y, z)))
        return out

    def blockers(self, p, q):
        """Atom images strictly inside the diametral sphere of the images
        p, q = (atom, cell, cart): frozenset of (atom, cell). Empty for a
        bonded (Gabriel) pair. Cached per periodic pair."""
        a, ca, pa = p
        b, cb, pb = q
        rel = (cb[0] - ca[0], cb[1] - ca[1], cb[2] - ca[2])
        got = self._blk.get((a, b, rel))
        if got is None:
            m = [(pa[j] + pb[j]) / 2 for j in range(3)]
            r2 = _d2(pa, pb) / 4
            got = frozenset(
                (x, (c[0] - ca[0], c[1] - ca[1], c[2] - ca[2]))
                for d2, x, c, _ in self.images(m, math.sqrt(r2))
                if d2 < r2 * (1 - TOL))
            self._blk[(a, b, rel)] = got
            self._blk[(b, a, (-rel[0], -rel[1], -rel[2]))] = frozenset(
                (x, (c[0] - rel[0], c[1] - rel[1], c[2] - rel[2]))
                for x, c in got)
        return frozenset((x, (c[0] + ca[0], c[1] + ca[1], c[2] + ca[2]))
                         for x, c in got)


def _edge_key(p, q):
    """Periodic identity of an image pair (independent of the cell both
    images are printed in)."""
    a, ca = p[0], p[1]
    b, cb = q[0], q[1]
    rel = (cb[0] - ca[0], cb[1] - ca[1], cb[2] - ca[2])
    if a > b or (a == b and rel < (0, 0, 0)):
        a, b = b, a
        rel = (-rel[0], -rel[1], -rel[2])
    return (a, b, rel)


def _cliques(nodes, adj):
    """All sets of three or more mutually adjacent nodes."""
    out = []

    def grow(cur, cand):
        if len(cur) >= 3:
            out.append(tuple(cur))
        for k, v in enumerate(cand):
            grow(cur + [v], [u for u in cand[k + 1:] if u in adj[v]])

    grow([], sorted(nodes))
    return out


def _attach_one(geo, el, centre, d1_hint):
    """Attachment of one WF centre.
      ("v", atom)                     own WF of the atom
      ("x", (atom, ...))              heteronuclear site: unshared in the flow
      ("u", [option, ...])            homonuclear site(s); option = tuple of
                                      (atom, cell) images, first = nearest
    A shared site can only win when its centre is nearer than the nearest
    atom (distance d1); with the empty-sphere condition this forces every
    atom of it to lie within 3*d1 of the WF, so that ball holds every
    possible winner (no enumeration cutoff)."""
    imgs = geo.images(centre, 3 * d1_hint * (1 + TOL) + TOL)
    imgs.sort(key=lambda t: t[0])
    d1sq = imgs[0][0]
    n = len(imgs)
    dist = [math.sqrt(t[0]) for t in imgs]
    pid = [(t[1], t[2], t[3]) for t in imgs]
    sites = []                      # (d2 of the centre, members)
    bonded_same = {}                # adjacency of same-element bonded pairs
    blocked = []                    # (dm2, i, j, first blocker) same element
    lim = d1sq * (1 - TOL)
    for i in range(n):
        p = imgs[i][3]
        for j in range(i + 1, n):
            q = imgs[j][3]
            mx, my, mz = (p[0] + q[0]) / 2, (p[1] + q[1]) / 2, (p[2] + q[2]) / 2
            dm2 = (mx - centre[0]) ** 2 + (my - centre[1]) ** 2 + (mz - centre[2]) ** 2
            if dm2 >= lim:
                continue
            # an atom inside the pair's sphere is within radius + |m - WF| of
            # the WF, hence in the ball; the nearest ones are tried first
            r2 = _d2(p, q) / 4
            reach = math.sqrt(r2) + math.sqrt(dm2)
            r2 *= 1 - TOL
            hit = None
            for k in range(n):
                if dist[k] > reach:
                    break
                z = imgs[k][3]
                if ((z[0] - mx) ** 2 + (z[1] - my) ** 2 + (z[2] - mz) ** 2 < r2
                        and k != i and k != j):
                    hit = k
                    break
            same = el[imgs[i][1]] == el[imgs[j][1]]
            if hit is None:
                sites.append((dm2, (i, j)))
                if same:
                    bonded_same.setdefault(i, set()).add(j)
                    bonded_same.setdefault(j, set()).add(i)
            elif same:
                blocked.append((dm2, i, j, hit))
    for cl in _cliques(list(bonded_same), bonded_same):
        cen = [sum(imgs[i][3][k] for i in cl) / len(cl) for k in range(3)]
        sites.append((_d2(cen, centre), cl))
    if not sites:
        near = sorted({t[1] for t in imgs if t[0] <= d1sq * (1 + TOL) + TOL * TOL})
        if len(near) == 1:
            return ("v", near[0])
        return ("x", tuple(near))
    dmin = min(s[0] for s in sites)
    tied = [s[1] for s in sites if s[0] <= dmin + TOL * d1sq]
    members = sorted({i for s in tied for i in s})
    if any(len({el[imgs[i][1]] for i in s}) > 1 for s in tied):
        return ("x", tuple(sorted({imgs[i][1] for i in members})))
    by_el = {}
    for i in members:
        by_el.setdefault(el[imgs[i][1]], []).append(i)
    # several elements can only tie exactly (two crossing bonded pairs with
    # one midpoint); the option whose atoms are nearer comes first
    ranked = sorted(by_el.items(), key=lambda kv: (
        round(sum(imgs[i][0] for i in kv[1]) / len(kv[1]) / (TOL * d1sq)), kv[0]))
    options = [tuple((imgs[i][1], imgs[i][2]) for i in grp) for _, grp in ranked]
    if len(options) == 1:
        # crossing pairs: pairs of another element that are not bonded only
        # because atoms of the winning site sit inside their sphere
        win = set(options[0])
        seen = set()
        for dm2, i, j, hit in sorted(blocked, key=lambda t: t[:3]):
            pair = ((imgs[i][1], imgs[i][2]), (imgs[j][1], imgs[j][2]))
            if (el[pair[0][0]] == ranked[0][0]
                    or (imgs[hit][1], imgs[hit][2]) not in win
                    or not geo.blockers(pid[i], pid[j]) <= win):
                continue
            key = _edge_key(pair[0], pair[1])
            if key not in seen:
                seen.add(key)
                options.append(pair)
    return ("u", options)


# ---------------------------------------------------------------- counting
def _links(el, cost_rows, att):
    """WF bond graph: two atoms are linked when they are the two atoms
    nearest to some WF (exact ties all count) or share the site of a shared
    WF. Local by construction: no WF has its two nearest atoms in different
    molecules. Distances are the minimum-image ones behind the shipped
    costs."""
    na = len(el)
    chi_max = max(eneg(e) for e in set(el)) or 1.0
    off = [BETA * (1.0 - eneg(e) / chi_max) for e in el]
    links = [set() for _ in range(na)]

    def join(atoms):
        for x in atoms:
            for y in atoms:
                if x != y:
                    links[x].add(y)

    for w, row in enumerate(cost_rows):
        if att[w][0] == "u":
            join(sorted({x for x, _ in att[w][1][0]}))
        if na < 2:
            continue
        d2 = sorted((row[a] - off[a], a) for a in range(na))
        k = 2                       # the two nearest, plus whatever ties with the second
        while k < na and d2[k][0] <= d2[1][0] * (1 + TOL) + TOL * TOL:
            k += 1
        if d2[1][0] <= d2[0][0] * (1 + TOL) + TOL * TOL:
            join([a for _, a in d2[:k]])
        else:
            for _, a in d2[1:k]:
                join([d2[0][1], a])
    return links


def _owners(att, cost_rows):
    """Unshared WFs: the atom of the WF's own site to which it is cheapest
    (an atom-centred WF: that atom; a heteronuclear bond WF: normally the
    more electronegative end), as dict atom -> share. None for shared WFs."""
    out = []
    for w, a in enumerate(att):
        if a[0] == "v":
            out.append({a[1]: Fraction(1)})
        elif a[0] == "x":
            k = [int(round(cost_rows[w][x] / TOL)) for x in a[1]]
            best = [x for x, v in zip(a[1], k) if v == min(k)]
            out.append({x: Fraction(1, len(best)) for x in best})
        else:
            out.append(None)
    return out


def _shares(case, wfs, att, owners, bounds):
    """Per WF: list of options, each a dict atom -> Fraction share (summing
    to 1), or None for an unshared WF. Also the graph signature per atom."""
    el = case.structure.elements
    na = len(el)
    nch = len(case.channels)
    half = Fraction(1, 2)
    touch = [[Fraction(0)] * na for _ in range(nch)]
    mult = [dict() for _ in range(nch)]       # same-element pair -> WFs on it
    own = [[Fraction(0)] * na for _ in range(nch)]
    part = [[[] for _ in range(na)] for _ in range(nch)]
    for w, a in enumerate(att):
        c = wfs[w].channel - 1
        if owners[w] is not None:
            for x, s in owners[w].items():
                touch[c][x] += s
                own[c][x] += s
            continue
        site = a[1][0]
        atoms = sorted({x for x, _ in site})
        for x in atoms:
            touch[c][x] += 1
            part[c][x].append((len(site), el[x]))
        if len(site) == 2:
            key = _edge_key(site[0], site[1])
            mult[c][key] = mult[c].get(key, 0) + 1
    # atoms over their valence shell give up shares in their multiple bonds
    n_multi = [[0] * na for _ in range(nch)]
    for c in range(nch):
        for (a, b, _), m in mult[c].items():
            if m >= 2 and a != b:
                n_multi[c][a] += m
                n_multi[c][b] += m
    red = [[Fraction(0)] * na for _ in range(nch)]
    for c in range(nch):
        for a in range(na):
            excess = touch[c][a] - Fraction(bounds[a][1], 2)
            if excess > 0 and n_multi[c][a]:
                red[c][a] = min(Fraction(1), excess / n_multi[c][a])
    shares = []
    for w, a in enumerate(att):
        if a[0] != "u":
            shares.append(None)
            continue
        c = wfs[w].channel - 1
        opts = []
        for k, site in enumerate(a[1]):
            sh = {}
            for x, _ in site:
                sh[x] = sh.get(x, 0) + Fraction(1, len(site))
            if k == 0 and len(site) == 2 and len(sh) == 2:
                p, q = site[0][0], site[1][0]
                if mult[c][_edge_key(site[0], site[1])] >= 2:
                    sp = half * (1 - red[c][p]) + half * red[c][q]
                    sh = {p: sp, q: 1 - sp}
            opts.append({x: s for x, s in sh.items() if s > 0})
        shares.append(opts)
    # spin labels carry no chemistry: the channels enter as a multiset
    sig = [(el[a], tuple(sorted((own[c][a], tuple(sorted(part[c][a])))
                                for c in range(nch)))) for a in range(na)]
    return shares, sig


def _occupations(case, wfs):
    """Depletion (full minus occupation, clipped at 0) of every WF, or None
    when the EOSCAR does not list exactly one occupation per WF."""
    det = case.eoscar.wf_detail
    if not det or len(det) != len(wfs):
        return None
    occ = {}
    for d in det:
        occ[(d.spin, d.index)] = d.occ
    if len(occ) != len(wfs):
        return None
    out = []
    for wf in wfs:
        o = occ.get((wf.channel, wf.index))
        if o is None:
            return None
        out.append(min(float(wf.supply), max(0.0, wf.supply - o)))
    return out


def _holes(case, links, wfs, owners, shares, depl, n_holes, trace=None):
    """Whole holes taken out of individual WFs before the flow:
    removed[w] = dict atom -> electrons (int) of WF w's part owned by atom."""
    el = case.structure.elements
    na = len(el)
    owners = [own if own is not None else shares[w][0]
              for w, own in enumerate(owners)]
    dsum = [0.0] * na
    for w, own in enumerate(owners):
        for x, s in own.items():
            dsum[x] += depl[w] * float(s)
    scale = max(dsum) if dsum else 0.0
    if n_holes <= 0 or scale <= 0.0:
        return [dict() for _ in wfs]
    # ionic approximation: hand an atom's depletion on to bonded, less
    # electronegative neighbours that carry more of it
    hole = [0.0] * na
    for x in range(na):
        if dsum[x] <= TOL * scale:
            continue
        chi = eneg(el[x])
        wts = {}
        for m in links[x]:
            if dsum[m] > dsum[x] * (1 + TOL) and eneg(el[m]) < chi:
                wts[m] = dsum[m]
        if wts:
            tot = sum(wts.values())
            for m, v in wts.items():
                hole[m] += dsum[x] * v / tot
        else:
            hole[x] += dsum[x]
    tot = sum(hole)
    budget = min(n_holes, int(math.floor(tot + Fraction(1, 2))))
    quota = [h * budget / tot for h in hole]
    count = [int(math.floor(q + TOL)) for q in quota]
    left = max(0, budget - sum(count))
    order = sorted(range(na), key=lambda a: -(quota[a] - count[a]))
    rem = [quota[a] - count[a] for a in order]
    # a cut through equal remainders is not decided here: those holes are
    # left to the flow cap and the class pooling
    while 0 < left < na and rem[left - 1] - rem[left] <= TOL:
        left -= 1
    for a in order[:left]:
        count[a] += 1
    if trace is not None:
        trace.update(depletion=dsum, holes=hole, hole_count=count)
    # take each atom's holes out of its own most depleted WFs
    removed = [dict() for _ in wfs]
    mine = [[] for _ in range(na)]
    for w, own in enumerate(owners):
        for x, s in own.items():
            if depl[w] > 0.0:
                mine[x].append((-depl[w] * float(s), w, s))
    for a in range(na):
        need = count[a]
        for _, w, s in sorted(mine[a], key=lambda t: (t[0], t[1])):
            if need <= 0:
                break
            room = int(math.floor(wfs[w].supply * s))
            take = min(need, room)
            if take > 0:
                removed[w][a] = take
                need -= take
    return removed


def electron_counts(case, trace=None):
    """Fractional electron count per atom (Fractions) plus what the rounding
    needs: groups of atoms that share electrons and a tie-break key."""
    st = case.structure
    el = st.elements
    wfs, cost_rows, bounds, target, na, nw = _build_costs(
        case, BETA, "allen", False, 0.0)
    geo = _Geometry(st)
    att = [_attach_one(geo, el, centre, wfs[k0 + p].d1)
           for k0, ch in _channel_offsets(case) for p, centre in enumerate(ch.centres)]
    owners = _owners(att, cost_rows)
    links = _links(el, cost_rows, att)
    shares, sig = _shares(case, wfs, att, owners, bounds)
    total = sum(wf.supply for wf in wfs)
    if target is None:
        target = total

    depl = _occupations(case, wfs)
    removed = [dict() for _ in wfs]
    if depl is not None:
        removed = _holes(case, links, wfs, owners, shares, depl,
                         total - target, trace)

    # ---- arcs of the flow, in electrons (Fractions) and shipped costs
    arcs = []                      # per WF: list of (atom, cap, cost, option)
    supply = []
    for w, wf in enumerate(wfs):
        s = wf.supply
        if shares[w] is None:
            s_left = s - sum(removed[w].values())
            arcs.append([(a, Fraction(s_left), cost_rows[w][a], None)
                         for a in range(na)])
            supply.append(Fraction(s_left))
            continue
        lst = []
        for k, opt in enumerate(shares[w]):
            cost = sum(float(sh) * cost_rows[w][a] for a, sh in opt.items())
            for a, sh in sorted(opt.items()):
                cap = s * sh - (removed[w].get(a, 0) if k == 0 else 0)
                if cap > 0:
                    lst.append((a, cap, cost, k))
        arcs.append(lst)
        supply.append(Fraction(s - sum(removed[w].values())))

    den = 1
    for lst in arcs:
        for _, cap, _, _ in lst:
            den = den * cap.denominator // math.gcd(den, cap.denominator)
    quantum = TOL
    big = int(round(BIG_NEG / quantum))
    S, T = 0, nw + na + 1
    mcf = MinCostFlow(nw + na + 2)
    for w in range(nw):
        mcf.add(S, 1 + w, int(supply[w] * den), 0)
        for a, cap, cost, _ in arcs[w]:
            mcf.add(1 + w, 1 + nw + a, int(cap * den), int(round(cost / quantum)))
    for a in range(na):
        n_min, n_max = bounds[a]
        if n_min > 0:
            mcf.add(1 + nw + a, T, n_min * den, -big)
        if n_max > n_min:
            mcf.add(1 + nw + a, T, (n_max - n_min) * den, 0)
    goal = min(target, int(sum(supply)))
    mcf.solve(S, T, max_flow=goal * den)

    # ---- conservation: electrons below the cap that their site could not
    # hold (its atoms are full) are shipped like those of an unshared WF
    sent_all = [[int(cap * den) - mcf.g[1 + w][1 + i][1]
                 for i, (_, cap, _, _) in enumerate(arcs[w])] for w in range(nw)]
    short = goal * den - sum(sum(x) for x in sent_all)
    spill = [0] * na
    if short > 0:
        inflow = [0] * na
        for w in range(nw):
            for (a, _, _, _), f in zip(arcs[w], sent_all[w]):
                inflow[a] += f
        rest = [w for w in range(nw) if int(supply[w] * den) > sum(sent_all[w])]
        mcf2 = MinCostFlow(len(rest) + na + 2)
        T2 = len(rest) + na + 1
        for k, w in enumerate(rest):
            cap = int(supply[w] * den) - sum(sent_all[w])
            mcf2.add(0, 1 + k, cap, 0)
            for a in range(na):
                mcf2.add(1 + k, 1 + len(rest) + a, cap,
                         int(round(cost_rows[w][a] / quantum)))
        for a in range(na):
            if bounds[a][1] * den > inflow[a]:
                mcf2.add(1 + len(rest) + a, T2, bounds[a][1] * den - inflow[a], 0)
        mcf2.solve(0, T2, max_flow=short)
        for k, w in enumerate(rest):
            for a in range(na):
                spill[a] += int(supply[w] * den) - sum(sent_all[w]) - mcf2.g[1 + k][1 + a][1]

    ne = [Fraction(x, den) for x in spill]
    near = [0.0] * na              # tie-break: cost of the atom's shared electrons
    groups = []                    # atom groups that share electrons
    ghost_hosts = set()
    for w in range(nw):
        sent = [Fraction(x, den) for x in sent_all[w]]
        shipped = sum(sent)
        if shipped < supply[w]:
            ghost_hosts.update(owners[w] if owners[w] is not None
                               else shares[w][0])
        if shares[w] is None:
            for (a, _, _, _), f in zip(arcs[w], sent):
                if f:
                    ne[a] += f
            continue
        for k in range(len(shares[w])):
            idx = [i for i, arc in enumerate(arcs[w]) if arc[3] == k]
            got = sum(sent[i] for i in idx)
            if not got:
                continue
            cap_k = sum(arcs[w][i][1] for i in idx)
            for i in idx:
                a = arcs[w][i][0]
                part = got * arcs[w][i][1] / cap_k
                ne[a] += part
                near[a] += float(part) * cost_rows[w][a]
            if len(idx) > 1:
                groups.append([arcs[w][i][0] for i in idx])
    if trace is not None:
        trace.update(att=att, owners=owners, shares=shares, removed=removed,
                     goal=goal, shipped=sum(ne), supply=supply, arcs=arcs,
                     spill=Fraction(sum(spill), den),
                     sent=[[Fraction(x, den) for x in row] for row in sent_all])
    return ne, groups, near, sig, ghost_hosts, links


def _channel_offsets(case):
    out, k = [], 0
    for ch in case.channels:
        out.append((k, ch))
        k += len(ch.centres)
    return out


def pool_classes(ne, sig, ghost_hosts):
    """Cap symmetry: atoms with one graph signature share their electrons
    equally when the class hosts a WF that the cap left (partly) empty."""
    classes = {}
    for a, s in enumerate(sig):
        classes.setdefault(s, []).append(a)
    out = list(ne)
    pooled = []
    for cls in classes.values():
        if len(cls) < 2 or not ghost_hosts.intersection(cls):
            continue
        mean = sum(ne[a] for a in cls) / len(cls)
        for a in cls:
            out[a] = mean
        pooled.append(cls)
    return out, pooled


def _components(na, groups):
    """Sets of atoms connected through shared WFs or a pooled class."""
    root = list(range(na))

    def find(i):
        while root[i] != i:
            root[i] = root[root[i]]
            i = root[i]
        return i

    for g in groups:
        for a in g[1:]:
            root[find(a)] = find(g[0])
    comp = {}
    for a in range(na):
        comp.setdefault(find(a), []).append(a)
    return [c for c in comp.values() if len(c) > 1]


def pool_equivalent(case, links, ne, comps):
    """Group symmetry: when the WFs shared inside a group of atoms leave
    some of them with a fractional count (the WFs are not all whole pairs
    between two atoms, so the group total has to be apportioned), atoms of
    the group that the structure's bond graph cannot tell apart - same
    element, same elements among their bonded neighbours - share equally."""
    el = case.structure.elements
    out = list(ne)
    for atoms in comps:
        if all(ne[a].denominator == 1 for a in atoms):
            continue
        classes = {}
        for a in atoms:
            key = (el[a], tuple(sorted(el[b] for b in links[a])))
            classes.setdefault(key, []).append(a)
        for cls in classes.values():
            mean = sum(ne[a] for a in cls) / len(cls)
            for a in cls:
                out[a] = mean
    return out


def integer_counts(case, ne, comps, near):
    """Whole electrons per atom. Fractions exist only inside groups of atoms
    that share WFs (or a pooled class); each group's total is apportioned by
    largest remainder, with atoms of equal count kept equal: the remainder
    goes to whole blocks of equal atoms when their sizes can add up to it
    (the choice with the largest total remainder; fewer electrons first).
    Only when no such choice exists is one block split: larger remainder,
    then fewer electrons, then the atom nearer to its shared WFs, then
    position."""
    st = case.structure
    el = st.elements
    na = len(el)
    out = [int(x) if x.denominator == 1 else None for x in ne]
    pos = [tuple(round((f - math.floor(f)) / TOL) for f in st.frac[a])
           for a in range(na)]

    def apportion(atoms, total):
        fl = {a: int(math.floor(ne[a])) for a in atoms}
        left = total - sum(fl.values())
        for a in atoms:
            out[a] = fl[a]
        if left <= 0:
            return
        blocks = {}
        for a in atoms:
            if ne[a] != fl[a]:
                blocks.setdefault(ne[a], []).append(a)
        order = sorted(blocks, key=lambda v: (-(v - math.floor(v)), v))
        best = {0: (Fraction(0), ())}       # electrons placed -> (remainder, blocks)
        for k, v in enumerate(order):
            size = len(blocks[v])
            gain = size * (v - math.floor(v))
            for placed, (val, pick) in list(best.items()):
                if placed + size > left:
                    continue
                new = (val + gain, pick + (k,))
                cur = best.get(placed + size)
                if cur is None or new[0] > cur[0] or (new[0] == cur[0]
                                                     and new[1] < cur[1]):
                    best[placed + size] = new
        if left in best:
            for k in best[left][1]:
                for a in blocks[order[k]]:
                    out[a] += 1
            return
        ranked = sorted(atoms, key=lambda a: (-(ne[a] - fl[a]), fl[a],
                                              round(near[a] / TOL), pos[a], a))
        for a in ranked[:left]:
            out[a] += 1

    loose = {}
    seen = set()
    for atoms in comps:
        seen.update(atoms)
        if all(ne[a].denominator == 1 for a in atoms):
            continue
        tot = sum(ne[a] for a in atoms)
        if tot.denominator == 1:
            apportion(atoms, int(tot))
        else:
            loose.setdefault(el[atoms[0]], []).extend(atoms)
    for a in range(na):
        if a not in seen and ne[a].denominator != 1:
            loose.setdefault(el[a], []).append(a)
    rest = []
    for atoms in loose.values():
        tot = sum(ne[a] for a in atoms)
        if tot.denominator == 1:
            apportion(atoms, int(tot))
        else:
            rest.extend(atoms)
    if rest:
        tot = sum(ne[a] for a in rest)
        apportion(rest, int(math.floor(tot + Fraction(1, 2))))
    return out


def assign(case):
    el = case.structure.elements
    zval = case.eoscar.zval
    ne, groups, near, sig, ghost_hosts, links = electron_counts(case)
    ne, pooled = pool_classes(ne, sig, ghost_hosts)
    comps = _components(len(el), groups + pooled)
    ne = pool_equivalent(case, links, ne, comps)
    n_int = integer_counts(case, ne, comps, near)
    return [round(zval[el[a]]) - n_int[a] for a in range(len(el))]
