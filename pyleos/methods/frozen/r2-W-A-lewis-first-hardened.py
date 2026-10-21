"""Lewis-first oxidation states, hardened cap path (with_table).

Every Wannier function (WF) is classified from geometry, electrons are
counted with the IUPAC ionic-approximation rules, and surplus WFs (WF supply
above NELECT) are removed by explicit rules before the Lewis structure is
built. A min-cost flow is used only when a capacity bound is broken.

1. Host set (geometry only). A WF starts on its nearest atom A (lone pair /
   core). It is a two-centre bond A-B when its centre is closer to the A-B
   midpoint than to A, and a same-element set keeps absorbing atoms while the
   centroid of the set keeps approaching the centre (multi-centre bond; a
   last atom of another element may join the same way and makes the set
   heteronuclear). Length scale from the structure itself: an image may
   join only if the WF lies inside the bond sphere of the enlarged set
   (distance to the centroid below half the largest member-member
   distance), no other atom lies inside that sphere, and the WF is closer
   to every member than any other atom is to it; a set never holds two
   images of one atom. Any other WF is a lone pair of its nearest atom.
2. IUPAC counting. Lone pair -> its atom; heteronuclear bond -> the partner
   with the higher Allen electronegativity; same-element (multi-centre) bond
   -> equal shares to all members.
3. Surplus WFs (supply above NELECT), decided before the Lewis structure:
   a. occupations available (EOSCAR #WF block covering every WF): a WF whose
      occupation rounds to zero electrons is absent. The remaining holes are
      placed on atoms, not on WFs: per-atom depletion sums D_a = sum over
      the atom's own WFs of (supply - occupation) are scaled to the hole
      count and rounded by largest remainder. Ionic approximation for empty
      states: only atoms that are the less electronegative partner of their
      closest heteronuclear contact (cation sites) receive holes while any
      of them can give up electrons (n_min of the capacity bounds).
   b. no occupations: electrons are removed one at a time, never below
      n_min: first from WFs counted in a valence shell that stays over
      n_max after the octet rule of step 4 (what the shell cannot hold is
      the surplus), then by the highest own
      Lewis-destination cost (shipped transport cost
      d^2 + beta (1 - chi / chi_max)), and only last from the one WF left
      in a bond per spin channel (the Lewis skeleton, kept unless every
      owner of the bond is over its shell). Atoms of one Lewis
      class (same element, same multiset of own WF kinds in their spin
      channels) then share the removed electrons equally (class mean, largest
      remainder), so equivalent atoms lose equally whatever the order of
      nearly tied costs.
4. Octet rule on the WFs that are left (shell capacity n_max of the shipped
   capacity bounds): an over-full atom turns its excess same-element pairs
   into lone pairs of a partner that is within its own shell, most polarised
   pair first, spread over bonds and spin channels, never a bond's last WF.
5. Integers. Equal shares are pooled per connected same-element fragment
   and rounded by largest remainder. A fragment that holds a multi-centre
   bond WF is a delocalised cluster whose members share its total equally.
6. Capacity. If a count breaks a capacity bound (or a remainder tie is
   left), one min-cost flow with integer costs decides: Lewis destinations
   at the shipped transport cost, the atoms around each WF as penalised
   alternatives, one capacity arc per atom.
"""
from __future__ import annotations

import math

from pyleos.capacity import capacity_bounds, eneg
from pyleos.parsers import inv3
from pyleos.transport import MinCostFlow, _min_image_d2

BETA = 10.0             # shipped transport weight
LEWIS_PENALTY = 100.0   # lexicographic: leaving a Lewis destination
SEMICORE_REWARD = 1.0e3  # shipped BIG_NEG: n_min shells filled first
SHELL_TOL = 1.0e-4      # Angstrom: numerical resolution of the method
MULTI_CENTRE = 3        # centres from which a bond WF counts as delocalised

_REVERSE_TIES = False   # audit switch: order in which exactly tied surplus
                        # WFs are emptied (index order or its reverse)

META = {
    "name": "lewis-first-hardened",
    "mechanism": (
        "Classify each WF from geometry as lone pair, 2-centre bond (centre "
        "closer to the bond midpoint than to the nearest atom) or "
        "multi-centre bond (centroid of the set keeps approaching the "
        "centre); a set is a bond only if the WF lies inside its bond sphere "
        "(diameter = largest member-member distance), no other atom does, "
        "and the WF is closer to every member than any other atom is to it. "
        "Count with IUPAC rules (heteronuclear -> higher Allen chi, same "
        "element -> equal shares). "
        "Surplus WFs are settled before the Lewis structure: with "
        "occupations, WFs that round to zero electrons are absent and the "
        "remaining holes go to cation-site atoms in proportion to their "
        "per-atom depletion sums (largest remainder); without occupations "
        "electrons are removed first from over-full valence shells, then "
        "by cost, a bond's last WF per spin channel last, and atoms of one "
        "Lewis class share the loss equally. "
        "Then octet repair (capacity_bounds n_max), per-fragment "
        "largest-remainder rounding, equal sharing inside delocalised "
        "clusters, and a min-cost flow with integer costs only when a "
        "capacity bound is broken."),
    "params": {
        "beta": [BETA, "shipped transport electronegativity weight; ranks "
                 "which surplus electrons are removed when the input carries "
                 "no occupations and prices the capacity flow"],
        "lewis_penalty": [LEWIS_PENALTY, "lexicographic priority, not tuned: "
                          "must exceed any d^2 + beta cost difference "
                          "between plausible hosts (~20 A^2) and stay below "
                          "the semicore reward"],
        "semicore_reward": [SEMICORE_REWARD, "shipped BIG_NEG of "
                            "pyleos.transport (n_min shells filled first)"],
        "shell_tol": [SHELL_TOL, "numerical resolution (A): distances equal "
                      "to this precision are one shell; also the quantum to "
                      "which flow costs are rounded (integer costs, exact "
                      "ties) and the guard of every float comparison; "
                      "wannier90 prints centres to 1e-6 A, real symmetry "
                      "partners differ by up to 1e-3 A and are handled by "
                      "the class rule, not by this number"],
        "multi_centre": [MULTI_CENTRE, "definition, not tuned: a bond WF "
                         "shared by three or more atoms of one element is a "
                         "multi-centre (delocalised cluster) bond; two "
                         "centres is an ordinary bond"],
    },
    "table_free": False,
}


# ---------------------------------------------------------------- geometry
class _Env:
    """Atom images around one WF centre, nearest first, as tuples
    (distance, atom, vector centre->image, image key). Images are enumerated
    exactly within a radius from reciprocal-vector bounds (valid for any cell
    shape) and the radius is extended on demand."""

    def __init__(self, fc, st, bnorm, radius):
        self.fc, self.st, self.bn = fc, st, bnorm
        self.R = -1.0
        self.items = []
        self.extend(radius)

    def extend(self, radius):
        if radius <= self.R:
            return
        st, fc, bn = self.st, self.fc, self.bn
        latt = st.lattice
        out = []
        for a, fa in enumerate(st.frac):
            df = [fa[k] - fc[k] for k in range(3)]
            df = [x - round(x) for x in df]
            rng = [range(math.ceil(-radius * bn[k] - df[k]),
                         math.floor(radius * bn[k] - df[k]) + 1) for k in range(3)]
            for n0 in rng[0]:
                for n1 in rng[1]:
                    for n2 in rng[2]:
                        f = (df[0] + n0, df[1] + n1, df[2] + n2)
                        v = [f[0] * latt[0][j] + f[1] * latt[1][j]
                             + f[2] * latt[2][j] for j in range(3)]
                        d = math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
                        if d <= radius:
                            out.append((d, a, v, (a, n0, n1, n2)))
        out.sort(key=lambda t: t[0])
        self.items, self.R = out, radius


def _reach(k, cur):
    """Triangle inequality: an image farther than (2k+1) x the current
    centroid distance of a k-member set cannot bring the centroid closer."""
    return (k + k + 1) * cur


def _gap(p, q):
    return math.sqrt(sum((p[2][j] - q[2][j]) ** 2 for j in range(3)))


def _span(members):
    """Largest member-member distance of a host set."""
    return max([_gap(p, q) for i, p in enumerate(members)
                for q in members[:i]] or [0.0])


def _is_bond(members, env):
    """Length scale of a host set, taken from the structure itself. The bond
    sphere is centred on the centroid of the member images and has the
    largest member-member distance as diameter. The set is a bond of the WF
    only if
      - the WF centre lies inside that sphere,
      - no other atom does (the members are neighbours with nothing
        between them), and
      - the WF is closer to every member than any other atom is to that
        member (a bond WF sits inside its atoms' own neighbour shells; this
        is what rules out a pair of atoms facing each other across
        vacuum)."""
    k = len(members)
    cen = [sum(m[2][j] for m in members) / k for j in range(3)]
    s = math.sqrt(sum(x * x for x in cen))
    rad = _span(members) / 2          # radius = half the diameter
    if not s < rad:
        return False
    far = max(m[0] for m in members)
    env.extend(max(s + rad, far + far))
    keys = set(m[3] for m in members)
    for t in env.items:
        if t[0] >= s + rad:
            break
        if t[3] not in keys and math.sqrt(sum(
                (t[2][j] - cen[j]) ** 2 for j in range(3))) < rad - SHELL_TOL:
            return False
    for m in members:
        for t in env.items:
            if t[0] >= m[0] + m[0]:
                break
            if t[3] != m[3] and _gap(t, m) < m[0] - SHELL_TOL:
                return False
    return True


def _grow(hosts, env, elements, el):
    """Images to add to a host set: the single image (of element el, or any
    element when el is None) that brings the centroid closest to the WF
    centre, plus images tied with it; None when no image qualifies. An image
    qualifies when it is not an image of a member atom, brings the centroid
    closer, and the enlarged set is a bond of the WF (_is_bond)."""
    k = len(hosts)
    sums = [sum(h[2][j] for h in hosts) for j in range(3)]
    cur = math.sqrt(sum(x * x for x in sums)) / k
    reach = _reach(k, cur)
    env.extend(reach)
    taken = set(h[1] for h in hosts)
    scored = []
    for c in env.items:
        if c[0] >= reach:
            break
        if c[1] in taken or (el is not None and elements[c[1]] != el):
            continue
        s = math.sqrt(sum((sums[j] + c[2][j]) ** 2 for j in range(3))) / (k + 1)
        if s < cur:
            scored.append((s, c))
    scored.sort(key=lambda t: t[0])
    best, add = None, []
    for s, c in scored:
        if best is not None and s - best > SHELL_TOL:
            break
        if _is_bond(hosts + [c], env):
            best = s if best is None else best
            add.append(c)
    if not add:
        return None
    if len(add) > 1:
        if len(set(c[1] for c in add)) < len(add):
            return None
        if not _is_bond(hosts + add, env):
            return None
    return add


def _host_set(env, elements):
    """Lewis host set of one WF: (list of images, is_bond).

    Lone pair unless the centre is closer to the midpoint of its nearest
    atom and some partner than to that atom and the pair is a bond of the WF
    (_is_bond); a same-element set then keeps absorbing same-element atoms
    one at a time while that brings the set's centroid closer, and finally
    one atom of another element if that still does. Atoms equidistant from
    the centre start together; if they are not a bond of the WF they share
    it as a lone pair (is_bond False)."""
    d1 = env.items[0][0]
    hosts, seen = [], set()
    for c in env.items:
        if c[0] - d1 > SHELL_TOL:
            break
        if c[1] not in seen:
            seen.add(c[1])
            hosts.append(c)
    if len(hosts) == 1:
        add = _grow(hosts, env, elements, None)
        if add is None:
            return hosts, False
        hosts = hosts + add
    elif not _is_bond(hosts, env):
        return hosts, False
    els = {elements[h[1]] for h in hosts}
    if len(els) > 1:
        return hosts, True              # heteronuclear: the pair decides
    el = next(iter(els))
    while True:
        add = _grow(hosts, env, elements, el)
        if add is None:
            break
        hosts = hosts + add
    # a same-element set that an atom of another element can still join is a
    # heteronuclear multi-centre bond: the most electronegative member wins,
    # whichever of its atoms happens to be nearest
    add = _grow(hosts, env, elements, None)
    return (hosts if add is None else hosts + add), True


def _occupations(case, counts):
    """Per-WF occupation in channel-major order, or None unless the EOSCAR
    #WF block covers every WF of every channel exactly once."""
    wd = getattr(case.eoscar, "wf_detail", None)
    if not wd or len(wd) != sum(counts):
        return None
    offs = [sum(counts[:c]) for c in range(len(counts))]
    occ = [None] * sum(counts)
    for w in wd:
        c, p = w.spin - 1, w.index - 1
        if not (0 <= c < len(counts) and 0 <= p < counts[c]):
            return None
        if occ[offs[c] + p] is not None:
            return None
        occ[offs[c] + p] = float(w.occ)
    return occ


def _analyse(case):
    """Geometry and IUPAC classification of every WF."""
    st = case.structure
    els = st.elements
    na = len(els)
    binv = inv3(st.lattice)
    zval = case.eoscar.zval
    chi = [eneg(e) for e in els]
    chi_max = max(chi) or 1.0
    bounds = [capacity_bounds(e, zval[e]) for e in els]
    s_unit = 2 if case.ispin == 1 else 1
    occ = _occupations(case, [len(ch.centres) for ch in case.channels])

    bnorm = [math.sqrt(sum(binv[i][j] ** 2 for i in range(3))) for j in range(3)]
    wfs = []
    for c, ch in enumerate(case.channels):
        for cen in ch.centres:
            fc = [sum(cen[i] * binv[i][j] for i in range(3)) for j in range(3)]
            d2 = [_min_image_d2(fc, st.frac[a], st.lattice) for a in range(na)]
            # every image that can matter for the first host-set step lies
            # within the reach of a one-member set at the nearest distance
            env = _Env(fc, st, bnorm,
                       _reach(1, math.sqrt(min(d2))) + SHELL_TOL + SHELL_TOL)
            hosts, bond = _host_set(env, els)
            i = len(wfs)
            wfs.append({"idx": i, "s": s_unit, "e": s_unit, "q": s_unit,
                        "full": s_unit, "hosts": hosts, "all": hosts,
                        "bond": bond, "d2": d2, "ch": c, "fc": fc,
                        "near": sorted(set(t[1] for t in env.items)),
                        "occ": None if occ is None
                        else min(float(s_unit), max(0.0, occ[i]))})

    for w in wfs:
        atoms = [h[1] for h in w["hosts"]]
        hel = {els[a] for a in atoms}
        if len(atoms) == 1:
            w["kind"] = "lp"
        elif len(hel) == 1:
            w["kind"] = "homo"
        else:
            w["kind"] = "hetero"
            top = max(chi[a] for a in atoms)
            winners = [h for h in w["hosts"] if chi[h[1]] == top]
            if len(winners) > 1:
                w["kind"] = "homo"      # electronegativity tie: share
            w["hosts"] = winners
    return st, els, na, chi, chi_max, bounds, wfs, zval, binv


def _image_key(st, binv, a, va, b, vb):
    """Integer lattice translation of image b relative to image a."""
    dv = [vb[j] - va[j] for j in range(3)]
    df = [sum(dv[i] * binv[i][j] for i in range(3)) for j in range(3)]
    return (b,) + tuple(round(df[k] - (st.frac[b][k] - st.frac[a][k]))
                        for k in range(3))


def _bond_key(st, binv, w):
    """Identity of the bond a WF sits in: the images of the atoms that share
    it (all members when a single atom wins it) up to a common lattice
    translation."""
    imgs = sorted(w["hosts"] if len(w["hosts"]) > 1 else w["all"],
                  key=lambda h: h[1])
    a0, v0 = imgs[0][1], imgs[0][2]
    return (w["ch"],) + tuple(sorted(_image_key(st, binv, a0, v0, h[1], h[2])
                                     for h in imgs))


def _int_costs(wfs, na, chi, chi_max):
    """Shipped transport cost d^2 + beta (1 - chi / chi_max) in integer units
    of SHELL_TOL^2 (distances rounded to SHELL_TOL), so equal costs are
    exactly equal and the flow arithmetic is exact."""
    unit = SHELL_TOL * SHELL_TOL
    eterm = [round(BETA * (1.0 - chi[a] / chi_max) / unit) for a in range(na)]
    for w in wfs:
        w["cost"] = [round(math.sqrt(w["d2"][a]) / SHELL_TOL) ** 2 + eterm[a]
                     for a in range(na)]
        hs = [h[1] for h in w["hosts"]]
        w["own"] = sum(w["cost"][a] for a in hs) / len(hs)


# ---------------------------------------------------------------- rounding
def _round_share(E, members, quota, rank):
    """Integer split of E electrons over members closest to the fractional
    quotas (largest remainder; rank breaks equal remainders). Returns
    (target, tied, need): tied members still compete for need electrons.
    rank values are tuples of integers compared lexicographically."""
    def grid(x):
        return round(x / SHELL_TOL)
    fl = {a: math.floor(quota[a] + SHELL_TOL) for a in members}
    R = E - sum(fl.values())
    target = dict(fl)
    if R <= 0:
        return target, [], 0
    key = {a: (grid(quota[a] - fl[a]),) + tuple(rank[a]) for a in members}
    order = sorted(members, key=lambda a: key[a], reverse=True)
    cut = key[order[min(R, len(order)) - 1]]
    above = [a for a in members if key[a] > cut]
    at = [a for a in members if key[a] == cut]
    for a in above:
        target[a] += 1
    need = R - len(above)
    if need >= len(at):
        for a in at:
            target[a] += 1
        return target, [], 0
    return target, at, need


# ---------------------------------------------------------------- surplus
def _shares(wfs, na, field):
    out = [0.0] * na
    for w in wfs:
        for h in w["hosts"]:
            out[h[1]] += w[field] / len(w["hosts"])
    return out


def _cap_by_cost(st, binv, na, bounds, wfs, holes):
    """No occupations: remove `holes` electrons one at a time. Order:
    (0) a WF counted in a valence shell that is over n_max even after the
        octet rule - what the shell cannot hold is the surplus;
    (1) otherwise the costliest own Lewis destination;
    (2) last, the only WF left in a bond (per spin channel): the Lewis
        skeleton is kept while anything else can be removed - unless every
        owner of that bond is over its shell, then it counts as (0).
    Inside a tier the costliest goes first; n_min is respected while
    possible. Sets w['e'] (electrons left on each WF)."""
    # valence shells as the octet rule would leave them: an excess that the
    # rule can turn into partner lone pairs (resonance) is not surplus
    trial = [dict(w) for w in wfs]
    _octet_repair(st, binv, na, bounds, trial)
    shell = [0] * na
    for w, t in zip(wfs, trial):
        w["own_set"] = sorted(set(h[1] for h in t["hosts"]))
        for a in w["own_set"]:
            shell[a] += w["s"]
    left = {}
    for w in wfs:
        w["bkey"] = _bond_key(st, binv, w) if w["bond"] else None
        if w["bond"]:
            left[w["bkey"]] = left.get(w["bkey"], 0) + 1
    sign = -1 if _REVERSE_TIES else 1
    order = sorted(wfs, key=lambda w: (-w["own"], sign * w["idx"]))
    held = _shares(wfs, na, "s")
    for respect in (True, False):
        while holes > 0:
            pick, tier = None, None
            for w in order:
                if w["e"] <= 0:
                    continue
                part = 1.0 / len(w["hosts"])
                if respect and any(held[h[1]] - part < bounds[h[1]][0] - SHELL_TOL
                                   for h in w["hosts"]):
                    continue
                over = [shell[a] > bounds[a][1] for a in w["own_set"]]
                if w["bond"] and left[w["bkey"]] == 1 and not all(over):
                    t = 2
                else:
                    t = 0 if any(over) else 1
                if tier is None or t < tier:
                    pick, tier = w, t
                if t == 0:
                    break
            if pick is None:
                break
            pick["e"] -= 1
            holes -= 1
            for a in pick["own_set"]:
                shell[a] -= 1
            for h in pick["hosts"]:
                held[h[1]] -= 1.0 / len(pick["hosts"])
            if pick["e"] == 0 and pick["bond"]:
                left[pick["bkey"]] -= 1
        if holes <= 0:
            break


def _classes(els, na, wfs):
    """Lewis classes: atoms of one element whose own WFs form the same
    multiset of kinds (lone pair / heteronuclear bond won / k-centre
    same-element set) in their spin channels. IUPAC counting is a function
    of exactly these data; the two channels are compared as an unordered
    pair, so antiferromagnetic partners are one class."""
    tok = [dict() for _ in range(na)]
    for w in wfs:
        k = len(w["hosts"])
        for a in set(h[1] for h in w["hosts"]):
            if k > 1:
                t = ("h", k)
            else:
                t = ("l",) if len(w["all"]) == 1 else ("w",)
            tok[a].setdefault(w["ch"], []).append(t)
    groups = {}
    for a in range(na):
        sig = (els[a], tuple(sorted(tuple(sorted(v)) for v in tok[a].values())))
        groups.setdefault(sig, []).append(a)
    return list(groups.values())


def _pool_classes(n, lost, classes):
    """Atoms of one Lewis class lose surplus electrons equally: each member
    gets (its count with nothing removed) - (class mean of the removal),
    rounded inside the class by largest remainder (class total conserved;
    equal remainders keep the cost ranking's own preference)."""
    n = list(n)
    for cls in classes:
        if len(cls) < 2:
            continue
        mean = sum(lost[a] for a in cls) / len(cls)
        quota = {a: n[a] + lost[a] - mean for a in cls}
        tgt, tied, need = _round_share(sum(n[a] for a in cls), cls, quota,
                                       {a: (n[a],) for a in cls})
        for a in sorted(tied)[:need]:
            tgt[a] += 1
        for a in cls:
            n[a] = tgt[a]
    return n


def _cation_sites(st, els, na, chi):
    """Atoms that are not the more electronegative partner of their closest
    heteronuclear contact (ionic approximation: an empty state shared by a
    bonded pair belongs to the less electronegative partner)."""
    out = []
    for a in range(na):
        near = [(_min_image_d2(st.frac[a], st.frac[b], st.lattice), b)
                for b in range(na) if els[b] != els[a]]
        if not near:
            out.append(True)
            continue
        d0 = math.sqrt(min(near)[0])
        out.append(all(chi[a] <= chi[b] for d2, b in near
                       if math.sqrt(d2) - d0 <= SHELL_TOL))
    return out


def _holes_by_depletion(holes, dep, n, bounds, cation):
    """Integer holes per atom from per-atom depletion sums: the sums of the
    eligible atoms are scaled to the hole count and rounded by largest
    remainder; cation sites first, any atom afterwards; an atom never drops
    below n_min. None when the depletion cannot place every hole."""
    na = len(n)
    h = [0] * na
    for tier in (True, False):
        while holes > 0:
            mem = [a for a in range(na) if (cation[a] or not tier)
                   and n[a] - h[a] > bounds[a][0] and dep[a] > SHELL_TOL]
            tot = sum(dep[a] for a in mem)
            if not mem or tot <= SHELL_TOL:
                break
            quota = {a: dep[a] * holes / tot for a in mem}
            full = [a for a in mem
                    if quota[a] > n[a] - h[a] - bounds[a][0] + SHELL_TOL]
            if full:                      # saturate, then share the rest
                for a in full:
                    room = n[a] - h[a] - bounds[a][0]
                    h[a] += room
                    holes -= room
                continue
            tgt, tied, need = _round_share(holes, mem, quota,
                                           {a: (round(dep[a] / SHELL_TOL),)
                                            for a in mem})
            for a in sorted(tied)[:need]:
                tgt[a] += 1
            for a in mem:
                h[a] += tgt[a]
            holes = 0
    return h if holes == 0 else None


# ---------------------------------------------------------------- Lewis
def _octet_repair(st, binv, na, bounds, wfs):
    """Convert excess same-element pairs of over-bonded atoms into partner
    lone pairs (Lewis octet rule with n_max from capacity_bounds).

    Valence shell D(a): every WF whose (IUPAC-reduced) host set contains a,
    counted with the electrons it holds (w['q']) - its lone pairs, the
    heteronuclear bonds it wins and its same-element bonds. Resonance only
    moves pi / lone-pair electrons, never the skeleton: a bond (per spin
    channel) keeps at least one WF. The receiving partner must itself be
    within its shell; then every eligible WF has exactly one over-full host
    and the processing order cannot matter."""
    D = [0.0] * na
    for w in wfs:
        for a in {h[1] for h in w["hosts"]}:
            D[a] += w["q"]
    over = [D[a] > bounds[a][1] + SHELL_TOL for a in range(na)]
    for a in range(na):
        if not over[a]:
            continue
        X = D[a] - bounds[a][1]
        groups = {}
        for i, w in enumerate(wfs):
            if w["kind"] != "homo" or w.get("demoted"):
                continue
            mine = [h for h in w["hosts"] if h[1] == a]
            if len(mine) != 1 or len(w["hosts"]) < 2:
                continue
            ha = mine[0]
            others = [h for h in w["hosts"] if h is not ha]
            if any(over[h[1]] for h in others):
                continue
            hb = min(others, key=lambda h: h[0])
            r = ha[0] / max(hb[0], SHELL_TOL)
            key = (w["ch"],) + tuple(sorted(
                _image_key(st, binv, a, ha[2], h[1], h[2]) for h in others))
            groups.setdefault(key, []).append((r, i))
        if not groups:
            continue
        glist = [sorted(g, key=lambda t: -t[0]) for g in groups.values()]
        glist.sort(key=lambda g: -g[0][0])
        rnd = 0
        # the last WF of each (bond, channel) group is the sigma skeleton
        while X > SHELL_TOL and any(rnd < len(g) - 1 for g in glist):
            for g in glist:
                if X <= SHELL_TOL:
                    break
                if rnd < len(g) - 1:
                    _r, i = g[rnd]
                    w = wfs[i]
                    w["hosts"] = [h for h in w["hosts"] if h[1] != a]
                    w["demoted"] = True
                    if len(w["hosts"]) == 1:
                        w["kind"] = "lp"
                    D[a] -= w["q"]
                    X -= w["q"]
            rnd += 1


def _lewis_counts(na, wfs):
    """Integer Lewis electron counts per atom. Shared WFs are pooled per
    connected fragment and their equal shares rounded by largest remainder
    (fragment total conserved). A fragment holding a multi-centre bond WF
    (more than two centres) is a delocalised cluster (handled by
    _equalise_clusters); elsewhere a remainder tie is ambiguous and goes to
    the flow. Returns (n, fragments, ambiguous)."""
    n = [0] * na
    parent = list(range(na))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for w in wfs:
        if w["kind"] == "homo":
            atoms = [h[1] for h in w["hosts"]]
            for b in atoms[1:]:
                ra, rb = find(atoms[0]), find(b)
                if ra != rb:
                    parent[max(ra, rb)] = min(ra, rb)
    share = [0.0] * na
    frag_e, frag_w, frag_multi = {}, {}, {}
    for i, w in enumerate(wfs):
        if w["kind"] == "homo":
            k = len(w["hosts"])
            for h in w["hosts"]:
                share[h[1]] += w["s"] / k
            r = find(w["hosts"][0][1])
            frag_e[r] = frag_e.get(r, 0) + w["s"]
            frag_w.setdefault(r, []).append(i)
            frag_multi[r] = frag_multi.get(r, False) or (
                w["bond"] and k >= MULTI_CENTRE)
        else:
            n[w["hosts"][0][1]] += w["s"]

    frags = []
    ambiguous = False
    for r, E in frag_e.items():
        members = [a for a in range(na) if find(a) == r and share[a] > 0]
        cluster = frag_multi[r]
        quota = {a: share[a] for a in members}
        target, tied, need = _round_share(E, members, quota,
                                          {a: () for a in members})
        if tied and cluster:
            # placeholder only: a cluster's count is re-shared equally later
            for a in sorted(tied)[:need]:
                target[a] += 1
            tied, need = [], 0
        ambiguous = ambiguous or bool(tied)
        # cluster tie-break: larger Lewis share first, then the member the
        # fragment's WF centres sit closest to (mean squared host distance)
        d2sum = {a: 0.0 for a in members}
        d2cnt = {a: 0 for a in members}
        for i in frag_w[r]:
            for h in wfs[i]["hosts"]:
                d2sum[h[1]] += h[0] * h[0]
                d2cnt[h[1]] += 1
        frags.append({"members": members, "floor": target, "tied": tied,
                      "need": need, "wfs": frag_w[r], "cluster": cluster,
                      "rank": {a: (round(share[a] / SHELL_TOL),
                                   -round(d2sum[a] / max(1, d2cnt[a])
                                          / (SHELL_TOL * SHELL_TOL)))
                               for a in members}})
        for a in members:
            n[a] += target[a]
    return n, frags, ambiguous


def _equalise_clusters(n, frags):
    """Delocalised same-element cluster (fragment with a multi-centre WF):
    which of its electrons sit in one-, two- or many-centre WFs is a gauge
    choice of the localisation, its total is not. Share the members' total
    equally; extras go to the members with the largest Lewis shares."""
    n = list(n)
    for fr in frags:
        if not fr["cluster"]:
            continue
        mem = fr["members"]
        P = sum(n[a] for a in mem)
        quota = {a: P / len(mem) for a in mem}
        tgt, tied, need = _round_share(P, mem, quota, fr["rank"])
        for a in sorted(tied, key=lambda a: (-n[a], a))[:need]:
            tgt[a] += 1
        for a in mem:
            n[a] = tgt[a]
    return n


# ---------------------------------------------------------------- flow
def _flow(na, wfs, frags, bounds, shift, wide, lifted):
    """Capacity flow: Lewis destinations at the shipped transport cost, the
    atoms around each WF (all atoms when wide) at that cost + LEWIS_PENALTY,
    one capacity arc per atom (n_min rewarded; n_max dropped when lifted).
    Costs are integers. Returns (electrons per atom, electrons not placed)."""
    unit = SHELL_TOL * SHELL_TOL
    pen = round(LEWIS_PENALTY / unit)
    rew = round(SEMICORE_REWARD / unit)
    nw, nf = len(wfs), len(frags)
    S = 0
    W0 = 1
    F0 = W0 + nw
    A0 = F0 + nf
    T = A0 + na
    mcf = MinCostFlow(T + 1)
    frag_of = {}
    for f, fr in enumerate(frags):
        for i in fr["wfs"]:
            frag_of[i] = f
    total = 0
    for i, w in enumerate(wfs):
        total += w["s"]
        mcf.add(S, W0 + i, w["s"], 0)
        if i in frag_of:
            mcf.add(W0 + i, F0 + frag_of[i], w["s"], 0)
        else:
            a = w["hosts"][0][1]
            mcf.add(W0 + i, A0 + a, w["s"], w["cost"][a])
        alts = range(na) if wide else sorted(
            set(w["near"]) | set(h[1] for h in w["all"]))
        for a in alts:
            mcf.add(W0 + i, A0 + a, w["s"], w["cost"][a] + pen)
    for f, fr in enumerate(frags):
        idx = fr["wfs"]
        for a in fr["members"]:
            ws = [i for i in idx if any(h[1] == a for h in wfs[i]["hosts"])]
            c = round(sum(wfs[i]["cost"][a] for i in ws) / max(1, len(ws)))
            cap = fr["floor"][a] + (1 if a in fr["tied"] else 0)
            if cap > 0:
                mcf.add(F0 + f, A0 + a, cap, c)
    for a in range(na):
        lo, hi = bounds[a][0] + shift[a], bounds[a][1] + shift[a]
        if lifted:
            hi = max(hi, total)
        if lo > 0:
            mcf.add(A0 + a, T, lo, -rew)
        if hi > lo:
            mcf.add(A0 + a, T, hi - max(lo, 0), 0)
    mcf.solve(S, T, max_flow=total)
    # electrons into atom a = flow on its atom->sink arcs (reverse caps)
    n = [sum(mcf.g[T][rev][1] for v, _cap, _c, rev in mcf.g[A0 + a] if v == T)
         for a in range(na)]
    return n, total - sum(n)


def _count(st, binv, na, bounds, wfs, shift):
    """Octet repair, Lewis counts, clusters and (only if a bound is broken
    or a remainder tie is left) the capacity flow, for WFs holding w['s']
    electrons. shift[a] = holes to be taken from atom a afterwards."""
    _octet_repair(st, binv, na, bounds, wfs)
    n, frags, ambiguous = _lewis_counts(na, wfs)
    n = _equalise_clusters(n, frags)
    ok = not ambiguous and all(
        bounds[a][0] <= n[a] - shift[a] <= bounds[a][1] for a in range(na))
    if not ok:
        for wide, lifted in ((False, False), (True, False), (True, True)):
            n, left = _flow(na, wfs, frags, bounds, shift, wide, lifted)
            if left == 0:
                break
        n = _equalise_clusters(n, frags)
    return n, frags


def _with_occupations(analysis, target):
    """Surplus settled from the #WF occupations: WFs that round to zero
    electrons are absent, the other holes are placed on atoms from their
    depletion sums. None when the occupations cannot account for the holes."""
    st, els, na, chi, chi_max, bounds, wfs, zval, binv = analysis
    _int_costs(wfs, na, chi, chi_max)
    live = [w for w in wfs if round(w["occ"]) > 0]
    holes = sum(w["s"] for w in live) - target
    if holes < 0:
        return None
    for w in live:
        w["q"] = w["occ"]
        w["gap"] = w["s"] - w["occ"]
    full, frags = _count(st, binv, na, bounds, live, [0] * na)
    h = _holes_by_depletion(holes, _shares(live, na, "gap"), full, bounds,
                            _cation_sites(st, els, na, chi))
    if h is None:
        return None
    if any(not bounds[a][0] <= full[a] - h[a] <= bounds[a][1] for a in range(na)):
        for wide, lifted in ((False, False), (True, False), (True, True)):
            full, left = _flow(na, live, frags, bounds, h, wide, lifted)
            if left == 0:
                break
    return _equalise_clusters([full[a] - h[a] for a in range(na)], frags)


def assign(case):
    analysis = _analyse(case)
    wfs = analysis[6]
    supply = sum(w["s"] for w in wfs)
    target = supply
    if case.eoscar.nelect_wann is not None:
        target = min(round(case.eoscar.nelect_wann), supply)
    n = None
    if supply > target and all(w["occ"] is not None for w in wfs):
        n = _with_occupations(analysis, target)
        if n is None:
            analysis = _analyse(case)
    st, els, na, chi, chi_max, bounds, wfs, zval, binv = analysis
    if n is None:
        _int_costs(wfs, na, chi, chi_max)
        lost = [0.0] * na
        classes = []
        if supply > target:
            # no occupations: cost-ranked removal, shared inside Lewis classes
            _cap_by_cost(st, binv, na, bounds, wfs, supply - target)
            for w in wfs:
                w["gap"] = w["s"] - w["e"]
            lost = _shares(wfs, na, "gap")
            classes = _classes(els, na, wfs)
        live = [w for w in wfs if w["e"] > 0]
        for w in live:
            w["s"] = w["q"] = w["e"]
        n, frags = _count(st, binv, na, bounds, live, [0] * na)
        n = _equalise_clusters(_pool_classes(n, lost, classes), frags)
    return [round(zval[els[a]]) - n[a] for a in range(na)]
