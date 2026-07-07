"""Capacity-constrained min-cost-flow assignment of Wannier electrons to atoms.

Graph: source -> WF nodes (supply 2e for ISPIN=1, 1e per spin channel for
ISPIN=2) -> atom nodes -> sink. Atom capacity [n_min, n_max] comes from the
wide physical bounds in capacity.py; the n_min part is enforced with a
strongly negative-cost edge segment (semicore shells must stay filled).
Integer supplies and capacities guarantee an integer optimal flow, so
oxidation states come out integer without a rounding step.
"""
from __future__ import annotations

import heapq
import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from .capacity import capacity_bounds, eneg
from .parsers import CaseData, inv3

# Reward for filling mandatory n_min capacity (semicore shells). Must beat any
# geometric detour (max d^2 + beta ~ 1e2) but stay small enough that potentials
# keep float headroom above the comparison tolerance: at 1e5 scale the ULP
# (~1e-11) reaches the old 1e-12 tolerance, so cancellation residue registered
# as fake improvements and Dijkstra stalled re-pushing the same nodes.
BIG_NEG = 1.0e3
STRAY_SPREAD = 3.0  # Angstrom^2; WFs above this are flagged as unlocalised
STRAY_DIST = 1.6    # Angstrom; centre farther than this from every atom
                    # (must clear homonuclear bond midpoints, e.g. As-As/2 = 1.28)


def _min_image_d2(fc: List[float], fa: List[float], latt: List[List[float]]) -> float:
    df = [fc[k] - fa[k] for k in range(3)]
    df = [x - round(x) for x in df]
    best = math.inf
    for n0 in (-1, 0, 1):
        for n1 in (-1, 0, 1):
            for n2 in (-1, 0, 1):
                f = (df[0] + n0, df[1] + n1, df[2] + n2)
                cart = [sum(f[k] * latt[k][j] for k in range(3)) for j in range(3)]
                d2 = cart[0] ** 2 + cart[1] ** 2 + cart[2] ** 2
                if d2 < best:
                    best = d2
    return best


class MinCostFlow:
    """Successive shortest paths with potentials; Bellman-Ford bootstrap
    so that a single layer of negative-cost edges (n_min rewards) is OK."""

    def __init__(self, n: int):
        self.n = n
        self.g: List[List[List[float]]] = [[] for _ in range(n)]
        self.pot: List[float] = []      # final SSP potentials (dual solution)

    def add(self, u: int, v: int, cap: int, cost: float) -> None:
        self.g[u].append([v, cap, cost, len(self.g[v])])
        self.g[v].append([u, 0, -cost, len(self.g[u]) - 1])

    def _bellman_ford(self, s: int) -> List[float]:
        dist = [math.inf] * self.n
        dist[s] = 0.0
        for _ in range(self.n):
            changed = False
            for u in range(self.n):
                if dist[u] == math.inf:
                    continue
                for v, cap, cost, _ in self.g[u]:
                    if cap > 0 and dist[u] + cost < dist[v] - 1e-9:
                        dist[v] = dist[u] + cost
                        changed = True
            if not changed:
                break
        return dist

    def solve(self, s: int, t: int, max_flow: Optional[int] = None) -> float:
        pot = self._bellman_ford(s)
        for v in range(self.n):
            if pot[v] == math.inf:
                pot[v] = 0.0
        total = 0.0
        shipped = 0
        supply = sum(e[1] for e in self.g[s])
        if max_flow is None:
            max_flow = supply
        max_rounds = 2 * supply + 100        # each augmentation ships >=1 unit
        rounds = 0
        while shipped < max_flow:
            rounds += 1
            if rounds > max_rounds:
                raise RuntimeError(
                    f"min-cost-flow failed to terminate within {max_rounds} "
                    f"augmentations (supply {supply}) - solver bug, not data")
            dist = [math.inf] * self.n
            dist[s] = 0.0
            prev: List[Optional[Tuple[int, int]]] = [None] * self.n
            pq: List[Tuple[float, int]] = [(0.0, s)]
            while pq:
                d, u = heapq.heappop(pq)
                if d > dist[u] + 1e-9:
                    continue
                for i, (v, cap, cost, _) in enumerate(self.g[u]):
                    if cap <= 0:
                        continue
                    nd = d + cost + pot[u] - pot[v]
                    if nd < dist[v] - 1e-9:
                        dist[v] = nd
                        prev[v] = (u, i)
                        heapq.heappush(pq, (nd, v))
            if prev[t] is None:
                break
            # Potential update must also lift unreached nodes (standard
            # Johnson/SSP treatment): reusing a stale potential when a node
            # becomes reachable again creates negative reduced costs and
            # degrades Dijkstra to exponential behaviour.
            dmax = max(d for d in dist if d < math.inf)
            for v in range(self.n):
                pot[v] += dist[v] if dist[v] < math.inf else dmax
            f = max_flow - shipped
            v = t
            while v != s:
                u, i = prev[v]
                f = min(f, self.g[u][i][1])
                v = u
            shipped += f
            v = t
            while v != s:
                u, i = prev[v]
                self.g[u][i][1] -= f
                self.g[v][self.g[u][i][3]][1] += f
                total += f * self.g[u][i][2]
                v = u
        self.pot = pot
        return total


@dataclass
class WfInfo:
    channel: int
    index: int              # 1-based within channel
    spread: float
    supply: int             # electrons carried
    dests: List[Tuple[int, int]] = field(default_factory=list)  # (atom, e)
    d1: float = 0.0
    nearest: int = -1


@dataclass
class AssignResult:
    ne: List[int]
    os_int: List[Optional[int]]     # None for atoms folded into a fragment
    wfs: List[WfInfo]
    stray: List[WfInfo]
    supply_total: int        # 2*n_wann (or n_wann per spin channel)
    target: int              # electrons actually shipped target (NELECT cap)
    ghost: int               # supply beyond NELECT excluded by the flow cap
    infeasible: int          # target electrons that found no capacity (0 = OK)
    transport_cost: float
    # Dual-margin reliability (counterpart of Comas-Vila's R(%)): how much
    # extra transport cost the cheapest single rerouting would incur.
    wf_margin: List[float] = field(default_factory=list)    # per WF (arc bound)
    atom_margin: List[float] = field(default_factory=list)  # site level (cost units)
    atom_conf: List[float] = field(default_factory=list)    # site level, 0-100
    case_conf: float = 100.0                                # min site R
    # Element-ledger level: cheapest +-1 rearrangement that CHANGES some
    # element's total electron count. Intra-element shuffles (homonuclear
    # bond splits, equivalent sites) are excluded - they cannot move the
    # element-mean OS that the dev benchmark scores.
    atom_margin_elem: List[float] = field(default_factory=list)
    case_conf_elem: float = 100.0
    # Manual fragments (delocalised multi-centre units, e.g. I3-): each entry
    # is the 0-based atom-index list of one merged warehouse. Members report
    # os_int None; the unit's total charge state lives in frag_os.
    fragments: List[List[int]] = field(default_factory=list)
    frag_ne: List[int] = field(default_factory=list)
    frag_os: List[int] = field(default_factory=list)


def supply_total_of(wfs: List["WfInfo"]) -> int:
    return sum(w.supply for w in wfs)


# Cost-gap below which tier-2 composition may arbitrate. Chosen as the widest
# zero-damage point on the 38-case #WP bench: tau=1.0 fixes hi3 + the Rh
# complex and breaks nothing; 1.5 starts overriding clear geometry (tmao).
ARBITER_TAU = 1.0
ARBITER_PMARGIN = 0.05  # composition lead required to overturn the geometry


def _build_costs(case: CaseData, beta: float, eneg_scale: str,
                 drop_stray: bool, gamma: float, tier2: str = "linear"):
    """WF geometry + per-(wf,atom) transport costs + capacity bounds + flow
    target. Split out so reliability_margins can re-solve with perturbed
    bounds without recomputing the (dominant) minimum-image geometry.

    tier2 modes (all need the EOSCAR v3 #WP block):
      linear   - historical v1 formula: c -= gamma * p_hat (p_hat per-WF
                 normalised). Fails systematically: heavy atoms carry more
                 and stronger projectors, so raw shares drag every electron
                 towards them.
      arbiter  - tier-2 information only breaks near-ties: when the two
                 cheapest hosts of a WF sit within ARBITER_TAU of each other
                 and the composition clearly prefers the runner-up, nudge it
                 just below the leader. Clear geometry is never overridden.
      baseline - per-atom renormalisation: an atom's share only counts to
                 the extent it exceeds that atom's own average share across
                 all WFs (removing the heavy-atom bias at its source), then
                 c -= gamma * excess."""
    st = case.structure
    na = len(st.elements)
    binv = inv3(st.lattice)
    if case.eoscar is None:
        raise ValueError(f"{case.path}: EOSCAR required for ZVAL")
    zval = case.eoscar.zval

    wfs: List[WfInfo] = []
    d2rows: List[List[float]] = []
    for ch_idx, ch in enumerate(case.channels):
        supply = 2 if case.ispin == 1 else 1
        for w, centre in enumerate(ch.centres):
            fc = [sum(centre[i] * binv[i][j] for i in range(3)) for j in range(3)]
            row = [_min_image_d2(fc, st.frac[a], st.lattice) for a in range(na)]
            wfs.append(WfInfo(ch_idx + 1, w + 1, ch.spreads[w], supply))
            d2rows.append(row)
            wfs[-1].d1 = math.sqrt(min(row))
            wfs[-1].nearest = row.index(min(row))

    if drop_stray:
        # Attribution mode: exclude unlocalised WFs from the supply so the
        # assignment quality of the chemically meaningful WFs can be judged
        # in isolation (their electrons are simply not distributed).
        keep = [i for i, w in enumerate(wfs)
                if w.spread <= STRAY_SPREAD and w.d1 <= STRAY_DIST]
        wfs = [wfs[i] for i in keep]
        d2rows = [d2rows[i] for i in keep]

    nw = len(wfs)
    chi_max = max(eneg(e, eneg_scale) for e in set(st.elements)) or 1.0

    phat_rows: List[Dict[int, float]] = []
    use_wp = case.eoscar.wp_detail and (gamma > 0.0 or tier2 == "arbiter")
    for wf in wfs:
        phat = {}
        if use_wp:
            pairs = case.eoscar.wp_detail.get((wf.channel, wf.index), [])
            tot = sum(x for _, x in pairs)
            if tot > 1e-8:
                phat = {at - 1: x / tot for at, x in pairs}
        phat_rows.append(phat)

    if tier2 == "baseline" and use_wp:
        # An atom's average share across all WFs is its "always there"
        # projector background; only the excess above it is chemical signal.
        base = [0.0] * na
        for phat in phat_rows:
            for a, x in phat.items():
                base[a] += x
        base = [b / max(1, nw) for b in base]
        phat_rows = [{a: max(0.0, x - base[a]) for a, x in phat.items()}
                     for phat in phat_rows]

    cost_rows: List[List[float]] = []
    for w in range(nw):
        geo = [d2rows[w][a]
               + beta * (1.0 - eneg(st.elements[a], eneg_scale) / chi_max)
               for a in range(na)]
        if tier2 in ("linear", "baseline"):
            row = [geo[a] - gamma * phat_rows[w].get(a, 0.0) for a in range(na)]
        elif tier2 == "arbiter":
            row = list(geo)
            if phat_rows[w] and na >= 2:
                order = sorted(range(na), key=geo.__getitem__)
                a1, a2 = order[0], order[1]
                if geo[a2] - geo[a1] < ARBITER_TAU:
                    p1 = phat_rows[w].get(a1, 0.0)
                    p2 = phat_rows[w].get(a2, 0.0)
                    if p2 > p1 + ARBITER_PMARGIN:
                        # Composition overrules a near-tie: put the runner-up
                        # just below the geometric leader.
                        row[a2] = geo[a1] - 0.2 * ARBITER_TAU
        else:
            raise ValueError(f"unknown tier2 mode {tier2!r}")
        cost_rows.append(row)

    bounds = [capacity_bounds(el, zval[el]) for el in st.elements]
    target = None
    if case.eoscar.nelect_wann is not None:
        target = min(round(case.eoscar.nelect_wann), supply_total_of(wfs))
    elif drop_stray:
        target = supply_total_of(wfs)
    return wfs, cost_rows, bounds, target, na, nw


def _solve_flow(wfs: List[WfInfo], cost_rows: List[List[float]],
                bounds: List[Tuple[int, int]], target: Optional[int],
                na: int, nw: int):
    """Build the flow network for the given bounds and solve. Returns the
    solver (for duals), per-atom electron counts, per-WF destination lists,
    raw cost, and the geometric cost with the n_min reward term removed."""
    S, T = 0, nw + na + 1
    mcf = MinCostFlow(nw + na + 2)
    for w, wf in enumerate(wfs):
        mcf.add(S, 1 + w, wf.supply, 0.0)
        for a in range(na):
            mcf.add(1 + w, 1 + nw + a, wf.supply, cost_rows[w][a])
    for a in range(na):
        n_min, n_max = bounds[a]
        if n_min > 0:
            mcf.add(1 + nw + a, T, n_min, -BIG_NEG)
        if n_max > n_min:
            mcf.add(1 + nw + a, T, n_max - n_min, 0.0)
    cost = mcf.solve(S, T, max_flow=target)
    ne = [0] * na
    dests: List[List[Tuple[int, int]]] = [[] for _ in wfs]
    for w, wf in enumerate(wfs):
        for v, cap_left, c, _ in mcf.g[1 + w]:
            # Forward wf->atom arcs are exactly those pointing into the atom
            # node range (a sign test on c breaks once gamma makes costs
            # negative, which silently dropped assignments in the tier-2 runs).
            if 1 + nw <= v <= nw + na:
                sent = wf.supply - cap_left
                if sent > 0:
                    a = v - 1 - nw
                    ne[a] += sent
                    dests[w].append((a, sent))
    # Undo the n_min reward so solutions with different bounds are comparable
    # on pure geometry+electronegativity cost.
    reward_flow = sum(min(ne[a], bounds[a][0]) for a in range(na))
    geo_cost = cost + BIG_NEG * reward_flow
    return mcf, ne, dests, cost, geo_cost


def _make_units(na: int, fragments: Optional[List[List[int]]]):
    """Assignment units: one per unmerged atom plus one per fragment (merged
    warehouse: capacities add, an electron pays the fare to its nearest
    member). Returns (units, unit_of_atom)."""
    frags = [sorted(f) for f in (fragments or [])]
    seen: set = set()
    for f in frags:
        for a in f:
            if a in seen:
                raise ValueError(f"atom {a+1} appears in two fragments")
            seen.add(a)
    units = [[a] for a in range(na) if a not in seen] + frags
    unit_of = {a: u for u, members in enumerate(units) for a in members}
    return units, unit_of


def _unit_problem(cost_rows, bounds, units, nw):
    """Collapse atom-level costs/bounds onto units. Identity when every unit
    is a single atom (list order then equals atom order)."""
    if all(len(m) == 1 for m in units):
        return cost_rows, bounds
    u_cost = [[min(cost_rows[w][a] for a in members) for members in units]
              for w in range(nw)]
    u_bounds = [(sum(bounds[a][0] for a in members),
                 sum(bounds[a][1] for a in members)) for members in units]
    return u_cost, u_bounds


def assign_case(case: CaseData, beta: float = 3.0,
                eneg_scale: str = "allen",
                drop_stray: bool = False,
                gamma: float = 0.0,
                fragments: Optional[List[List[int]]] = None,
                tier2: str = "linear") -> AssignResult:
    """gamma > 0 enables the tier-2 cost: c -= gamma * p_hat(w,a), where
    p_hat is the per-WF-normalised projector-norm composition from the
    EOSCAR v3 #WP block (chemical fingerprint; complements the centre
    geometry and stays informative for stray WFs).

    fragments: 0-based atom-index groups treated as single merged warehouses
    (delocalised multi-centre units such as I3- or metal clusters); the group
    reports one total charge state instead of per-member integers."""
    st = case.structure
    zval = case.eoscar.zval if case.eoscar else {}
    wfs, cost_rows, bounds, target, na, nw = _build_costs(
        case, beta, eneg_scale, drop_stray, gamma, tier2)
    for f in (fragments or []):
        for a in f:
            if a < 0 or a >= na:
                raise ValueError(f"fragment atom index {a+1} out of range")
    units, _unit_of = _make_units(na, fragments)
    nu = len(units)
    u_cost, u_bounds = _unit_problem(cost_rows, bounds, units, nw)
    mcf, ne_u, dests, cost, _geo = _solve_flow(wfs, u_cost, u_bounds, target,
                                               nu, nw)
    for w, wf in enumerate(wfs):
        wf.dests = dests[w]        # unit indices when fragments are active

    routed = sum(ne_u)
    supply_total = sum(w.supply for w in wfs)
    tgt = target if target is not None else supply_total
    ne = [0] * na
    os_int: List[Optional[int]] = [None] * na
    frag_list: List[List[int]] = []
    frag_ne: List[int] = []
    frag_os: List[int] = []
    for u, members in enumerate(units):
        if len(members) == 1:
            a = members[0]
            ne[a] = ne_u[u]
            os_int[a] = round(zval[st.elements[a]]) - ne_u[u]
        else:
            frag_list.append(members)
            frag_ne.append(ne_u[u])
            frag_os.append(round(sum(zval[st.elements[a]] for a in members))
                           - ne_u[u])
    stray = [w for w in wfs
             if w.spread > STRAY_SPREAD or w.d1 > STRAY_DIST]

    # Per-WF ambiguity margin from the duals: for an UNUSED wf->atom arc the
    # reduced cost rc = c + pot[wf] - pot[atom] lower-bounds the extra cost of
    # rerouting one unit along it (used arcs sit at rc ~ 0 by complementary
    # slackness). A small wf_margin marks an electron pair with a near-tie
    # alternative host - the fragment-suggestion diagnostic feeds on this.
    # NOTE: this is a per-arc lower bound only; for atom/case-level
    # reliability use reliability_margins(), which re-solves with forced
    # +-1 electron counts and therefore prices the full chain rearrangement.
    pot = mcf.pot
    wf_margin: List[float] = []
    for w, wf in enumerate(wfs):
        used = {u for u, _ in wf.dests}
        m = math.inf
        for v, _cap, c, _ in mcf.g[1 + w]:
            if not (1 + nw <= v <= nw + nu):
                continue
            u = v - 1 - nw
            if u in used:
                continue
            rc = max(0.0, c + pot[1 + w] - pot[v])
            if rc < m:
                m = rc
        wf_margin.append(m)
    return AssignResult(ne, os_int, wfs, stray, supply_total, tgt,
                        supply_total - tgt, tgt - routed, cost,
                        wf_margin, fragments=frag_list, frag_ne=frag_ne,
                        frag_os=frag_os)


def reliability_margins(case: CaseData, beta: float = 3.0,
                        eneg_scale: str = "allen",
                        drop_stray: bool = False,
                        gamma: float = 0.0,
                        res: Optional[AssignResult] = None,
                        fragments: Optional[List[List[int]]] = None,
                        tier2: str = "linear") -> AssignResult:
    """R(%)-counterpart reliability: for every assignment unit (atom or
    fragment), the cheapest total transport-cost increase of forcing its
    electron count to ne+1 or ne-1 (full re-solve, so chain rearrangements
    are priced in, unlike the single-arc wf_margin bound). Fragment members
    share their unit's margin. Fills atom_margin/atom_conf/case_conf (and
    the element-ledger variants) on the AssignResult and returns it.

    Cost units are ~ Angstrom^2 of squared-distance detour; the 0-100
    normalisation R = 100*m/(m+1) maps a 1 A^2 margin to R = 50. A unit
    whose count cannot change at all (capacity-infeasible both ways) gets
    margin inf / R = 100."""
    wfs, cost_rows, bounds, target, na, nw = _build_costs(
        case, beta, eneg_scale, drop_stray, gamma, tier2)
    units, unit_of = _make_units(na, fragments)
    nu = len(units)
    u_cost, u_bounds = _unit_problem(cost_rows, bounds, units, nw)
    _, ne0, _, _, geo0 = _solve_flow(wfs, u_cost, u_bounds, target, nu, nw)
    if res is None:
        res = assign_case(case, beta=beta, eneg_scale=eneg_scale,
                          drop_stray=drop_stray, gamma=gamma,
                          fragments=fragments, tier2=tier2)
    for u, members in enumerate(units):
        if len(members) == 1:
            assert res.ne[members[0]] == ne0[u], \
                "reliability re-solve disagrees with base solve"

    elements = case.structure.elements
    # Ledger label per unit: element symbol, or a pseudo-element per fragment
    # (moving charge between a fragment and anything else always counts as a
    # ledger change; shuffles inside the merged unit are invisible anyway).
    labels = [elements[m[0]] if len(m) == 1 else f"FRAG{u}"
              for u, m in enumerate(units)]

    def _ledger(ne: List[int]) -> Dict[str, int]:
        s: Dict[str, int] = {}
        for x, lb in zip(ne, labels):
            s[lb] = s.get(lb, 0) + x
        return s

    ledger0 = _ledger(ne0)

    def _try(u: int, forced: int) -> Optional[Tuple[float, bool]]:
        n_min, n_max = u_bounds[u]
        # n_min is a PHYSICAL floor (semicore shells must stay filled):
        # raiding it is not a legitimate alternative assignment, so a unit
        # sitting at its floor simply cannot lose an electron.
        if forced < n_min or forced > n_max:
            return None
        b2 = list(u_bounds)
        if forced > ne0[u]:
            b2[u] = (forced, n_max)            # reward edge pulls one in
        else:
            b2[u] = (n_min, forced)            # hard cap pushes one out
        _, ne1, _, _, geo1 = _solve_flow(wfs, u_cost, b2, target, nu, nw)
        if ne1[u] != forced or sum(ne1) != sum(ne0):
            return None                        # force did not take: infeasible
        # Semicore integrity: reject a solution that raids another unit's
        # mandatory n_min shell (rewards cancel in geo, so check explicitly).
        for b in range(nu):
            if b != u and ne1[b] < min(u_bounds[b][0], ne0[b]):
                return None
        delta = geo1 - geo0
        # The perturbed problem is a restriction of the base feasible set, so
        # its optimum can only be costlier; anything below the float noise
        # floor indicates a solver/bookkeeping bug.
        assert delta > -1e-6, f"negative margin {delta} forcing unit {u}"
        return delta, _ledger(ne1) != ledger0

    unit_margin: List[float] = []
    unit_margin_elem: List[float] = []
    for u in range(nu):
        site, elem = [], []
        for forced in (ne0[u] + 1, ne0[u] - 1):
            r = _try(u, forced)
            if r is None:
                continue
            delta, changed = r
            site.append(delta)
            if changed:
                elem.append(delta)
        unit_margin.append(max(0.0, min(site)) if site else math.inf)
        # inf here means: no single-unit rearrangement can move any ledger
        # total - the element-mean OS is immune to +-1 perturbations.
        unit_margin_elem.append(max(0.0, min(elem)) if elem else math.inf)

    def _conf(ms: List[float]) -> List[float]:
        return [100.0 if m == math.inf else 100.0 * m / (m + 1.0) for m in ms]

    res.atom_margin = [unit_margin[unit_of[a]] for a in range(na)]
    res.atom_conf = _conf(res.atom_margin)
    res.case_conf = min(res.atom_conf) if res.atom_conf else 100.0
    res.atom_margin_elem = [unit_margin_elem[unit_of[a]] for a in range(na)]
    conf_elem = _conf(res.atom_margin_elem)
    res.case_conf_elem = min(conf_elem) if conf_elem else 100.0
    return res
