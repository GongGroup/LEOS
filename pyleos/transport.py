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
    os_int: List[int]
    wfs: List[WfInfo]
    stray: List[WfInfo]
    supply_total: int        # 2*n_wann (or n_wann per spin channel)
    target: int              # electrons actually shipped target (NELECT cap)
    ghost: int               # supply beyond NELECT excluded by the flow cap
    infeasible: int          # target electrons that found no capacity (0 = OK)
    transport_cost: float


def supply_total_of(wfs: List["WfInfo"]) -> int:
    return sum(w.supply for w in wfs)


def assign_case(case: CaseData, beta: float = 3.0,
                eneg_scale: str = "allen",
                drop_stray: bool = False) -> AssignResult:
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

    # nodes: 0 source | 1..nw WFs | nw+1..nw+na atoms | nw+na+1 sink
    S, T = 0, nw + na + 1
    mcf = MinCostFlow(nw + na + 2)
    for w, wf in enumerate(wfs):
        mcf.add(S, 1 + w, wf.supply, 0.0)
        for a in range(na):
            c = d2rows[w][a] + beta * (1.0 - eneg(st.elements[a], eneg_scale) / chi_max)
            mcf.add(1 + w, 1 + nw + a, wf.supply, c)
    for a in range(na):
        n_min, n_max = capacity_bounds(st.elements[a], zval[st.elements[a]])
        if n_min > 0:
            mcf.add(1 + nw + a, T, n_min, -BIG_NEG)
        if n_max > n_min:
            mcf.add(1 + nw + a, T, n_max - n_min, 0.0)
    # Ghost-WF guard: num_wann can exceed the truly occupied manifold, and
    # the full-occupancy assumption then injects electrons that do not exist.
    # Capping total flow at NELECT(Wannier) keeps the min-cost k-flow optimal
    # (SSP ships cheapest units first), so ghost WFs are excluded naturally.
    target = None
    if case.eoscar.nelect_wann is not None:
        target = min(round(case.eoscar.nelect_wann), supply_total_of(wfs))
    elif drop_stray:
        target = supply_total_of(wfs)
    cost = mcf.solve(S, T, max_flow=target)

    ne = [0] * na
    for w, wf in enumerate(wfs):
        for v, cap_left, c, _ in mcf.g[1 + w]:
            if 1 + nw <= v <= nw + na and c >= 0:
                sent = wf.supply - cap_left
                if sent > 0:
                    a = v - 1 - nw
                    ne[a] += sent
                    wf.dests.append((a, sent))

    routed = sum(ne)
    supply_total = sum(w.supply for w in wfs)
    tgt = target if target is not None else supply_total
    os_int = [round(zval[st.elements[a]]) - ne[a] for a in range(na)]
    stray = [w for w in wfs
             if w.spread > STRAY_SPREAD or w.d1 > STRAY_DIST]
    return AssignResult(ne, os_int, wfs, stray, supply_total, tgt,
                        supply_total - tgt, tgt - routed, cost)
