#!/usr/bin/env python3
"""Prototype: capacity-constrained min-cost-flow (transport) assignment of
Wannier centres to atoms, tier-1 geometric cost. Stdlib only.

Validates the transport formulation on a LEOS pathological case by comparing
  (a) current greedy result read from EOSCAR,
  (b) integer min-cost-flow with per-element valence-shell capacities,
against the chemically expected oxidation states.

Usage: python3 ot_assign_proto.py <case_dir> [--omax N_e per element, e.g. S=8]
"""
from __future__ import annotations

import heapq
import sys
from pathlib import Path

# Valence-shell electron capacity (wide/physical bound) and PAW ZVAL used here.
CAPACITY = {"H": 2, "B": 8, "C": 8, "N": 8, "O": 8, "F": 8,
            "S": 8, "As": 8, "Se": 8, "Cl": 8, "I": 8}
ZVAL = {"H": 1, "B": 3, "C": 4, "N": 5, "O": 6, "F": 7,
        "S": 6, "As": 5, "Se": 6, "Cl": 7, "I": 7}
ENEG = {"H": 2.20, "B": 2.04, "C": 2.55, "N": 3.04, "O": 3.44, "F": 3.98,
        "S": 2.58, "As": 2.18, "Se": 2.55, "Cl": 3.16, "I": 2.66}


def mat_inv3(a):
    d = (a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
         - a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
         + a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))
    cof = [[(a[1][1]*a[2][2]-a[1][2]*a[2][1]), -(a[0][1]*a[2][2]-a[0][2]*a[2][1]), (a[0][1]*a[1][2]-a[0][2]*a[1][1])],
           [-(a[1][0]*a[2][2]-a[1][2]*a[2][0]), (a[0][0]*a[2][2]-a[0][2]*a[2][0]), -(a[0][0]*a[1][2]-a[0][2]*a[1][0])],
           [(a[1][0]*a[2][1]-a[1][1]*a[2][0]), -(a[0][0]*a[2][1]-a[0][1]*a[2][0]), (a[0][0]*a[1][1]-a[0][1]*a[1][0])]]
    return [[cof[i][j]/d for j in range(3)] for i in range(3)]


def parse_poscar(path):
    lines = path.read_text().splitlines()
    scale = float(lines[1])
    latt = [[float(x)*scale for x in lines[i].split()] for i in (2, 3, 4)]
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    assert lines[7].strip().lower().startswith("d"), "expect Direct coords"
    elems, fracs = [], []
    row = 8
    for sp, n in zip(species, counts):
        for _ in range(n):
            fracs.append([float(x) for x in lines[row].split()[:3]])
            elems.append(sp)
            row += 1
    return latt, elems, fracs


def parse_wout_final(path):
    centres, spreads, in_final = [], [], False
    for line in path.read_text().splitlines():
        if "Final State" in line:
            in_final = True
            continue
        if in_final:
            if "WF centre and spread" not in line:
                if centres:
                    break
                continue
            body = line.split("(")[1]
            xyz, tail = body.split(")")
            centres.append([float(x) for x in xyz.split(",")])
            spreads.append(float(tail.split()[0]))
    return centres, spreads


def parse_eoscar(path):
    rows = []
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 6 and parts[0].isdigit():
            rows.append((int(parts[0]), parts[1], float(parts[3]), int(parts[5])))
    return rows  # (ion, el, ne_assigned, os_int)


def min_image_d2(frac_c, frac_a, latt):
    best = None
    df = [frac_c[k] - frac_a[k] for k in range(3)]
    df = [x - round(x) for x in df]
    for n0 in (-1, 0, 1):
        for n1 in (-1, 0, 1):
            for n2 in (-1, 0, 1):
                f = (df[0]+n0, df[1]+n1, df[2]+n2)
                cart = [sum(f[k]*latt[k][j] for k in range(3)) for j in range(3)]
                d2 = sum(c*c for c in cart)
                if best is None or d2 < best:
                    best = d2
    return best


class MinCostFlow:
    def __init__(self, n):
        self.n = n
        self.g = [[] for _ in range(n)]

    def add(self, u, v, cap, cost):
        self.g[u].append([v, cap, cost, len(self.g[v])])
        self.g[v].append([u, 0, -cost, len(self.g[u]) - 1])

    def solve(self, s, t):
        n, INF = self.n, float("inf")
        pot = [0.0]*n
        total = 0.0
        while True:
            dist = [INF]*n
            dist[s] = 0.0
            prev = [None]*n
            pq = [(0.0, s)]
            while pq:
                d, u = heapq.heappop(pq)
                if d > dist[u] + 1e-12:
                    continue
                for i, e in enumerate(self.g[u]):
                    v, cap, cost, _ = e
                    if cap <= 0:
                        continue
                    nd = d + cost + pot[u] - pot[v]
                    if nd < dist[v] - 1e-12:
                        dist[v] = nd
                        prev[v] = (u, i)
                        heapq.heappush(pq, (nd, v))
            if prev[t] is None:
                break
            for v in range(n):
                if dist[v] < INF:
                    pot[v] += dist[v]
            f = INF
            v = t
            while v != s:
                u, i = prev[v]
                f = min(f, self.g[u][i][1])
                v = u
            v = t
            while v != s:
                u, i = prev[v]
                self.g[u][i][1] -= f
                rev = self.g[u][i][3]
                self.g[v][rev][1] += f
                total += f * self.g[u][i][2]
                v = u
        return total


def main():
    case = Path(sys.argv[1])
    latt, elems, fracs = parse_poscar(case / "POSCAR")
    centres, spreads = parse_wout_final(case / "wannier90.wout")
    eoscar = parse_eoscar(case / "EOSCAR")
    binv = mat_inv3(latt)
    nw, na = len(centres), len(elems)
    print(f"case={case.name}  WFs={nw}  atoms={na}  electrons={2*nw}")

    # WF centres (cartesian) -> fractional
    wf_frac = [[sum(c[i]*binv[i][j] for i in range(3)) for j in range(3)] for c in centres]

    d2 = [[min_image_d2(wf_frac[w], fracs[a], latt) for a in range(na)] for w in range(nw)]

    # nodes: 0=source, 1..nw = WFs, nw+1..nw+na = atoms, nw+na+1 = sink
    beta = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
    chi_max = max(ENEG[e] for e in elems)
    S, T = 0, nw + na + 1
    mcf = MinCostFlow(nw + na + 2)
    for w in range(nw):
        mcf.add(S, 1 + w, 2, 0.0)               # each closed-shell WF carries 2 e
        for a in range(na):
            # +beta*(1 - chi/chi_max) is the electronegativity pull, shifted
            # to keep edge costs non-negative (constant offset per unit flow).
            c = d2[w][a] + beta * (1.0 - ENEG[elems[a]] / chi_max)
            mcf.add(1 + w, 1 + nw + a, 2, c)
    for a in range(na):
        mcf.add(1 + nw + a, T, CAPACITY[elems[a]], 0.0)
    cost = mcf.solve(S, T)

    ne = [0]*na
    for w in range(nw):
        for e in mcf.g[1 + w]:
            v, cap_left, c, rev = e
            if 1 + nw <= v <= nw + na and c >= 0:
                ne[v - 1 - nw] += 2 - cap_left

    assigned = sum(ne)
    print(f"assigned electrons: {assigned} / {2*nw}   transport cost = {cost:.3f} A^2")

    greedy = {ion: (el, ne_g, os_g) for ion, el, ne_g, os_g in eoscar}
    print(f"\n{'ion':>4} {'el':>3} {'ne_flow':>8} {'OS_flow':>8} {'ne_greedy':>10} {'OS_greedy':>10}")
    from collections import defaultdict
    elem_flow, elem_greedy = defaultdict(list), defaultdict(list)
    for a in range(na):
        el = elems[a]
        os_flow = ZVAL[el] - ne[a]
        g_el, g_ne, g_os = greedy.get(a + 1, (el, float("nan"), 0))
        print(f"{a+1:>4} {el:>3} {ne[a]:>8} {os_flow:>+8} {g_ne:>10.3f} {g_os:>+10}")
        elem_flow[el].append(os_flow)
        elem_greedy[el].append(g_os)
    print("\nper-element OS summary (min..max, mean):")
    for el in sorted(elem_flow):
        f, g = elem_flow[el], elem_greedy[el]
        print(f"  {el}: flow {min(f):+d}..{max(f):+d} mean {sum(f)/len(f):+.2f}   "
              f"greedy {min(g):+d}..{max(g):+d} mean {sum(g)/len(g):+.2f}")


if __name__ == "__main__":
    main()
