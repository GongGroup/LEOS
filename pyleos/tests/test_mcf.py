"""Unit tests for the min-cost-flow solver and capacity bounds. Stdlib only.

Run: python3 pyleos/tests/test_mcf.py
"""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from pyleos.capacity import capacity_bounds            # noqa: E402
from pyleos.transport import MinCostFlow, BIG_NEG      # noqa: E402


class TestMinCostFlow(unittest.TestCase):
    def test_simple_assignment(self):
        # 2 WFs x 2 atoms; WF0 cheap to atom0, WF1 cheap to atom1.
        # nodes: 0 s, 1-2 wf, 3-4 atoms, 5 t
        mcf = MinCostFlow(6)
        mcf.add(0, 1, 2, 0); mcf.add(0, 2, 2, 0)
        mcf.add(1, 3, 2, 1.0); mcf.add(1, 4, 2, 9.0)
        mcf.add(2, 3, 2, 9.0); mcf.add(2, 4, 2, 1.0)
        mcf.add(3, 5, 4, 0); mcf.add(4, 5, 4, 0)
        cost = mcf.solve(0, 5)
        self.assertAlmostEqual(cost, 4.0)               # 2e*1.0 + 2e*1.0

    def test_capacity_forces_split(self):
        # One atom capped at 1 electron: second electron must go elsewhere.
        mcf = MinCostFlow(5)                            # s, wf, a0, a1, t
        mcf.add(0, 1, 2, 0)
        mcf.add(1, 2, 2, 1.0)
        mcf.add(1, 3, 2, 5.0)
        mcf.add(2, 4, 1, 0)                             # cap 1
        mcf.add(3, 4, 8, 0)
        cost = mcf.solve(0, 4)
        self.assertAlmostEqual(cost, 6.0)               # 1*1 + 1*5

    def test_nmin_reward_edge(self):
        # Mandatory n_min on the far atom wins despite distance cost.
        mcf = MinCostFlow(5)
        mcf.add(0, 1, 2, 0)
        mcf.add(1, 2, 2, 1.0)                           # near atom
        mcf.add(1, 3, 2, 5.0)                           # far atom, n_min=1
        mcf.add(2, 4, 8, 0)
        mcf.add(3, 4, 1, -BIG_NEG)
        mcf.add(3, 4, 7, 0)
        mcf.solve(0, 4)
        sent_far = 2 - mcf.g[1][1][1]                   # edge 1->3 remaining cap
        self.assertEqual(sent_far, 1)

    def test_integer_flow(self):
        # Symmetric homonuclear bond: 2e split 1+1 or 2+0, both cost-optimal;
        # flow must still be integral and conserve electrons.
        mcf = MinCostFlow(5)
        mcf.add(0, 1, 2, 0)
        mcf.add(1, 2, 2, 1.0)
        mcf.add(1, 3, 2, 1.0)
        mcf.add(2, 4, 8, 0); mcf.add(3, 4, 8, 0)
        cost = mcf.solve(0, 4)
        self.assertAlmostEqual(cost, 2.0)
        total = (2 - mcf.g[1][0][1]) + (2 - mcf.g[1][1][1])
        self.assertEqual(total, 2)


class TestCapacity(unittest.TestCase):
    def test_oxygen(self):
        self.assertEqual(capacity_bounds("O", 6.0), (4, 8))   # OS +2..-2

    def test_semicore_na(self):
        # Na_pv ZVAL=7 (3s2 3p6 4s1 minus... 2s2p semicore in valence):
        # OS in [-1, +1] -> n in [6, 8]; semicore stays filled via n_min.
        self.assertEqual(capacity_bounds("Na", 7.0), (6, 8))

    def test_sulfur_wide(self):
        n_min, n_max = capacity_bounds("S", 6.0)
        self.assertEqual((n_min, n_max), (0, 8))              # +6..-2


if __name__ == "__main__":
    unittest.main(verbosity=2)
