"""cases/INDEX.tsv, cases/STRUCTURES.tsv and the run directories agree.

Run: python3 pyleos/tests/test_registry.py
"""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from pyleos.cases import check_registry, min_distance          # noqa: E402
from pyleos.parsers import Structure                          # noqa: E402


class RegistryTests(unittest.TestCase):
    def test_registry_consistent(self):
        problems = check_registry()
        self.assertEqual(problems, [], "\n".join(problems[:20]))

    def test_min_distance_sees_periodic_images(self):
        st = Structure([[3.0, 0, 0], [0, 3.0, 0], [0, 0, 3.0]], ["Li", "Li"],
                       [[0.05, 0, 0], [0.95, 0, 0]])
        self.assertAlmostEqual(min_distance(st), 0.3, places=6)


if __name__ == "__main__":
    unittest.main()
