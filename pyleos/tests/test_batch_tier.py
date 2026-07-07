#!/usr/bin/env python3
"""Unit tests for the tiered SCDM screening logic in batch.py.

The trigger rules must capture the six known SCDM-only mismatches from the
38-case bench: WCl4/5/6 + MoCl2 via the heavy-TM element rule (pre-escalated
to full MLWF at stage time), GaSe/CSe2-type cases via runtime quality flags
(stray spread / dNe drift) picked up by escalate."""
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from pyleos.batch import (HEAVY_TM, eoscar_quality, poscar_elements,   # noqa: E402
                          escalate, load_state, save_state)

POSCAR_WCL6 = """WCl6
1.0
6 0 0
0 6 0
0 0 6
W Cl
1 6
Direct
0 0 0
""" + "0.2 0 0\n" * 6

POSCAR_CO2 = """CO2
1.0
6 0 0
0 6 0
0 0 6
C O
1 2
Direct
0 0 0
0.1 0 0
0.9 0 0
"""

EOSCAR_CLEAN = """# LEOS
# NELECT(Wannier) =      16.000000
# |dNe| =   1.2E-14
#WF  1     1  2.000000   0.0 0.0 0.0   0.8519     1  2.0000     2  0.0000
    1  C    4.00      0.000000     4.000       4     0.000
    2  O    6.00      8.000000    -2.000      -2     0.000
"""

EOSCAR_STRAY = EOSCAR_CLEAN.replace(
    "#WF  1     1  2.000000   0.0 0.0 0.0   0.8519",
    "#WF  1     1  2.000000   0.0 0.0 0.0   5.4000")

EOSCAR_DNE = EOSCAR_CLEAN.replace("1.2E-14", "3.1E-02")


class TierTests(unittest.TestCase):
    def test_heavy_tm_rule_catches_4d5d_mismatches(self):
        # the four element-rule cases from the 38-case bench
        for els in ({"W", "Cl"}, {"Mo", "Cl"}):
            self.assertTrue(els & HEAVY_TM)
        self.assertFalse({"C", "Se"} & HEAVY_TM)   # CSe2 must rely on flags
        self.assertFalse({"Ga", "Se"} & HEAVY_TM)  # GaSe likewise

    def test_poscar_elements(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "POSCAR"
            p.write_text(POSCAR_WCL6)
            self.assertEqual(poscar_elements(p), {"W", "Cl"})

    def test_eoscar_quality_flags(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "EOSCAR"
            p.write_text(EOSCAR_CLEAN)
            self.assertEqual(eoscar_quality(p), [])
            p.write_text(EOSCAR_STRAY)
            self.assertEqual(eoscar_quality(p), ["stray:1"])
            p.write_text(EOSCAR_DNE)
            self.assertEqual(eoscar_quality(p), ["dNe"])

    def test_escalate_rewrites_and_requeues(self):
        with tempfile.TemporaryDirectory() as td:
            dest = Path(td)
            d = dest / "gase"
            d.mkdir()
            (d / "EOSCAR").write_text(EOSCAR_STRAY)
            (d / "INCAR").write_text("LEOS = .TRUE.\nEOS_SCDM_ONLY = .TRUE.\n")
            save_state(dest, {"gase": {"status": "OK", "jobid": "1",
                                       "retries": 0, "tier": "scdm"}})
            escalate(dest)
            state = load_state(dest)
            self.assertEqual(state["gase"]["status"], "RETRY")
            self.assertEqual(state["gase"]["tier"], "full")
            self.assertIn("EOS_SCDM_ONLY = .FALSE.",
                          (d / "INCAR").read_text())
            self.assertTrue((d / "EOSCAR.scdm").is_file())
            self.assertFalse((d / "EOSCAR").is_file())

    def test_escalate_keeps_clean_screens(self):
        with tempfile.TemporaryDirectory() as td:
            dest = Path(td)
            d = dest / "co2"
            d.mkdir()
            (d / "EOSCAR").write_text(EOSCAR_CLEAN)
            (d / "INCAR").write_text("EOS_SCDM_ONLY = .TRUE.\n")
            save_state(dest, {"co2": {"status": "OK", "jobid": "1",
                                      "retries": 0, "tier": "scdm"}})
            escalate(dest)
            state = load_state(dest)
            self.assertEqual(state["co2"]["status"], "OK")
            self.assertTrue((d / "EOSCAR").is_file())


if __name__ == "__main__":
    unittest.main(verbosity=1)
