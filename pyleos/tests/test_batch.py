#!/usr/bin/env python3
"""Unit tests for the batch driver's failure handling (issue #5).

Covers the three confirmed defects: a FAIL_W90 retry must really change the
INCAR (NBANDS_SCALE resolved to a concrete NBANDS, SCDM turned off), a
diagnostics-only EOSCAR carrying the LEOS_ABORTED marker must triage as
FAIL_LEOS instead of OK, and a failed squeue query must return None and
poison no case state (no re-triage, no submission)."""
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from pyleos.batch import (load_state, retry, save_state, squeue_ids,  # noqa: E402
                          status, submit, triage)

REPO = Path(__file__).resolve().parent.parent.parent
ABORTED_CASE = REPO / "data" / "v3" / "regress_v31" / "mumnox_abort2"


def _make_case(dest: Path, name: str, incar: str, status_: str) -> Path:
    d = dest / name
    d.mkdir(parents=True)
    (d / "INCAR").write_text(incar)
    state = load_state(dest)
    state[name] = {"status": status_, "jobid": None, "retries": 0,
                   "tier": "full"}
    save_state(dest, state)
    return d


class BatchFailureHandlingTests(unittest.TestCase):
    def test_w90_retry_really_changes_incar(self):
        with tempfile.TemporaryDirectory() as td:
            dest = Path(td)
            # NBANDS in the INCAR: scaled in place
            d1 = _make_case(dest, "nb_incar",
                            "LEOS = .TRUE.\nNBANDS = 100\n", "FAIL_W90")
            # NBANDS only in OUTCAR: picked up from there
            d2 = _make_case(dest, "nb_outcar", "LEOS = .TRUE.\n", "FAIL_W90")
            (d2 / "OUTCAR").write_text(
                "   k-points  NKPTS =   1   number of bands    NBANDS=     64\n")
            # NBANDS nowhere: fall back to a tag that still changes the run
            d3 = _make_case(dest, "nb_none", "LEOS = .TRUE.\n", "FAIL_W90")
            before = {d: (d / "INCAR").read_text() for d in (d1, d2, d3)}

            retry(dest)

            for d in (d1, d2, d3):
                after = (d / "INCAR").read_text()
                self.assertNotEqual(before[d], after,
                                    f"{d.name}: retry left the INCAR unchanged")
                self.assertIn("EOS_SCDM_ONLY = .FALSE.", after)
            self.assertIn("NBANDS = 130", (d1 / "INCAR").read_text())  # 100 x 1.3
            self.assertIn("NBANDS = 83", (d2 / "INCAR").read_text())   # 64 x 1.3
            self.assertNotIn("NBANDS", (d3 / "INCAR").read_text())
            state = load_state(dest)
            for name in ("nb_incar", "nb_outcar", "nb_none"):
                self.assertEqual(state[name]["status"], "RETRY")
                self.assertEqual(state[name]["retries"], 1)

    def test_triage_flags_aborted_eoscar_as_fail_leos(self):
        # real diagnostics-only EOSCAR from the NELECT-abort regression case
        self.assertTrue((ABORTED_CASE / "EOSCAR").is_file(),
                        f"fixture missing: {ABORTED_CASE}")
        self.assertEqual(triage(ABORTED_CASE), "FAIL_LEOS")
        # and a normal EOSCAR must still triage as OK
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            (d / "EOSCAR").write_text(
                "# LEOS\n# |dNe| =   1.2E-14\n    1  C    4.00\n")
            self.assertEqual(triage(d), "OK")

    def test_squeue_failure_returns_none_and_changes_nothing(self):
        fail = mock.Mock(returncode=1, stdout="", stderr="slurm_load_jobs error")
        with mock.patch("pyleos.batch.subprocess.run", return_value=fail):
            self.assertIsNone(squeue_ids())
        with mock.patch("pyleos.batch.subprocess.run",
                        side_effect=FileNotFoundError):
            self.assertIsNone(squeue_ids())

        with tempfile.TemporaryDirectory() as td:
            dest = Path(td)
            _make_case(dest, "running", "LEOS = .TRUE.\n", "SUBMITTED")
            _make_case(dest, "waiting", "LEOS = .TRUE.\n", "STAGED")
            state_before = (dest / "batch_state.json").read_text()
            with mock.patch("pyleos.batch.subprocess.run",
                            return_value=fail) as run:
                status(dest)     # must not re-triage "running" as failed
                submit(dest, 8)  # must not sbatch "waiting"
                sbatch_calls = [c for c in run.call_args_list
                                if c.args[0][0] == "sbatch"]
                self.assertEqual(sbatch_calls, [])
            self.assertEqual((dest / "batch_state.json").read_text(),
                             state_before)
            state = load_state(dest)
            self.assertEqual(state["running"]["status"], "SUBMITTED")
            self.assertEqual(state["waiting"]["status"], "STAGED")


if __name__ == "__main__":
    unittest.main(verbosity=1)
