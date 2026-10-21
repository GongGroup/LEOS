#!/usr/bin/env python3
"""Per-spin wannier90 .win generator for the LEOS offline-w90 workflow.

A VASP export run (LEOS=.TRUE. + LWRITE_MMN_AMN=.TRUE. + LWANNIER90_RUN=.FALSE.)
writes per-spin wannier90.{1,2}.{mmn,amn,eig} but only ONE shared wannier90.win
carrying max(num_wann) over spins and the shared frozen window: the v3.1
per-spin ghost fix rewrites the spin .win only inside the library-mode path,
which the export run skips.  Running the minority channel offline with the
shared values reproduces the ghost bias this fix removed.

This tool replicates the VASP-side v3.1 logic from the export site and writes
wannier90.1.win / wannier90.2.win with per-spin num_wann and dis_froz_max.

Inputs read from the case directory:
  wannier90.win   shared template (required)
  PROCAR          eigenvalues + occupations, 8 decimals (preferred source)
  EIGENVAL        fallback source, 6 decimals
  OUTCAR          cross-check only: the "LEOS: per-spin num_wann" print.
                  (OUTCAR band energies carry 4 decimals - never used as the
                  eigenvalue source, F14.7 fidelity is impossible from them.)

Replicated algorithm (constants lifted from eos.F / mlwf.F, patch v3.1):
  num_wann(s)    = max over k of count(occ(b,k,s) > 1e-8)        [EOS_OCC_TOL]
  dis_froz_max(s): E_OCC_MAX  = max_k E(nw_s,   k, s)
                   E_NEXT_MIN = min_k E(nw_s+1, k, s)
                   gap > 0.05 [GAP_TOL] -> 0.5*(E_OCC_MAX + E_NEXT_MIN)
                   else                 -> min_k E(nw_s, k, s) - 1e-5
  line formats:   'num_wann = ' + I0        'dis_froz_max =' + F14.7
A spin at the shared max is NOT rewritten: MLWF_WANNIER90_RUN only calls the
rewrite under IF (NW /= MLWF%NUM_WANN), so that channel is a verbatim copy of
the shared wannier90.win (whose dis_froz_max is the cross-spin shared window).
Guards mirrored from MLWF_WIN_SET_NUM_WANN (rewrite is all-or-nothing per
spin): a 'begin projections' block, or no valid per-spin window while the
template has a dis_froz_max line, or a template num_wann that does not match
the computed max (explicit INCAR NUM_WANN = shared-upper-bound semantics)
-> the affected spin keeps the shared template verbatim.

The exported .amn keeps the shared max projection columns; library mode
feeds wannier90 the first nw_s columns of the in-memory A matrix, so for a
trimmed spin the offline equivalent is dropping the trailing projection
columns.  --trim-amn writes wannier90.<s>.amn.spin<s> next to the wins for
every trimmed spin (the original exported .amn is never touched).

Usage:
  offline_w90_prep.py CASE_DIR                      # wins into CASE_DIR
  offline_w90_prep.py CASE_DIR --out-dir D          # wins into D
  offline_w90_prep.py CASE_DIR --out-dir D --check CASE_DIR
                                # byte/numeric compare vs reference wins
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

# Constants lifted from the Fortran (eos.F / mlwf.F, patch v3.1).
EOS_OCC_TOL = 1e-8
GAP_TOL = 5e-2
FROZ_FALLBACK_SHIFT = 1e-5


class PrepError(RuntimeError):
    pass


# ---------------------------------------------------------------- parsers

def parse_procar(path: Path):
    """eig[s][k][b], occ[s][k][b] from PROCAR (8-decimal energies)."""
    band_re = re.compile(
        r"^band\s+\d+\s+#\s+energy\s+([-\d.]+)\s+#\s+occ\.\s+([-\d.]+)")
    eig, occ = [], []
    with path.open(errors="replace") as fh:
        for line in fh:
            if line.startswith("# of k-points:"):
                eig.append([])          # new spin section
                occ.append([])
            elif line.startswith(" k-point"):
                eig[-1].append([])
                occ[-1].append([])
            else:
                m = band_re.match(line)
                if m:
                    eig[-1][-1].append(float(m.group(1)))
                    occ[-1][-1].append(float(m.group(2)))
    if not eig or not eig[0] or not eig[0][0]:
        raise PrepError(f"no bands parsed from {path}")
    return eig, occ


def parse_eigenval(path: Path):
    """eig[s][k][b], occ[s][k][b] from EIGENVAL (6-decimal energies)."""
    lines = path.read_text(errors="replace").splitlines()
    try:
        ispin = int(lines[0].split()[3])
        _, nk, nb = (int(x) for x in lines[5].split()[:3])
    except (IndexError, ValueError) as exc:
        raise PrepError(f"unrecognised EIGENVAL header in {path}") from exc
    eig = [[[0.0] * nb for _ in range(nk)] for _ in range(ispin)]
    occ = [[[0.0] * nb for _ in range(nk)] for _ in range(ispin)]
    i = 6
    for k in range(nk):
        while i < len(lines) and not lines[i].strip():
            i += 1
        i += 1                                    # k-vector + weight line
        for b in range(nb):
            t = lines[i].split()
            i += 1
            if ispin == 2:
                eig[0][k][b], eig[1][k][b] = float(t[1]), float(t[2])
                occ[0][k][b], occ[1][k][b] = float(t[3]), float(t[4])
            else:
                eig[0][k][b] = float(t[1])
                occ[0][k][b] = float(t[2])
    return eig, occ


def read_eigenvalues(case: Path):
    """Best available eigenvalue source: (source_name, eig, occ)."""
    procar = case / "PROCAR"
    if procar.is_file():
        try:
            eig, occ = parse_procar(procar)
            return "PROCAR", eig, occ
        except PrepError as exc:
            print(f"WARNING: {exc}; falling back to EIGENVAL", file=sys.stderr)
    eigenval = case / "EIGENVAL"
    if eigenval.is_file():
        eig, occ = parse_eigenval(eigenval)
        return "EIGENVAL", eig, occ
    raise PrepError(f"neither PROCAR nor EIGENVAL found in {case}")


def outcar_leos_num_wann(path: Path):
    """(shared, [per-spin]) from the LEOS OUTCAR prints, or (None, None).

    'per-spin num_wann = A (up) / B (down)' is only printed when the two
    channels differ; with only the 'auto set' line both spins equal shared."""
    shared, spins, explicit = None, None, False
    if not path.is_file():
        return None, None, False
    with path.open(errors="replace") as fh:
        for line in fh:
            m = re.search(r"LEOS: auto set NUM_WANN = (\d+)", line)
            if m:
                shared = int(m.group(1))
            m = re.search(r"LEOS: per-spin num_wann = (\d+) \(up\) / (\d+) \(down\)",
                          line)
            if m:
                spins = [int(m.group(1)), int(m.group(2))]
            if "explicit NUM_WANN is the shared upper bound" in line:
                explicit = True
    if shared is not None and spins is None:
        spins = [shared, shared]
    return shared, spins, explicit


# ---------------------------------------------------------------- algorithm

def per_spin_params(eig, occ, occ_tol=EOS_OCC_TOL, gap_tol=GAP_TOL):
    """[(num_wann, dis_froz_max-or-None), ...] per spin channel.

    Mirrors EOS_PREPARE_MLWF (band count) + the per-spin branch of
    ESTIMATE_AUTO_WINDOW (frozen window).  dis_froz_max is None when the
    channel occupies every band (no window can be defined)."""
    out = []
    for s, (es, os_) in enumerate(zip(eig, occ)):
        nb = len(es[0])
        nw = max(sum(1 for f in kocc if f > occ_tol) for kocc in os_)
        if nw <= 0:
            raise PrepError(f"spin {s + 1}: no occupied bands found")
        if nw >= nb:
            out.append((nw, None))
            continue
        emax_nw = max(ek[nw - 1] for ek in es)      # ELIM(2, nw)
        emin_nw = min(ek[nw - 1] for ek in es)      # ELIM(1, nw)
        emin_up = min(ek[nw] for ek in es)          # ELIM(1, nw+1)
        gap = emin_up - emax_nw
        if gap > gap_tol:
            froz = 0.5 * (emax_nw + emin_up)
        else:
            froz = emin_nw - FROZ_FALLBACK_SHIFT
        out.append((nw, froz))
    return out


def fmt_froz(froz: float) -> str:
    """Fortran WRITE('(A,F14.7)') 'dis_froz_max =', froz."""
    return "dis_froz_max =" + f"{froz:14.7f}"


def fmt_num_wann(nw: int) -> str:
    """Fortran WRITE('(A,I0)') 'num_wann = ', nw."""
    return f"num_wann = {nw}"


def template_num_wann(lines) -> int | None:
    for ln in lines:
        if ln.lstrip().startswith("num_wann"):
            try:
                return int(ln.split("=", 1)[1].strip())
            except (IndexError, ValueError):
                return None
    return None


def render_win(lines, nw, froz):
    """Per-spin win lines, or None when the Fortran guards refuse a rewrite
    (caller keeps the shared template verbatim, exactly like the fallback in
    MLWF_WIN_SET_NUM_WANN)."""
    has_froz = any(ln.lstrip().startswith("dis_froz_max") for ln in lines)
    if any("begin projections" in ln for ln in lines):
        return None                       # warm-start projections block
    if has_froz and froz is None:
        return None                       # window present, no per-spin value
    out = []
    for ln in lines:
        stripped = ln.lstrip()
        if stripped.startswith("num_wann"):
            out.append(fmt_num_wann(nw))
        elif stripped.startswith("dis_froz_max") and froz is not None:
            out.append(fmt_froz(froz))
        else:
            out.append(ln.rstrip())       # Fortran writes TRIM(line)
    return out


def trim_amn(src: Path, dst: Path, nw: int) -> None:
    """Keep the first nw projection columns of an exported .amn.

    Library mode passes wannier90 the leading nw columns of the shared-size
    A matrix, so dropping rows with iw > nw is the exact offline equivalent."""
    with src.open(errors="replace") as fin, dst.open("w") as fout:
        title = fin.readline()
        header = fin.readline().split()
        nb, nk, nw_shared = (int(x) for x in header[:3])
        if nw > nw_shared:
            raise PrepError(f"{src}: cannot trim to {nw} > shared {nw_shared}")
        fout.write(title)
        fout.write(f"{nb:12d}{nk:12d}{nw:12d}\n")
        for line in fin:
            t = line.split(None, 2)
            if len(t) >= 2 and int(t[1]) <= nw:
                fout.write(line)


# ---------------------------------------------------------------- driver

def generate(case: Path, out_dir: Path | None = None, overwrite: bool = False,
             do_trim_amn: bool = False, crosscheck: bool = True) -> dict:
    """Generate per-spin win files for an export case directory.

    Returns a summary dict; raises PrepError on any inconsistency."""
    case = Path(case)
    out_dir = Path(out_dir) if out_dir else case
    template = case / "wannier90.win"
    if not template.is_file():
        raise PrepError(f"missing shared template {template}")
    tmpl_lines = template.read_text(errors="replace").splitlines()

    source, eig, occ = read_eigenvalues(case)
    ispin = len(eig)
    nk, nb = len(eig[0]), len(eig[0][0])
    if ispin == 1:
        return {"ispin": 1, "source": source, "spins": [],
                "note": "ISPIN=1: shared wannier90.win is already complete"}

    params = per_spin_params(eig, occ)
    shared_nw = max(nw for nw, _ in params)

    tmpl_nw = template_num_wann(tmpl_lines)
    shared_semantics = tmpl_nw is not None and tmpl_nw != shared_nw

    oc_shared, oc_spins, oc_explicit = outcar_leos_num_wann(case / "OUTCAR")
    if crosscheck and oc_spins is not None and not shared_semantics:
        mine = [nw for nw, _ in params]
        if mine != oc_spins[:ispin]:
            raise PrepError(
                f"per-spin num_wann mismatch: computed {mine} from {source}, "
                f"OUTCAR prints {oc_spins[:ispin]} - refusing to write")
    if oc_explicit:
        shared_semantics = True

    summary = {"ispin": ispin, "nk": nk, "nb": nb, "source": source,
               "shared_num_wann": shared_nw, "template_num_wann": tmpl_nw,
               "outcar_num_wann": oc_spins, "spins": []}

    tmpl_text = template.read_text(errors="replace")

    out_dir.mkdir(parents=True, exist_ok=True)
    for s, (nw, froz) in enumerate(params, start=1):
        target = out_dir / f"wannier90.{s}.win"
        if target.exists() and not overwrite:
            raise PrepError(f"{target} exists - use --overwrite or --out-dir")
        if nw == shared_nw and not shared_semantics:
            # VASP only rewrites under IF (NW /= NUM_WANN): the max channel
            # is a verbatim COPYFILE of the shared win (cross-spin window).
            target.write_text(tmpl_text)
            summary["spins"].append({"spin": s, "num_wann": nw,
                                     "dis_froz_max": None, "win": str(target),
                                     "kept_shared": "at shared max "
                                     "(VASP copies wannier90.win unchanged)"})
            continue
        if shared_semantics:
            lines, kept = None, "template num_wann mismatch (explicit NUM_WANN?)"
        else:
            lines = render_win(tmpl_lines, nw, froz)
            kept = "rewrite guard (projections block / missing window)"
        if lines is None:
            target.write_text(tmpl_text)
            summary["spins"].append({"spin": s, "num_wann": tmpl_nw,
                                     "dis_froz_max": None, "win": str(target),
                                     "kept_shared": kept})
            continue
        target.write_text("\n".join(lines) + "\n")
        rec = {"spin": s, "num_wann": nw, "dis_froz_max": froz,
               "win": str(target), "kept_shared": None}
        if do_trim_amn and nw < shared_nw:
            src_amn = case / f"wannier90.{s}.amn"
            if src_amn.is_file():
                dst_amn = out_dir / f"wannier90.{s}.amn.spin{s}"
                trim_amn(src_amn, dst_amn, nw)
                rec["trimmed_amn"] = str(dst_amn)
            else:
                rec["trimmed_amn"] = f"SKIPPED ({src_amn} missing)"
        summary["spins"].append(rec)
    return summary


def check(summary: dict, ref_dir: Path, tol: float) -> int:
    """Compare generated wins against reference per-spin wins. Returns rc."""
    rc = 0
    for rec in summary["spins"]:
        s = rec["spin"]
        gen = Path(rec["win"])
        ref = Path(ref_dir) / f"wannier90.{s}.win"
        if not ref.is_file():
            print(f"[check] spin {s}: reference {ref} missing - skipped")
            continue
        gtext, rtext = gen.read_text(), ref.read_text()
        byte_ok = gtext == rtext
        gnw = template_num_wann(gtext.splitlines())
        rnw = template_num_wann(rtext.splitlines())
        gfz = rfz = None
        for ln in gtext.splitlines():
            if ln.lstrip().startswith("dis_froz_max"):
                gfz = float(ln.split("=", 1)[1])
        for ln in rtext.splitlines():
            if ln.lstrip().startswith("dis_froz_max"):
                rfz = float(ln.split("=", 1)[1])
        dfz = (abs(gfz - rfz) if gfz is not None and rfz is not None
               else (0.0 if gfz == rfz else float("inf")))
        ok = byte_ok or (gnw == rnw and dfz <= tol)
        print(f"[check] spin {s}: byte-identical={byte_ok} num_wann {gnw} vs "
              f"{rnw} | dis_froz_max delta={dfz:.3e} (tol {tol:g}) -> "
              f"{'OK' if ok else 'MISMATCH'}")
        if not ok:
            rc = 1
    return rc


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case", type=Path, help="VASP export case directory")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="write wins here (default: the case directory)")
    ap.add_argument("--overwrite", action="store_true",
                    help="allow replacing existing wannier90.<s>.win")
    ap.add_argument("--trim-amn", action="store_true",
                    help="also write column-trimmed .amn for trimmed spins")
    ap.add_argument("--no-crosscheck", action="store_true",
                    help="skip the OUTCAR per-spin num_wann cross-check")
    ap.add_argument("--check", type=Path, metavar="REF_DIR",
                    help="compare generated wins against REF_DIR/wannier90.<s>.win")
    ap.add_argument("--check-tol", type=float, default=0.0,
                    help="allowed |delta dis_froz_max| in --check (default 0)")
    a = ap.parse_args()
    try:
        summary = generate(a.case, a.out_dir, a.overwrite, a.trim_amn,
                           crosscheck=not a.no_crosscheck)
    except PrepError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(f"source={summary['source']} ispin={summary['ispin']}", end="")
    if summary["ispin"] == 1:
        print(f"  ({summary['note']})")
        return 0
    print(f" nk={summary['nk']} nb={summary['nb']} "
          f"shared_num_wann={summary['shared_num_wann']}")
    for rec in summary["spins"]:
        froz = rec["dis_froz_max"]
        froz_s = f"{froz:.7f}" if froz is not None else "kept-shared"
        extra = f"  [{rec['kept_shared']}]" if rec["kept_shared"] else ""
        trim = f"  amn->{rec['trimmed_amn']}" if rec.get("trimmed_amn") else ""
        print(f"  spin {rec['spin']}: num_wann={rec['num_wann']} "
              f"dis_froz_max={froz_s} -> {rec['win']}{extra}{trim}")
    if a.check:
        return check(summary, a.check, a.check_tol)
    return 0


if __name__ == "__main__":
    sys.exit(main())
