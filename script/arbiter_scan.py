#!/usr/bin/env python3
"""Arbiter TAU x PMARGIN sensitivity map (v3.1 task 3, 2026-07-10).

Boss-approved scope: produce the table ONLY - no code change is made here.
pyleos/transport.py is untouched; the two module constants are overridden
in-process per grid point. The decision (move defaults / forbid H flips /
leave as is) is taken by the boss after reading the table (HANDOFF_V31).

Per grid point (TAU, PMARGIN), tier2="arbiter" is run over
  (a) every #WP-carrying benchmark case (same set as full_bench --only-wp),
      scored with the same dev<=0.25 criterion, and
  (b) the six legacy_rerun directories, scored per-atom on hydrogen against
      the verified patterns (RESULTS 3.14): all H +1, except ce_fs/ce_ts
      where exactly one hydride H(-1) is expected.
Baseline column: tier-1 pure geometry (tier2=linear, gamma=0), the
production default. Note the baseline itself carries 1 legacy H error
(one vircef H at 0 instead of +1) - a known ghost-bias residue slated
for the v3.1 num_wann fix (RESULTS 3.14/5.2), not an arbiter artefact.

A point is SAFE when bench passes >= baseline AND it introduces no NEW
legacy H errors (herr <= baseline herr).

Usage: python3 script/arbiter_scan.py [--out data/v3/arbiter_scan_20260710.tsv]
       [--root DIR]   # alternate data tree mirroring the repo layout
                      # (data/v3/m4_wave*/..., data/v3/legacy_rerun/...);
                      # default: the repo itself, byte-identical behaviour
"""
from __future__ import annotations

import argparse
import functools
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "script"))

import full_bench                                  # noqa: E402
import pyleos.transport as T                       # noqa: E402
from pyleos.parsers import load_case               # noqa: E402
from pyleos.transport import assign_case           # noqa: E402

print = functools.partial(print, flush=True)

BETA = 10.0          # production value (cli.py / full_bench defaults)
TAUS = [0.25, 0.5, 1.0, 1.5, 2.0, 3.0]
PMARGINS = [0.02, 0.05, 0.10, 0.20, 0.30]

# Verified hydrogen patterns per legacy case (RESULTS 3.14): value -> count.
LEGACY_H = {
    "mumnox": None,   # None = every H must be +1 (count read from the case)
    "vircef": None,
    "lasweh": None,
    "ce_is":  None,
    "ce_fs":  "one_hydride",   # exactly one H at -1, the rest +1
    "ce_ts":  "one_hydride",
}


def load_bench():
    """(name, case, ref_map) for every #WP benchmark case."""
    out = []
    for _batch, name, d, ref_map in full_bench.collect():
        case = load_case(d)
        if case.eoscar and case.eoscar.wp_detail:
            out.append((name, case, ref_map))
    return out


def load_legacy(root=ROOT):
    """(name, case, expected_sorted_H) for the six legacy directories."""
    out = []
    for name, mode in LEGACY_H.items():
        case = load_case(root / "data/v3/legacy_rerun" / name)
        n_h = sum(1 for el in case.structure.elements if el == "H")
        expected = sorted([-1] + [1] * (n_h - 1)) if mode == "one_hydride" \
            else [1] * n_h
        out.append((name, case, expected))
    return out


def bench_dev(case, ref_map, res):
    means = full_bench.element_means(case.structure.elements, res.os_int)
    return max((abs(means[el] - t) for el, t in ref_map.items() if el in means),
               default=None)


def h_errors(case, res, expected):
    got = sorted(v for el, v in zip(case.structure.elements, res.os_int)
                 if el == "H")
    return sum(1 for a, b in zip(got, expected) if a != b) + \
        abs(len(got) - len(expected))


def run_point(bench, legacy, tier2, tau=None, pmargin=None,
              baseline_os=None):
    """One grid point (or the baseline when tier2='linear')."""
    if tau is not None:
        T.ARBITER_TAU = tau
    if pmargin is not None:
        T.ARBITER_PMARGIN = pmargin
    ok, names_ok = 0, set()
    for name, case, ref_map in bench:
        res = assign_case(case, beta=BETA, tier2=tier2)
        dev = bench_dev(case, ref_map, res)
        if dev is not None and dev <= 0.25:
            ok += 1
            names_ok.add(name)
    herr, flips, legacy_os = 0, 0, {}
    for name, case, expected in legacy:
        res = assign_case(case, beta=BETA, tier2=tier2)
        herr += h_errors(case, res, expected)
        legacy_os[name] = list(res.os_int)
        if baseline_os is not None:
            flips += sum(1 for el, a, b in zip(case.structure.elements,
                                               res.os_int, baseline_os[name])
                         if el != "H" and a != b)
    return ok, names_ok, herr, flips, legacy_os


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path,
                    default=ROOT / "data/v3/arbiter_scan_20260710.tsv")
    ap.add_argument("--root", type=Path, default=None,
                    help="alternate data tree root mirroring the repo layout "
                         "(default: the repo root; collection rules and "
                         "scoring are unchanged)")
    args = ap.parse_args()
    root = args.root.resolve() if args.root else ROOT
    full_bench.ROOT = root      # collect() resolves case dirs against this

    t0 = time.time()
    bench, legacy = load_bench(), load_legacy(root)
    print(f"loaded {len(bench)} #WP bench cases + {len(legacy)} legacy cases "
          f"({time.time()-t0:.1f}s); grid {len(TAUS)}x{len(PMARGINS)} "
          f"+ baseline, beta={BETA}")

    t0 = time.time()
    base_ok, base_names, base_herr, _, base_os = run_point(
        bench, legacy, tier2="linear")
    print(f"[baseline tier-1] bench {base_ok}/{len(bench)} "
          f"legacy_H_err={base_herr} ({time.time()-t0:.1f}s)")

    rows, n_pt = [], len(TAUS) * len(PMARGINS)
    for i, tau in enumerate(TAUS):
        for j, pm in enumerate(PMARGINS):
            t1 = time.time()
            ok, names, herr, flips, _ = run_point(
                bench, legacy, tier2="arbiter", tau=tau, pmargin=pm,
                baseline_os=base_os)
            gained = ",".join(sorted(names - base_names)) or "-"
            lost = ",".join(sorted(base_names - names)) or "-"
            safe = "SAFE" if (ok >= base_ok and herr <= base_herr) else "unsafe"
            rows.append((tau, pm, ok, herr, flips, safe, gained, lost))
            k = i * len(PMARGINS) + j + 1
            print(f"[{k:2d}/{n_pt}] tau={tau:<5} pm={pm:<5} "
                  f"bench {ok}/{len(bench)} H_err={herr} "
                  f"nonH_flips={flips} {safe} (+{gained} -{lost}) "
                  f"({time.time()-t1:.1f}s)")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    src = f"; root={root}" if root != ROOT else ""
    with args.out.open("w") as fh:
        fh.write(f"# arbiter sensitivity scan {time.strftime('%Y-%m-%d')}; "
                 f"beta={BETA}; baseline tier-1: bench {base_ok}/{len(bench)}, "
                 f"legacy_H_err={base_herr}{src}\n")
        fh.write("tau\tpmargin\tbench_ok\tlegacy_H_err\tnonH_flips\t"
                 "verdict\tgained\tlost\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    print(f"\n== summary grid (bench_ok/H_err), baseline {base_ok} ==")
    print("tau\\pm " + " ".join(f"{pm:>7}" for pm in PMARGINS))
    for tau in TAUS:
        cells = [next(r for r in rows if r[0] == tau and r[1] == pm)
                 for pm in PMARGINS]
        print(f"{tau:<6} " + " ".join(
            f"{c[2]:>3}/{c[3]:<3}" for c in cells))
    n_safe = sum(1 for r in rows if r[5] == "SAFE")
    print(f"\n{n_safe}/{len(rows)} points SAFE "
          f"(bench >= {base_ok} and legacy H errors <= baseline "
          f"{base_herr}) -> {args.out}")


if __name__ == "__main__":
    main()
