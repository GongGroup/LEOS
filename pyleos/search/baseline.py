"""Clean-bench plugin: the shipped transport assignment (beta=10, Allen, tier 1).

Template for search candidates: a plugin defines assign(case) returning one
integer oxidation state per atom in POSCAR order, and may declare META with
its free parameters. See script/clean_bench.py for what a plugin may use.
"""
from pyleos.transport import assign_case

META = {"name": "transport baseline", "params": {"beta": 10.0}}


def assign(case):
    return assign_case(case, beta=10.0).os_int
