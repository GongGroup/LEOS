"""Clean-bench plugin: transport with the electronegativity term off (beta=0).

Reference point for table-free candidates: hosts are chosen by WF-centre
geometry under the capacity bounds. The cost builder still evaluates the
electronegativity table (multiplied by zero), so the run reports it as used.
"""
from pyleos.transport import assign_case

META = {"name": "geometry-only transport", "params": {"beta": 0.0}}


def assign(case):
    return assign_case(case, beta=0.0).os_int
