"""
STEP 4 (canonical-order variant). Reconstruction of ALL optimal reference polygons.
Same task as best-fence-estimate-hexa-step4-count.py, but the forward pass of the dynamic program
fixes the order of identical pieces (port of best-fence-estimate-hexa-step3-count2.py, see
best_fence_estimate_hexa_dp_canon.py).  The DAG then has ONE path per oriented polygon
(instead of one path per ordering of identical pieces), so enumeration is feasible.
=================================
Usage: pypy3 best-fence-step4.py <n_max> <goal> [--only SET_NUMBER] [--max-solutions N]
XX (arguments as in best-fence-estimate-hexa-step4-count.py; needs the same Step 3 result file.)
Note: the forward pass is re-run here with the canonical-order tables, which are NOT the same as
those of Step 3 (no sum-normalisation, extra last-octant component).  The optimum must agree.
"""
import sys
import argparse

import best_fence_estimate_hexa_dp as dp
import best_fence_estimate_hexa_dp_canon as canon

_parser = argparse.ArgumentParser(
    description="Step 4: count and generate all optimal reference polygons.")
_parser.add_argument("n_max", type=int, help="size of the polyominoes (6 for hexominoes)")
_parser.add_argument("goal", type=int,
                     help="target area (0.001 is subtracted internally); must match the Step 1/2/3 runs")
_parser.add_argument("--only", type=int, metavar="SET_NUMBER", default=None,
                     help="process only this vector set (1-based number in the Step 2 file)")
_parser.add_argument("--max-solutions", type=int, metavar="N", default=10000,
                     help="stop after generating N different solutions (default: 10000)")
_parser.add_argument("--storage-directory", type=str, metavar="DIR", default=".",
                     help='directory for temporary files for intermediate store of the DP table (default: ".")')
_parser.add_argument("--compress", action='store_false',
                     help='use compression for temporary files (not implemented) (default: off)')
args = _parser.parse_args()
n_max, goal = args.n_max, args.goal
target_area = goal - 0.001

CANDIDATES = dp.load_step_result(3, n_max, goal).RESULT
if args.only is None:
    SELECTED = [c for c in CANDIDATES if c["reaches_goal"]]
    if not SELECTED:
        sys.exit(f"No vector set of the Step 3 file reaches the goal {goal}. Nothing to reconstruct. "
                 f"(Use --only to reconstruct a feasible vector set that does not reach the goal.)")
else:
    SELECTED = [c for c in CANDIDATES if c["set_number"]==args.only]
    if not SELECTED:
        _parser.error(f"--only: vector set #{args.only} is not feasible. Feasible vector sets "
                      f"according to Step 3: {[c['set_number'] for c in CANDIDATES]}")

print(f"Canonical-order DAG, {target_area=}")
print(f"Reconstruct {len(SELECTED)} of {len(CANDIDATES)} feasible vector sets: "
      f"{[c['set_number'] for c in SELECTED]}")




for c in SELECTED:
    r = canon.solve_vector_set(c["set_number"], c["vectors"], c["neg_penalty"], target_area, keep_tables=True)
    if r.enclosed_area is None or abs(r.enclosed_area-c["enclosed_area"])>1e-6:
        sys.exit(f"Inconsistent with the Step 3 file: vector set #{c['set_number']} has enclosed area "
                 f"{c['enclosed_area']} there, but {r.enclosed_area} was computed now. "
                 f"Rerun Step 3 with the same arguments.")
    print(f"Number of optimal solutions according to the forward pass: {r.count}")
    SOLS = canon.reconstruct_dag_canonical(r.Vecs, r.TABLES, r.finals)
    num_dag = sum(e[0] for e in SOLS[len(r.Vecs)].values())
    if num_dag != r.count:
        sys.exit(f"Internal inconsistency: DAG has {num_dag} paths, forward pass counted {r.count}.")
    canon.enumerate_canonical_solutions(SOLS, len(r.Vecs), r.shift, args.max_solutions)

print(f"Reconstructed {len(SELECTED)} vector set(s).")
