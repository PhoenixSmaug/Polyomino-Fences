"""
STEP 4 (variant with counting). Reconstruction of ALL optimal reference polygons.
Same task as best-fence-estimate-hexa-step4.py --all, but less memory-intensive:
instead of carrying along all partial solutions, the program builds a DAG
(directed acyclic graph) of the optimal partial solutions. This allows to COUNT all solutions
without generating them, and to generate them afterwards (also: random sampling would be possible,
although not quite uniformly, since symmetries might introduce a bias).
=================================
Usage: pypy3 best-fence-estimate-hexa-step4-count.py <n_max> <goal> [--only SET_NUMBER] [--max-solutions N]

  n_max       size of the polyominoes (6 for hexominoes)
  goal        target area; must be the same as in the Step 1, 2 and 3 runs
              that produced the input file.
  --only SET_NUMBER
              process only the given vector set, instead of the default selection.
              The number is 1-based, in the order of the Step 2 file, and is the same as
              the "Vector set #k" headings printed by Step 3 (and listed in the README).
              The set must be feasible, i.e. listed in the Step 3 file. It need not reach the goal.
              Default: all vector sets that reach the goal.
  --max-solutions N
              stop generating solutions when N different solutions (up to symmetry)
              have been generated. The counting is not affected. (default: 10000)
"""
import sys
import argparse

import best_fence_estimate_hexa_dp as dp

_parser = argparse.ArgumentParser(
    description="Step 4 (counting variant): count and generate all optimal reference polygons "
                "of the feasible vector sets of Step 3, using a DAG.")
_parser.add_argument("n_max", type=int, help="size of the polyominoes (6 for hexominoes)")
_parser.add_argument("goal", type=int,
                     help="target area (0.001 is subtracted internally); must match the Step 1/2/3 runs")
_parser.add_argument("--only", type=int, metavar="SET_NUMBER", default=None,
                     help="process only this vector set (1-based number in the Step 2 file)")
_parser.add_argument("--max-solutions", type=int, metavar="N", default=10000,
                     help="stop after generating N different solutions (default: 10000)")
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

print(f"Counting all solutions via a DAG, {target_area=}")
print(f"Reconstruct {len(SELECTED)} of {len(CANDIDATES)} feasible vector sets: "
      f"{[c['set_number'] for c in SELECTED]}")

for c in SELECTED:
    r = dp.solve_vector_set(c["set_number"], c["vectors"], c["neg_penalty"], target_area, keep_tables=True)
    if r.enclosed_area is None or abs(r.enclosed_area-c["enclosed_area"])>1e-6:
        sys.exit(f"Inconsistent with the Step 3 file: vector set #{c['set_number']} has enclosed area "
                 f"{c['enclosed_area']} there, but {r.enclosed_area} was computed now. "
                 f"Rerun Step 3 with the same arguments.")
    SOLS = dp.reconstruct_dag(r.Vecs, r.TABLES, r.finals)
    dp.enumerate_dag_solutions(SOLS, len(r.Vecs), r.shift, args.max_solutions)

print(f"Reconstructed {len(SELECTED)} vector set(s).")
