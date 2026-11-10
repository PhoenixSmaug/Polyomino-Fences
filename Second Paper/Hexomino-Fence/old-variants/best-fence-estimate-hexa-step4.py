"""
STEP 4. Reconstruction of the optimal reference polygons.
Take the feasible vector sets found by Step 3
(File best_fence_candidate_solutions_step_3_n<n>t<goal>.py),
run the dynamic program on them again, this time keeping all tables,
and reconstruct one or all optimal reference polygons.
=================================
Usage: pypy3 best-fence-estimate-hexa-step4.py <n_max> <goal> [--only SET_NUMBER] [--all]

  n_max       size of the polyominoes (6 for hexominoes)
  goal        target area; must be the same as in the Step 1, 2 and 3 runs
              that produced the input file.
  --only SET_NUMBER
              process only the given vector set, instead of the default selection.
              The number is 1-based, in the order of the Step 2 file, and is the same as
              the "Vector set #k" headings printed by Step 3 (and listed in the README).
              The set must be feasible, i.e. listed in the Step 3 file. It need not reach the goal.
              Default: all vector sets that reach the goal.
  --all       reconstruct ALL optimal reference polygons (also accepted as "-all").
              Default: reconstruct just one.

Reconstruction keeps all tables in memory. With --all, in addition, all partial solutions are
stored, which needs much more memory. The program best-fence-estimate-hexa-step4-count.py
computes the same set of solutions with less memory, and it counts them.
"""
import sys
import argparse
from math import atan2

import best_fence_estimate_hexa_dp as dp

_parser = argparse.ArgumentParser(
    description="Step 4: reconstruct the optimal reference polygons of the feasible vector sets of Step 3.")
_parser.add_argument("n_max", type=int, help="size of the polyominoes (6 for hexominoes)")
_parser.add_argument("goal", type=int,
                     help="target area (0.001 is subtracted internally); must match the Step 1/2/3 runs")
_parser.add_argument("--only", type=int, metavar="SET_NUMBER", default=None,
                     help="process only this vector set (1-based number in the Step 2 file)")
_parser.add_argument("--all", "-all", dest="all_solutions", action="store_true",
                     help="reconstruct all optimal reference polygons (default: just one)")
args = _parser.parse_args()
n_max, goal = args.n_max, args.goal
target_area = goal - 0.001
FIND_ALL_SOLUTIONS = args.all_solutions

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

print(f"Constructing all solutions: {FIND_ALL_SOLUTIONS}, {target_area=}")
print(f"Reconstruct {len(SELECTED)} of {len(CANDIDATES)} feasible vector sets: "
      f"{[c['set_number'] for c in SELECTED]}")

for c in SELECTED:
    r = dp.solve_vector_set(c["set_number"], c["vectors"], c["neg_penalty"], target_area, keep_tables=True)
    if r.enclosed_area is None or abs(r.enclosed_area-c["enclosed_area"])>1e-6:
        sys.exit(f"Inconsistent with the Step 3 file: vector set #{c['set_number']} has enclosed area "
                 f"{c['enclosed_area']} there, but {r.enclosed_area} was computed now. "
                 f"Rerun Step 3 with the same arguments.")
    Vecs, Corners, shift = r.Vecs, r.Corners, r.shift

    if FIND_ALL_SOLUTIONS:
        SOLS = dp.reconstruct_all_lists(Vecs, r.TABLES, r.finals)
    else:
        bestsol, bestvol = r.finals[-1]
        SOLS = [dp.reconstruct_one(Vecs, r.TABLES, bestsol, bestvol)]
    print(len(SOLS),"solutions")
    ALL_SOLS = set()
    for SOL in SOLS:
        dp.restore_corners(SOL, Corners, shift)
        if not FIND_ALL_SOLUTIONS:
            print(f"{SOL=}")
            print(f"{sum(x for x,y in SOL)=}; {sum(y for x,y in SOL)=}", flush=True)
        assert sum(x for x,y in SOL)==sum(y for x,y in SOL)==0

        ALL_SOLS.add(dp.canonical(SOL))
    if FIND_ALL_SOLUTIONS:
        for i,S in enumerate(ALL_SOLS):
            SS = [(x,y) for _,x,y in sorted((atan2(y,x),x,y) for x,y in S)]
            print(i,SS)

print(f"Reconstructed {len(SELECTED)} vector set(s).")
