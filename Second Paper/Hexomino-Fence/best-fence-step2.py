"""STEP 2. Take the list of vectors computed without caring about rotations and parity (Step 1),
(File best_fence_candidate_solutions_step_1_n<n>t<goal>.py)
and take into account the rotations and the parity.
find all potential vector sequences
===================================
Usage: pypy3 best-fence-step2.py <n> <goal> [--progress] [--verbose]

  n           size of the polyominoes (6 for hexominoes)
  goal        target area;
              must be the same as in the Step 1 run that produced the input file.
  --verbose   list, for every piece, all reference vectors that fit the vector classes
              used in the surviving candidates (off by default). The consistency checks
              behind this listing are always performed; only the printout is optional.
  --progress  show progress counters while Step 0 generates the reference vectors
              (off by default)

Reads  best_fence_candidate_solutions_step_1_n<n>t<goal>.py   (from Step 1)
Writes best_fence_candidate_solutions_step_2_n<n>t<goal>.py   (for Step 3)
"""

from collections import defaultdict, Counter
from itertools import combinations
import importlib
import math
import sys

import argparse
_parser = argparse.ArgumentParser(
    description="Step 2: take rotation and parity into account and reduce the Step 1 candidates "
                "to the remaining vector sets.")
_parser.add_argument("n", type=int, help="size of the polyominoes (6 for hexominoes)")
_parser.add_argument("goal", type=int,
                     help="target area; must match the Step 1 run")
_parser.add_argument("--progress", action="store_true",
                     help="show progress counters while Step 0 generates the reference vectors")
_parser.add_argument("--verbose", action="store_true",
                     help="list all fitting reference vectors for every piece")
args = _parser.parse_args()
n_max, goal = args.n, args.goal
target = goal

# step 0 reads its parameters from sys.argv when it is imported
sys.argv = [sys.argv[0], str(n_max)] + (["--progress"] if args.progress else [])

from best_fence_step0 import vectors, vectors_precise, print_piece
importname = f"best_fence_candidate_solutions_step_1_n{n_max}t{goal}"
print("reading candidates from file", importname+".py")
try:
    _cand_mod = importlib.import_module(importname)
except ModuleNotFoundError as e:
    if e.name != importname:
        raise
    sys.exit(f'Input file "{importname}.py" not found. Run Step 1 first with the same '
             f'arguments (n={n_max}, goal={goal}) and run this program from the directory '
             f'that contains the Step 1 result file and the scripts.')
candidates = _cand_mod.B
print(len(candidates), "candidates read.")

def combine(u,v):
    "mixed product L1-norm x Lmax-norm. (The vectors u and v can come in any order.)"
    ux,uy=map(abs,u)
    vx,vy=map(abs,v)
    return max(max(ux,uy)*(vx+vy), (ux+uy)*max(vx,vy))/16

def F_M(vectors): # "Zauberformel" F_M, formula (8) in the paper
    return sum(combine(p,q) for p in vectors for q in vectors)

print("\nStart with unpruned vectors")

def update_max(Dict,key,value):
    if key not in Dict or value>Dict[key]:
        Dict[key] = value

def best_rot4m_correction(vectors):
    # dynamic programming to find the largest sum of corrections
    # subject to the condition that the total rotation is 360° (=4),
    # and the total parity is even.
    sum_best = {(0,0):0}
    for ref_vectors in vectors:
        #print("s",vectors.index(ref_vectors),len(ref_vectors),end="")
        sum_best_new = dict()
        for ((x,y),rot,par),pen in ref_vectors:
            corr = -pen
            for (nr,npar),v in sum_best.items():
                update_max(sum_best_new, (nr+rot, (npar+par)%2), v+corr)
        sum_best = sum_best_new
    return sum_best.get((4,0),-10000)

def vector_fingerprint(sol): # unique identifier for a vector multiset
    return tuple(sorted(Counter(v for v,p in sol).items()))

vectors_per_piece = [set() for _ in vectors] # statistics: how each piece can be used.

num_rem = 0
Vector_sets = defaultdict(int)    

for i,bound,solution in candidates:
    # compute the feasible set U_j for each piece j.
    vecs_with_rot = [ [((v2,rot,par),p2)  for (v2,rot,par),p2 in rr.items() if v2==v]
                       for (v,p),rr in zip(solution, vectors)]
    
    rot_corr = best_rot4m_correction(vecs_with_rot)
    # largest correction considering rotation and parity
    # = Delta^rot in the paper
    if rot_corr==-10000:
        continue # no feasible solution at all
    correction = -sum(pen for v,pen in solution) # correction is ADDED to the area
    vecs =  [x for x,p in solution]
    Z = F_M(vecs)
    
    assert abs(Z+correction - bound)<1e-9 # rounding errors
    bound = Z+correction # only multiplication and addition with half-integers.
                         # This should be exact with float.
    diff = rot_corr-correction
    if bound+diff<target:
        continue                                 

    
    for (vpp,(v,p)) in zip(vectors_per_piece,solution):
        vpp.add(v) # maintain overview statistics

    num_rem += 1
    print(f"#{num_rem} input {i}, {bound=}", #m4bound,
          f"{correction=}, ∆_rot={rot_corr}, {diff=}, {bound+diff=}")
    
    Vector_sets[vector_fingerprint(solution),bound,rot_corr] += 1

# Analyze the used vectors per piece
print("\nChecking the used vectors per piece:")
print("Check that only rotations +1 and -1 and 0 occur, and no weird situations.")
"""
The check fails for parameters n=5 and t=128 (optimal pentominoes)
because some weird situations that the program cannot handle cannot be
excluded.
The following example shows the W-pentomino with rotation +90° (+1),
where the reference vector is (2.5,-1.5) when the incoming attachment
edge points vertically downward:

   XXX
 ->XXX
   XXX
   XXXXXX
   XXXXXX ^ 
   XXXXXX |
      XXXXXX
      XXXXXX
      XXXXXX

The program is not prepared for this unusual situation, where the
reference vector has a negative component.  After reflecting the
pentomino at a -45° line, with the same attachment edges in opposite
roles, one gets the reference vector (-1.5,2.5), with rotation +90°.

The P pentomino allows a similar situation with reference vector
(2.5,-0.5), and rotation +90°.

Another bad case that the program cannot handle is a reference vector
with $|y|>x$ with rotation 0°. The X-pentomino allows a reference vector
(1,2), and the L and N pentominoes can have (2,3).
(For hexominoes, such "weird" possibilities of using a piece exist in
principle, but vector combinations that use a piece in that way are
excluded as not good enough by the bounding formulas.)
"""

VERBOSE = args.verbose
for i,(vecs,vp) in enumerate(zip(vectors_per_piece,vectors_precise), start=1):
    if VERBOSE: print("piece",i,"*********")
    for v0 in vecs:
        num_fit = 0
        if VERBOSE: print(" piece",i,f"used vector: {v0}. All vectors for piece",i,"that fit:")
        (x0,y0) = v0
        is_corner = x0!=int(x0)
        for ((x,y),rot,_),pen in vp.items():
            if tuple(sorted(map(abs,(x,y))))==v0:
                num_fit += 1
                if VERBOSE: print("  ",(x,y),rot,pen)
                if is_corner:
                    assert rot in (-1,1) # cannot deal with strange situations
                else:
                    assert rot in (0,)
                if rot==0:
                    assert not is_corner
                    assert x>0
                    assert x>=abs(y)
                if rot==+1:
                    assert is_corner
                    assert x>0
                    assert y>0
                if rot==-1:
                    assert is_corner
                    assert x>0
                    assert y<0
        if is_corner:
            if x0==y0:
                assert num_fit == 2
            else:
                assert num_fit == 4
        else:
            if x0==0:
                assert num_fit == 1
            else:
                assert num_fit == 2
print("Checks passed.\n")
                    
        
A = set(list(Vector_sets)[0][0])
for i,j,_ in Vector_sets:
    A.intersection_update(i)
A = tuple(sorted(A))
print("Common parts to all potential solutions:  ",A)
Vector_sets_list = sorted((bound,sol,Delta_rot,m) for (sol,bound,Delta_rot),m in Vector_sets.items())
print(len(Vector_sets), "vector sets.")


Result = []
Barely = []

for running,(bound,sol,Delta_rot,m) in enumerate(Vector_sets_list, start=1):
    perimeter = sum(math.sqrt(x**2+y**2)*mu for (x,y),mu in sol)  # perimeter only for information
    print(f"\nVectors #{running}:", A + tuple(s for s in sol if s not in A),f"{perimeter=:.4f} ({m} times)")

    solution = [v for v,m in sol for _ in range(m)] # expand multiplicities
    cornercandidates = [(x,y) for (x,y) in solution if x!=int(x)]
    rest = [(x,y) for x,y in solution if x==int(x)]
    ncorners = len(cornercandidates)
    assert ncorners >= 4 and ncorners%2 == 0
    if ncorners > 4:
        print("Number of corner candidates is larger than 4. Try all sets of reflex/convex combinations")
    Cornersets = set(frozenset(x for x in combinations(cornercandidates, (4+ncorners)//2)))
    # |corners|-|reflex| == 4
    Surviving_vector_set = False # refers to the vector set as a whole
    Barely_this_set = []
    for corners in Cornersets:
    # in the standard case with only 4 corner candidates, there is only one choice
        reflex = list(cornercandidates)
        for el in corners:
            reflex.remove(el)
        if ncorners > 4:
            print ("Try positive corner pieces:", *corners,"/ negative (reflex) pieces:",*reflex)
        Exclude_corner_set = barely_Exclude_corner_set = False
        for shift in 1.5,1,0.5,0,-0.5,-1,-1.5:
            # "inside" seems to lead to stronger bounds??
            corners_new = [(x+shift,y+shift) for x,y in corners] #move the reference point
            reflex_new = [(x-shift,y-shift) for x,y in reflex]

            if any (x<=0 or y<=0 for (x,y) in corners_new+reflex_new):
            #if any (x<0 or y<0 for (x,y) in corners_new+reflex_new):
                continue # (but we might also allow 0 coordinates)
            shifttext = f"shifting reference points by {shift}"
            if shift>0:
                shifttext += " (outwards)"
            if shift<0:
                shifttext += " (i.e., inwards)"

            correction_shifted = Delta_rot # to be added to the area of the reference polygon
            for (x,y) in corners:
                correction_shifted -= (x+shift)*(y+shift)/2 - x*y/2
            for (x,y) in reflex:
                correction_shifted += (x-shift)*(y-shift)/2 - x*y/2
                # taken with opposite sign because enclosed area is on the other side.
            for (x,y) in rest:
                assert x<=y
                correction_shifted -= y*shift # y is the larger coordinate

            vecsnew = rest + corners_new + reflex_new
            Z = F_M(vecsnew)
            bound_shifted = Z + correction_shifted
            if bound_shifted < target:
                sign, comparison = "*<"
                Exclude_corner_set = True
            elif bound_shifted == target:
                sign, comparison = "=="
                barely_Exclude_corner_set = True
                Barely_this_set.append ((running, list(corners), reflex, rest))
            else:
                sign, comparison = " >"
            print("%-50s"% (f"{sign}  {Z:.5f}{correction_shifted:+7.2f}="
                f"{Z+correction_shifted:.5f} {comparison} {target=}"),
                  shifttext)

        if Exclude_corner_set:
            print("Excluded.")
        elif barely_Exclude_corner_set:
            print("Barely excluded.")
            Barely += Barely_this_set
        else:
            Surviving_vector_set = True
            print("Not excluded.")
    if Surviving_vector_set:
        Result.append((A + tuple(s for s in sol if s not in A), Delta_rot))
        print(f"Surviving vector set. (#{len(Result)})")

    if ncorners!=4 and Surviving_vector_set:
        print("The number of corner candidates is not 4, but the vectors cannot be excluded.")
        print("The program for Step 3 is not prepared for this. Exit.")
        exit()

sol_file_name = f"best_fence_candidate_solutions_step_2_n{n_max}t{goal}.py"

with open(sol_file_name,'w') as sol_file:
    print(f"""# generated by the program {sys.argv[0]} with {target=}
RESULT=[""", file = sol_file)
    for s in Result:
        print(s,",", file = sol_file)
    print("]", file = sol_file)
print(f'{len(Result)} vector sets were written to file "{sol_file_name}".')
if Barely:
    print("\n*** The following vector sets were excluded under the assumption that")
    print("they cannot form a symmetric polygon with eight axes-preserving symmetries.\n")

    for runningcount, corners, reflex, rest in Barely:
        print(f"set #{runningcount}:",
              *tuple(sorted(Counter(corners+reflex+rest).items())))
        print("reflex:",reflex)
