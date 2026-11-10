"""
Faster version: only the bound ZF (Zauberformel) is used.
Recursive generation, taking into account all choices.
No preliminary pruning of reference-vector lists.

A) Absolute rotations are considered when generating reference vectors. (Rotations not relevant for this program)
B) Penalty is distributed equally in ZF-computation
C) Reference vectors remain attached to midpoint and are not moved to inner endpoint of attaching segment.
   (Otherwise the vectors with rot=-1 get an unfair advantage in length.)
=================================
Usage: pypy3 best-fence-step1.py <n> <goal> [--cache] [--progress] [--check]

  n           size of the polyominoes (6 for hexominoes)
  goal        target area; the program subtracts 0.001 to guard against rounding errors
  --cache     load/save the precomputed vector lists and ZF_i_and_after in a pickle file
              (<tempdir>/step-1-<n>.pkl); off by default
  --progress  show progress counters (off by default)
  --check     slow self-test: compare the incremental bound with the classical
              computation ZF(...) at every node (off by default)
"""

import math
import sys
from watch_memory import memory_checkpoint
import argparse
_parser = argparse.ArgumentParser(description="Step 1: enumerate all vector sets whose upper bound reaches the goal.")
_parser.add_argument("n", type=int, help="size of the polyominoes (6 for hexominoes)")
_parser.add_argument("goal", type=int, help="target area (0.001 is subtracted internally)")
_parser.add_argument("--cache", action="store_true", help="load/save precomputed vector lists in a pickle file")
_parser.add_argument("--progress", action="store_true", help="show progress counters")
_parser.add_argument("--check", action="store_true", help="slow self-test of the incremental bound")
args = _parser.parse_args()
n_max, goal = args.n, args.goal
target = goal - 0.001

# step 0 reads its parameters from sys.argv when it is imported
sys.argv = [sys.argv[0], str(n_max)] + (["--progress"] if args.progress else [])

from best_fence_step0 import vectors, Letters6, Pieces, vmax, Kmax

#print("target area", target)

max_rot = max(rot for rr in vectors for (x,y),rot,b in rr)
min_rot = min(rot for rr in vectors for (x,y),rot,b in rr)
#print(f"Rotation ranges from {min_rot} to {max_rot}.")

## print some vector sets, for illustration
i=0
for p,v,rr in zip(Pieces,vmax,vectors):
    if i==3:
        break
    if i==0:
        #print("\nSome pieces, and their vectors (for illustration)")
        pass
    #print ("-"*30,f"\nVectors for piece {i} ({Letters6[i]}): {p}, max. length {v}")
    for r in sorted(set(rot for (x,y),rot,b in rr)):
        for ((x,y),rot,b),pen in rr.items():
            if rot==r:
                #print (f" {math.sqrt(x**2+y**2):5.3f} {(x,y)}, {rot=}, penalty {pen}")
                pass
    i += 1


def combine(u,v):
    "mixed product L1-norm x Lmax-norm"
    ux,uy=map(abs,u)
    vx,vy=map(abs,v)
    return max(max(ux,uy)*(vx+vy), (ux+uy)*max(vx,vy))/16

def combine2(up,vp):
    "mixed product L1-norm x Lmax-norm, combined with penalties"
    u,pen_x = up
    v,pen_y = vp
    ux,uy=map(abs,u)
    vx,vy=map(abs,v)
    return max(max(ux,uy)*(vx+vy), (ux+uy)*max(vx,vy))/8 +(pen_x+pen_y)/Nvecs
def combine2_self(up):
    "mixed product L1-norm x Lmax-norm, combined with penalties"
    u,pen_x = up
    ux,uy=map(abs,u)
    return max(ux,uy)*(ux+uy)/16 +pen_x/Nvecs

def ZF(vectors,allowed_rotations=None): # Zauberformel from Minkowski-mean
    sum_best = 0
    N = len(vectors)
    for i,ref_vectors in  enumerate(vectors):
        for j,ref_vectors2 in  enumerate(vectors):
            if i==j:
                best = max(combine(p,p) - pen/N for (p,rot),pen in ref_vectors.items()
                           if not allowed_rotations or rot in allowed_rotations)
            else:
                best = max(combine(p,p2)-(pen+pen2)/(2*N) for (p,rot),pen in ref_vectors.items()
                           for (p2,rot2),pen2 in ref_vectors2.items()
                           if not allowed_rotations or
                           (rot in allowed_rotations and rot2 in allowed_rotations))
            # penalties are distributed over all NxN terms
            sum_best += best
    return sum_best

print("\nStart with unpruned vector lists")

def make_no_rot(refvecs):
    best_pen = dict()
    for (v,rot,b),pen in refvecs.items():
        if v not in best_pen or pen<best_pen[v]:
            best_pen[v] = pen # minimum penalty for each vector, regardless of rotation and parity
    return sorted(best_pen.items())

import pickle, tempfile, os

_cache_path = os.path.join(tempfile.gettempdir(), f"step-1-{n_max}.pkl")
if args.cache and os.path.exists(_cache_path):
    print(f"Loading cached vectors_no_rot and ZF_i_and_after from {_cache_path}")
    with open(_cache_path, "rb") as _f:
        vectors_no_rot, Nvecs, ZF_i_and_after = pickle.load(_f)
else:
    vectors_no_rot = [make_no_rot(r) for r in vectors]
    print("vector lists without regard to rotation:")
    print([len(r) for r in vectors_no_rot])
    print("total length of lists", sum(len(r) for r in vectors_no_rot))

    ## FASTER VERSION only ZF
    Nvecs = len(vectors) # number of pieces (35 for hexominoes)
    ZF_i_and_after = [0]*(Nvecs+1)
    print(f"Computing ZF_i_and_after for {Nvecs} pieces...", flush=True)
    _report_every = max(1, Nvecs // 200)
    for i in range(Nvecs-1,-1,-1):
        if args.progress and (Nvecs-1-i) % _report_every == 0:
            pct = 100*(Nvecs-i)//Nvecs
            print(f"\r  ZF_i_and_after: {Nvecs-i}/{Nvecs} ({pct}%)", end="", flush=True)
        ref_vectors = vectors_no_rot[i]
        sum_best = max(combine(p,p) - pen/Nvecs for p,pen in ref_vectors)
        for j in range(i+1,Nvecs):
            ref_vectors2 = vectors_no_rot[j]
            best = max(2*combine(p,p2)-(pen+pen2)/Nvecs
                       for p,pen in ref_vectors
                       for p2,pen2 in ref_vectors2)
            # penalties are distributed over all NxN terms
            sum_best += best
        ZF_i_and_after[i] = ZF_i_and_after[i+1] + sum_best
    if args.progress:
        print(f"\r  ZF_i_and_after: done.                    ", flush=True)

    if args.cache:
        print(f"Saving cache to {_cache_path}")
        with open(_cache_path, "wb") as _f:
            pickle.dump((vectors_no_rot, Nvecs, ZF_i_and_after), _f)


    
#print(ZF_i_and_after)

solution = [None]*Nvecs
solution_indices = [-1]*Nvecs

ZF_after = [ [0]*(Nvecs+1) for i in range(Nvecs) ]
ZF_after_accumulated = [ [0]*(Nvecs) for i in range(Nvecs+1) ]

def recursive_set_r_norot(level, ZF_so_far):
    global nr_survive,nr_calls
    nr_calls +=1
    if args.progress and (level<=3 or nr_calls%100000==0):
        print("P",*solution_indices[:level],"/",nr_calls, flush=True)
    if level==Nvecs:
        nr_survive += 1
        print(nr_survive,"ZF=",f"{ZF_so_far:.3f}", solution, flush=True)
        print("A(",(nr_survive,ZF_so_far, solution),')', file = sol_file, flush=True)
        return
    for i,(v,pen) in enumerate(vectors_no_rot[level]):
        solution[level] = (v,pen) # currently selected vector for this level
        solution_indices[level] = i
        # combination of current with all previous (and itself)
        ZF_contribution_before = ( combine(v,v)-pen/Nvecs
                            + sum(2*combine(v,v2)-(pen+pen2)/Nvecs
                                    for v2,pen2 in solution[:level]) )
        # combination of current with all following:
        for j in range(level+1,Nvecs):
            contrib_i_j = max(2*combine(v,v2)-(pen+pen2)/Nvecs
                                for v2,pen2 in vectors_no_rot[j])
            ZF_after[level][j] = contrib_i_j
            ZF_after_accumulated[level][j] = ZF_after_accumulated[level-1][j] + contrib_i_j
            # ( for level = 0: ZF_after_accumulated[level-1] refers to the last row, which is 0
        ZF_contribution_later = sum(ZF_after_accumulated[level][level:])
        ZF_total = ZF_so_far + ZF_contribution_before + ZF_contribution_later + ZF_i_and_after[level+1]
        if args.check: # TEST: compare with classical computation of ZF: (seems to work)
            vectors = [{(p,0): pen} for p,pen in solution[:level+1]] + [
                {(p,0): pen for p,pen in rr} for rr in vectors_no_rot[level+1:] ]
            assert(abs(ZF_total-ZF(vectors))<1e-9)
        if ZF_total >= target:
            recursive_set_r_norot(level+1, ZF_so_far + ZF_contribution_before)
            
sol_file_name = f"best_fence_candidate_solutions_step_1_n{n_max}t{goal}.py"
with open(sol_file_name,'w') as sol_file:
    print(f"""# generated by the program {sys.argv[0]} with target area {target}
B=[]
def A(x):
  B.append(x)
""", file = sol_file)

    nr_calls = nr_survive = 0
    recursive_set_r_norot(0, 0)
    print(nr_survive, "solutions found")
    print("# generation completed.", file = sol_file)
print(f'Solutions were written to file "{sol_file_name}".')
print("target area", target)
print("runtime",memory_checkpoint())
