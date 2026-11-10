# General-offset version of exact_reflex_bound.py.
# Reference points are moved by s along every attachment edge (s<0 inwards, s>0 outwards,
# s=-1/2 inner endpoint, s=+1/2 outer endpoint). For precise vector (x,y) in the folded class:
#   straight (rot 0, parallel edges, x>=|y|):  v unchanged,          correction += -s*x
#   corner   (rot +1):                         v -> (x+s, y+s),      correction += -(s*(x+y)+s*s)/2
#   reflex   (rot -1):                         v -> (x-s, y-s),      correction += (-s*(x+y)+s*s)/2
# Any common s gives a valid upper bound (Lemma 3 works for any common point on the line of
# each attachment edge), so the minimum over s is valid too.
import sys, importlib, collections
sys.path[:0] = [sys.argv[1], sys.argv[2]]
OFFSETS = [float(a) for a in sys.argv[3].split(',')] if len(sys.argv) > 3 else [k/20 for k in range(-30, 31)]
sys.argv = ['x', '6']
from best_fence_step0 import vectors
B = importlib.import_module('best_fence_candidate_solutions_step_1_n6t1597').B
def combine(u,v):
    ux,uy=map(abs,u); vx,vy=map(abs,v)
    return max(max(ux,uy)*(vx+vy),(ux+uy)*max(vx,vy))/16
def ZF(vs): return sum(combine(p,q) for p in vs for q in vs)
target = 1596.999
FA=tuple(sorted({(0.0,5.0):2,(0.0,6.0):1,(1.0,5.0):3,(0.0,4.0):3,(1.0,3.0):4,(1.0,4.0):11,(2.0,3.0):1,(2.0,4.0):4,(2.5,2.5):4,(2.5,3.5):2}.items()))
FB=tuple(sorted({(0.0,5.0):2,(0.0,6.0):1,(1.0,5.0):3,(0.0,4.0):3,(1.0,3.0):4,(1.0,4.0):11,(2.0,4.0):5,(2.5,2.5):5,(2.5,3.5):1}.items()))
#FA can be done, FB is hard
def bound(sol, half, refl, s):
    par_best = {0: 0.0}; vs = []
    for k,(v,p) in enumerate(sol):
        x,y = v
        if k == refl:   want, adj, w = -1, (-s*(x+y)+s*s)/2, (x-s, y-s)
        elif k in half: want, adj, w =  1, -(s*(x+y)+s*s)/2, (x+s, y+s)
        else:           want, adj, w =  0, -s*max(x,y), (x, y)
        vs.append(w)
        opts = {}
        for (v2,rot,par),pen in vectors[k].items():
            if v2==v and rot==want:
                opts[par] = max(opts.get(par,-1e9), -pen+adj)
        new = {}
        for q,val in par_best.items():
            for par,o in opts.items():
                new[(q+par)%2] = max(new.get((q+par)%2,-1e9), val+o)
        par_best = new
    return ZF(vs) + par_best[0] if 0 in par_best else -1e9
worst = -1e9; worst_info = None; n = 0
for nr, b, sol in B:
    half = [k for k,(v,p) in enumerate(sol) if v[0]!=int(v[0])]
    if len(half) != 6: continue
    fp = tuple(sorted(collections.Counter(v for v,p in sol).items()))
    if fp not in (FA, FB): continue
    n += 1
    for refl in half:
        vals = {s: bound(sol, half, refl, s) for s in OFFSETS}
        s_best = min(vals, key=vals.get)
        if vals[s_best] > worst: worst, worst_info = vals[s_best], (nr, sol[refl][0], s_best, vals.get(-0.5), vals.get(0.5))
print(n, "sequences; offsets tried:", OFFSETS[0], "...", OFFSETS[-1], f"({len(OFFSETS)} values)")
print("largest bound (max over sequences and reflex choices of min over offsets):", worst,
      ">= target" if worst >= target else "< target  => all excluded")
print("attained at sequence", worst_info[0], "reflex class", worst_info[1], "best offset", worst_info[2])
