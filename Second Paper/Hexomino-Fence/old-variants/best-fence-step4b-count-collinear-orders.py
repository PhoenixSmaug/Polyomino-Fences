import ast, sys, math, collections
from fractions import Fraction
from math import factorial, gcd
polys = []
for line in open(sys.argv[1]):
    s = line.strip()
    if not s or not s[1].isdigit(): continue
    #i, poly = ast.literal_eval(s)
    polys.append(ast.literal_eval(s))
print(len(polys), "polygons")

target = collections.Counter({(0,5):2,(0,6):1,(1,5):3,(0,4):3,(1,3):4,(1,4):11,(2,3):1,(2,4):6,(2.5,2.5):4})

def fold(v): return tuple(sorted(map(abs, v))) # fold to second octant
def direction(v):
    "direction normalized to max-norm 1"
    x, y = Fraction(v[0]), Fraction(v[1])
    L = max(abs(x), abs(y))
    return (x/L, y/L)
piece_factor = factorial(2)*factorial(3)*factorial(3)*factorial(11)*factorial(4)*factorial(6)*factorial(5)
bad = 0; total_orders_factor = 0
canon = set()
# for i, P in polys: ## old format, without correction
for i, correction, P in polys:
    assert collections.Counter(fold(v) for v in P) == target, i
    sx = sum(Fraction(v[0]) for v in P); sy = sum(Fraction(v[1]) for v in P)
    assert (sx, sy) == (0, 0), i
    # area via shoelace on the given edge order
    x = y = Fraction(0); A2 = Fraction(0)
    for vx, vy in P:
        nx, ny = x+Fraction(vx), y+Fraction(vy); A2 += x*ny - nx*y; x, y = nx, ny
    assert A2 == 2*1663, (i, A2/2)
    # convex / sorted by angle is not checked. #  angles = [math.atan2(v[1], v[0]) for v in P]
    # collinear groups
    groups = collections.defaultdict(list)
    for v in P: groups[direction(v)].append(v)
    mult = 1
    for d, vs in groups.items():
        c = collections.Counter(vs)
        if len(c) > 1:
            n = len(vs); m = factorial(n)
            for k in c.values(): m //= factorial(k)
            mult *= m
    if mult > 1: bad += 1
    total_orders_factor += mult
    canon.add(tuple(sorted(P)))
print("vector multiset:",target)
print("all closed, area 1663, multiset = set #10: OK;", len(canon), "distinct multisets of oriented vectors")
print("polygons with two different collinear vectors on at least one side:", bad)
print("sum over polygons of collinear-reordering factors:", total_orders_factor, "(paper originally used 244)")
print("orders, paper (old):", piece_factor*244)
print("orders, corrected:  ", piece_factor*total_orders_factor)
fences = piece_factor*total_orders_factor*2**10
print("fences, corrected:  ", fences, f"= {fences:.4e}")
# Lemma 5 condition 3: each axis direction present
miss = [i for i,_,P in polys if not all(any(direction(v)==d for v in P) for d in [(1,0),(0,1),(-1,0),(0,-1)])]
print("polygons missing an axis direction:", miss)
if not miss:
    print("Lemma 5, condition 3 is satisfied: All axis directions are present in each reference polygon.")
else:
    print("**** Lemma 5, condition 3 is not satisfied! ****")
dist = collections.Counter(); sides = collections.Counter()
for i, _, P in polys:
    groups = collections.defaultdict(list)
    for v in P: groups[direction(v)].append(v)
    mult = 1; description = []
    for d, vs in groups.items():
        # vs is the list of parallel vectors
        c = collections.Counter(vs)
        n = len(vs); m = factorial(n) # multinomial coefficient
        for k in c.values(): m //= factorial(k)
        mult *= m
        if d[0]*d[1]!=0: # not axis-parallel
            assert len(c)==1 # no distinct parallel vectors among those
        else:
            description.append(tuple(sorted(fold(v)[1] for v in vs)))
    dist[mult,tuple(sorted(description))] += 1
    for dd in description: sides[dd] += 1
    if any(d not in [(1,0),(0,1),(-1,0),(0,-1)] and len(set(vs))>1 for d,vs in groups.items()): print("non-axis mixed in", i)
print("factor distribution:")
for (f,x),k in sorted(dist.items()):
    print(f"{k:3} x factor {f}, edges", *[list(e) for e in x if len(e)==1], *[list(e) for e in x if len(e)!=1])
print("axis-parallel sides overall:", dict(sides))
