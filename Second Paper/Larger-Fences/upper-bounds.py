import sys
import itertools
from collections import defaultdict
import math
import gc
from decimal import Decimal, getcontext, ROUND_UP

n_max = int(sys.argv[1])
print(f"n_max = {n_max}")

DIRECTIONS = (0,-1),(1,0),(0,1),(-1,0) # 0=South, 1=East, 2=North, 3=West
# counterclockwise directions

best_known_area = {  # from greedy heuristic solver
    5: 128,
    6: 1597,
    7: 19_718, 
    8: 279_105,
    9: 4_003_811,
    10: 60_832_848,
    11: 934_201_230,
    12: 14_625_233_629,
    13: 229_830_309_620,
    14: 3_635_789_660_265,
    15: 57_656_581_776_211
}

target = best_known_area[n_max] + 0.999


def rotateL(P):
    return [(-y,x) for x,y in P]
def rotateR(P):
    return [(y,-x) for x,y in P]
def rotate2(P):
    return [(-x,-y) for x,y in P]
def rotate0(P): # identity
    return [(x,y) for x,y in P]
def rotate(P,dir):
    "rotate counterclockwise by 'dir' multiples of 90°"
    return (rotate0,rotateL,rotate2,rotateR)[dir%4](P)

def reflect(P): # reflect at x=y line
    return [(y,x) for x,y in P]
def normalize(P): # lowest coordinates become 0
    xmin = min(x for x,y in P)
    ymin = min(y for x,y in P)
    return [(x-xmin,y-ymin) for x,y in P]

def orient(a,b,c): # orientation predicate for three points
    (x1,y1),(x2,y2),(x3,y3) = a,b,c
    return (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)

def det(a,b): # area of parallelogram spanned by two vectors
    (x1,y1),(x2,y2) = a,b
    return x1*y2-y1*x2


## COMPUTE POLYOMINOES (Memory Optimized)
# Keep only the active generation in memory during expansion
current_generation = [((0,0),)]
for n in range(2, n_max+1):
    New = set()
    for P in current_generation:
        for x,y in P:
            for q in ((x+1,y),(x,y+1),(x-1,y),(x,y-1)):
                if q not in P:
                    P2 = P+(q,)
                    # normalize under rotations and reflections.
                    Q = min(tuple(sorted(normalize(P4))) for
                            P3 in (P2,rotateL(P2),rotateR(P2),rotate2(P2))
                            for P4 in (P3,reflect(P3)))
                    New.add(Q)
    
    current_generation = sorted(tuple(sorted(reflect(Q))) for Q in New)
    print(f"{n=}. There are {len(current_generation)} {n}-ominoes.")

Pieces = current_generation
num_polyominoes = len(Pieces)
gc.collect()


def compute_ref_vectors(P):
    """Generates the reference vector dictionary for a single piece."""
    boundary_successor = defaultdict(list)
    for (x,y) in P:
        for (dx,dy) in DIRECTIONS:
            xx,yy = x+dx,y+dy
            if (xx,yy) not in P:
                perpx,perpy = -dy/2,dx/2
                start,end=(((x+xx)/2+perpx, (y+yy)/2+perpy),
                           ((x+xx)/2-perpx, (y+yy)/2-perpy))
                boundary_successor[start].append(end)
                
    start = current = min(boundary_successor)
    path = [start]
    assert len(boundary_successor[current]) == 1
    
    current = boundary_successor[current][0]
    while current != start:
        path.append(current)
        if len(boundary_successor[current]) == 1:
            current = boundary_successor[current][0]
        else:
            dx,dy = current[0]-path[-2][0], current[1]-path[-2][1]
            dx,dy = dy,-dx # right by 90°
            next = current[0]+dx, current[1]+dy
            assert next in boundary_successor[current]
            current = next
            
    L = len(path)
    path0 = path
    ref_vectors = dict()
    
    for path in (path0, reflect(reversed(path0))):
        path2 = path + path
        area_under = []
        directions = []
        for (x1,y1),(x2,y2) in zip(path,path2[1:]):
            area_under.append((x2-x1)*(y1+y2))
            directions.append(DIRECTIONS.index(((x2-x1),(y2-y1))))
        area_under += area_under
        sum_area_under = [0]
        s = 0
        for a in area_under:
            s += a
            sum_area_under.append(s)
        assert sum_area_under[L] >= 2*n_max
        
        directions *= 2
        absolute_direction = [0]
        dir = 0
        for dir1,dir2 in zip(directions,directions[1:]):
            dir += {3:-1, 0:0, 1:1}[(dir2-dir1) % 4]
            absolute_direction.append(dir)
        assert absolute_direction[L] == -4
        assert len(absolute_direction) == 2*L
        
        for i in range(L):
            for j in range(i+1, i+L):           
                (x1,y1),(x2,y2) = path2[i:i+2]
                xS,yS = (x1+x2)/2, (y1+y2)/2
                dirS = DIRECTIONS.index(((x1-x2),(y1-y2)))
                (x1,y1),(x2,y2) = path2[j:j+2]
                xE,yE = (x1+x2)/2, (y1+y2)/2
                dirE = DIRECTIONS.index(((x2-x1),(y2-y1)))
                along_boundary = (area_under[i]+area_under[j])/2 + sum_area_under[j]- sum_area_under[i+1]
                a_under_ref_vector = (xE-xS)*(yE+yS)
                penalty = (along_boundary-a_under_ref_vector)/2
                rot_old = (dirE-dirS) % 4
                rot = absolute_direction[j]-absolute_direction[i] + 2
                assert (rot-rot_old)%4 == 0
                [v] = rotate([(xE-xS, yE-yS)], -dirS)
                v = tuple(sorted(map(abs,v)))
                
                key = (v,rot)
                if key not in ref_vectors or penalty < ref_vectors[key]:
                    ref_vectors[key] = penalty
                    
    return ref_vectors


# PASS 1: Calculate global vmax and Kmax without storing dictionaries
print("\nPASS 1: Calculating boundary values...")
vmax = []
Kmax = []
total_unpruned_vectors = 0

for P in Pieces:
    ref_vectors = compute_ref_vectors(P)
    total_unpruned_vectors += len(ref_vectors)
    vmax.append(max(math.sqrt(x**2+y**2) for (x,y),rot in ref_vectors))
    Kmax.append(-min(ref_vectors.values()))

# Rigorous upper bound for u_max with exact integer arithmetic.
# Every longest vector has integer coordinates (Lemma 4), so |v|^2 = m is an integer.
# sqrt(m) is rounded UP to DIGITS decimal places, and the sum is exact.
DIGITS = 20
SCALE = 10**DIGITS
assert num_polyominoes/SCALE < 1e-6
# addition of a number for each peace should give a total error <1e-6

def sqrt_up_scaled(m):
    "smallest integer >= sqrt(m)*SCALE"
    t = m * SCALE * SCALE
    r = math.isqrt(t)
    return r if r*r == t else r+1
perimeter = sum(vmax)
squares = [round(l*l) for l in vmax]
assert all(abs(math.sqrt(m)-l) < 1e-9 for m,l in zip(squares,vmax))   # really integers
u_up_scaled = sum(sqrt_up_scaled(m) for m in squares)      # exact int, >= u_max*SCALE
getcontext().prec = 60
getcontext().rounding = ROUND_UP
perimeter_precise = Decimal(u_up_scaled) / Decimal(SCALE)   # exact (power of 10)
print(f"{perimeter=}, precise upper bound =", perimeter_precise)
sumKmax = sum(Kmax) # exact: multiples of 1/8
print("total number of refvectors:", total_unpruned_vectors)
print(f"u_max={perimeter=} {sumKmax=}")
F_const = 1/(math.tan(math.pi/num_polyominoes)*4*num_polyominoes)
# five floating-point approximations: pi, division, tan, multiplication (*4 is exact), division "1/"
F_tan = math.tan(math.pi/num_polyominoes)
# three floating-point approximations: pi, division, tan
F_u = perimeter**2 * F_const
# u^2 is computed in Decimal with 60 digits (rounded up), the float constant
# gets the relative safety factor below; the final additions are done in Decimal.
F_u_precise = perimeter_precise**2 / (Decimal(F_tan)*4*num_polyominoes)
print(f"F{num_polyominoes}(u_max)=", F_u, "precise =",F_u_precise)
print(f"F{num_polyominoes}(u_max)+K_max=", F_u + sumKmax, "precise =",F_u_precise + Decimal(sumKmax))
fudge_factor = 1+1e-15 # rather generous, instead of directed rounding
# machine rounding error ≈ 1.11e-16
fudge = Decimal(fudge_factor)
print("with fudge factor for tan and pi", F_u_precise*fudge + Decimal(sumKmax))

gc.collect()


def compute_linear_function(perimeter, K_max):
    F_u = perimeter**2 * F_const
    u_crit = math.sqrt((target-K_max)/F_const)
    F_crit = u_crit**2*F_const
    print (f"u0={u_crit}, F{num_polyominoes}(u0)={F_crit}")
    assert math.isclose(F_crit, target-K_max, rel_tol=1e-12)

    A0=(F_u-F_crit)/(perimeter-u_crit)
    B0 = F_u-perimeter*A0
    assert math.isclose(A0*u_crit+B0, F_crit, rel_tol=1e-12)
    print(f"{A0=} {B0=}")
    u_mid = (u_crit+perimeter)/2
    print("max. approximation gap of the linear function (in the middle of the range):",
          (F_crit+F_u)/2- u_mid**2*F_const)
    return A0,B0

A0, B0 = compute_linear_function(perimeter, sumKmax)
original_A0 = A0


# PASS 2: Calculate S3 and local slack
print("\nPASS 2: Calculating S3 and initial slack...")
sum_best = sum_pen = 0
ALWAYS_LONGEST = True
for P,vm in zip(Pieces, vmax):
    ref_vectors = compute_ref_vectors(P)
    best = max(math.sqrt(x**2+y**2)*A0-pen for ((x,y),rot),pen in ref_vectors.items())
    best2 = max(math.sqrt(x**2+y**2)*A0-pen for ((x,y),rot),pen in ref_vectors.items() if
                math.sqrt(x**2+y**2)==vm)
    best_pen = min(pen for ((x,y),rot),pen in ref_vectors.items() if  math.sqrt(x**2+y**2)==vm)
    sum_pen += best_pen
    sum_best += best
    if best != best2:
        ALWAYS_LONGEST = False

S3 = sum_best + B0
print(f"{S3=}. Still to go to {target}: {S3-target}")
slack = S3 - target

if ALWAYS_LONGEST:
    print("The choice was always a longest vector.")
    print("S3 alternative (for comparison) =", perimeter*A0-sum_pen+B0)
    print("S3 direct (for comparison) =", F_u-sum_pen)
    print("S3 precise:", F_u_precise - Decimal(sum_pen))
    print("S3 precise with safety:", F_u_precise*fudge - Decimal(sum_pen))
else:
    print("The choice was NOT always a longest vector. Aborting")
    exit(1)        


exit() # The remainder were some trials for other purposes, which I forgot.


gc.collect()


# PASS 3: Apply vector reduction on-the-fly and collect statistics
print("\nNOW REDUCE THE VECTORS:")
print("PASS 3: Calculating reduced parameters...")

vmax2 = []
Kmax2 = []
tot_ref_vectors = 0
totl2 = 0
sumk2 = 0

for idx, P in enumerate(Pieces):
    ref_vectors = compute_ref_vectors(P)
    best_L = vmax[idx]
    best_K = Kmax[idx]
    
    # Prune vectors on-the-fly
    ref_vectors2 = {
        ((x,y),rot): pen
        for ((x,y),rot),pen in ref_vectors.items()
        if (best_L - math.sqrt(x**2+y**2))*A0 + (best_K-(-pen)) <= slack
    }
    
    tot_ref_vectors += len(ref_vectors2)
    vmax2.append(max(math.sqrt(x**2+y**2) for (x,y),rot in ref_vectors2))
    Kmax2.append(-min(ref_vectors2.values()))
    
    l2 = sum(1 for ((x,y),rot),pen in ref_vectors2.items() if rot==0)
    totl2 += l2
    
    k2_list = [-pen for ((x,y),rot),pen in ref_vectors2.items() if rot==0]
    if k2_list:
        sumk2 += max(k2_list)

perimeter2 = sum(vmax2)
sumKmax2 = sum(Kmax2)

print("total number of refvectors:", tot_ref_vectors)
print(f"u_max2={perimeter2=} {sumKmax2=} sum(kmax with rot=0): {sumk2}")
print("total number of refvectors with rot=0:", totl2)

gc.collect()


# PASS 4: Recalculate secondary linear function and S3_new
print("\nPASS 4: Performing secondary linear calculation...")
A0, B0 = compute_linear_function(perimeter2, sumKmax2)

sum_best = 0
for idx, P in enumerate(Pieces):
    ref_vectors = compute_ref_vectors(P)
    best_L = vmax[idx]
    best_K = Kmax[idx]
    
    # Regenerate pruned set using original coefficients
    ref_vectors2 = {
        ((x,y),rot): pen
        for ((x,y),rot),pen in ref_vectors.items()
        if (best_L - math.sqrt(x**2+y**2))*original_A0 + (best_K-(-pen)) <= slack
    }
    
    best = max(math.sqrt(x**2+y**2)*A0-pen for ((x,y),rot),pen in ref_vectors2.items())
    sum_best += best

S3_new = sum_best + B0
print(f"{S3_new=}. Still to go to {target}: {S3_new-target}")
