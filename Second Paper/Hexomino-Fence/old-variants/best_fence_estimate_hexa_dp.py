"""
Shared code of Step 3 and Step 4 (not run directly, it is imported by
best-fence-estimate-hexa-step3.py, best-fence-estimate-hexa-step4.py
and best-fence-estimate-hexa-step4-count.py).

Contents
  * load_step_result    reads the result file of an earlier step
  * solve_vector_set    runs the dynamic program on one vector set (forward pass)
  * reconstruct_one / reconstruct_all_lists / reconstruct_dag
                        three ways to reconstruct reference polygons from the
                        tables of the forward pass (used by Step 4)
The dynamic-programming code itself is unchanged with respect to the former
all-in-one Step 3 program.
"""
import sys
import importlib
import itertools
import math
from collections import defaultdict
from math import sqrt,atan2
from types import SimpleNamespace

# Dynamic programming from some given list of vectors plus penalty

# Pruning by area that can be reached with the current four sides (stored as "key")
# and the remaining vectors, forming the largest polygon with these edge lengths.
# That is: Accept only if the value "val" is decent enough.

# No artificial or heuristic pruning of the search space. Yields the optimal
# solution with the given reference vectors, provided "target_area"
# is not too high.

# For improvement of runtime (fewer calls of pruning function):
# We check the pruning criterium only AFTER the table is formed,
# and throw out the entries that fail.
# (Some memory could be saved by checking the pruning criterium with every insertion.)

# Possible variation or addition: use the Zauberformel (Minkowski-mean)
# This has been tried, but it yields very poor bounds.

# A weak form of meet-in-the middle:
# Compute reachable vector sums from the end (a quite moderate number at most W**2)
# Reduction of reachable set to first octant is not done, because the set is small anyway.

# In addition: rotate vector sum in solution into first
# octant to normalize (saves a factor 8 in most cases)
# The reconstruction procedure of the solution has been adapted accordingly.

# Another idea. Throw away the tables more than, say, 10 steps back.
# saves space. Aber dann müsste man sich die Lösung in mehreren Läufen
# rekonstruieren, und Schritt für Schritt zurückhangeln.

# Another idea: Push the reference vectors even more INSIDE the area.
# Have to think carefully whether this is valid.
# This may kill some more options when the corner pieces are not exactly diagonal.
# Empirically, it has no influence when all four corner vectors are
# equal (and equal to (2.5,2.5).)

# state used by the pruning functions (set by run_forward_dp / forward_dp)
angle_area_store = {}
area_threshold = None
sum_fixed_vectors = None
sum_x_flat = None
bin_search_calls = bin_search_decisions = 0


def step_module_name(step, n_max, goal):
    return f"best_fence_candidate_solutions_step_{step}_n{n_max}t{goal}"

def load_step_result(step, n_max, goal):
    """import the result file of Step 2 or Step 3; returns the module (its list is .RESULT)"""
    name = step_module_name(step, n_max, goal)
    try:
        return importlib.import_module(name)
    except ModuleNotFoundError as e:
        if e.name != name:
            raise
        sys.exit(f'Input file "{name}.py" not found. Run Steps 1 to {step} first with the same '
                 f'arguments (n_max={n_max}, goal={goal}) and run this program from the directory '
                 f'that contains the result files and the scripts.')


def rotateL(P):
    return [(-y,x) for x,y in P]
def rotateR(P):
    return [(y,-x) for x,y in P]
def rotate2(P):
    return [(-x,-y) for x,y in P]
def rotate0(P): # identity
    return [(x,y) for x,y in P]
def reflect(P): # reflect at x=y line
    return [(y,x) for x,y in P]

#### 2x2 matrices (a11,a12,a21,a22)
IDENTITY = (1,0,0,1)
REFL_diag = (0,1,1,0)
REFL_hori = (-1,0,0,1)
REFL_vert = (1,0,0,-1)
def matmul(A,B):
    a11,a12,a21,a22 = A
    b11,b12,b21,b22 = B
    return (a11*b11+a12*b21, a11*b12+a12*b22, a21*b11+a22*b21, a21*b12+a22*b22)
def apply_mat(A,x): # matrix times vector A*x
    a11,a12,a21,a22 = A
    x1,x2 = x
    return (a11*x1+a12*x2, a21*x1+a22*x2)

def det(a,b): # signed area of parallelogram spanned by two vectors
    (x1,y1),(x2,y2) = a,b
    return x1*y2-y1*x2

ALL_SYMMETRIES = list(itertools.product((rotate0,reflect),
                                        (rotate0,rotateL,rotate2,rotateR)))

def find_min_radius(vectors, eps=1e-9):
    """radius of a circle in which a polygon with
    edge lengths 'vectors' can be inscribed"""
    perimeter = sum(vectors)
    R_low = perimeter/(2*math.pi)
    R_high = perimeter/4
    while R_high-R_low > eps: # binary search
        R = (R_high+R_low)/2
        angle = sum(half_angle(R,v) for v in vectors)
        if angle>=math.pi: # R too small
            R_low = R
        else: # R too big
            R_high = R
    return R

def prune_solution(fixed_vectors,variable_vectors, additional_area, eps=1e-4):
    # the area unit for additional_area is 1.
    #global area_threshold, R_min, sum_fixed_vectors
    global bin_search_calls, bin_search_decisions
    bin_search_calls += 1
    
    # for the longest edge, the central angle might conceivably be bigger than pimeter.
    perimeter = sum_fixed_vectors + sum(variable_vectors)
    R_LOWER = perimeter/(2*math.pi)
    R_UPPER = perimeter/4.01 # bound on reasonable values

    R_2 = max(variable_vectors)/2
    if R_2>R_LOWER: # extra long vector exists:
        if sum(half_angle(R_2,v) for v in fixed_vectors+variable_vectors) <= math.pi:
            # overstretched: center of circle lies outside the inscribed polygon
            # Must be a weird situation
            assert R_2**2*math.pi/2 <=area_threshold # area of semicircle
            return True # prune
        else:
            R_LOWER = R_2
    R_low  = math.floor(R_LOWER)
    R_high= R_low+32
    # integers (or multiples of 1/16 or the like) so that caching gets a chance!
    
    while R_high-R_low > eps: # binary search
        R = (R_high+R_low)/2
        if R<R_LOWER:
            R_low = R
            continue
        if R>R_UPPER:
            R_high = R
            continue
        angle,area = angle_and_area(R,fixed_vectors)
        angle += sum(half_angle(R,v) for v in variable_vectors)
        area  += sum(triangle(R,v)   for v in variable_vectors) + additional_area
        bin_search_decisions += 1
        if angle<=math.pi: # R too big
            if area < area_threshold-1e-9:
                return True # prune: this cannot possibly lead to an optimal solution
            R_high = R
        else: # R too small
            if area>=area_threshold:
                return False # don't prune: this might still lead to an optimal solution
            R_low = R
    return False

"""As the computation proceeds, the radius R will get bigger"""

def angle_and_area(R,fixed_vectors):
    # fixed_vectors is supposed to remain unchanged between calls
    # while the cache angle_area_store is not reset.
    try:
        return angle_area_store[R]
    except KeyError:
        pass
    angle = sum(half_angle(R,v) for v in fixed_vectors)
    area = sum(triangle(R,v) for v in fixed_vectors)
    angle_area_store[R] = angle,area
    return angle,area
    
            
def triangle(R,s):
    "area of a triangle spanned by a chord of length s in a circle of radius R"
    return s/2 * sqrt(R**2-s**2/4)

def half_angle(R,s):
    "half the central angle spanned by a chord of length s in a circle of radius R"
    return math.asin(s/(2*R))

def combine(u,v):
    "mixed product L1-norm x Lmax-norm"
    ux,uy=map(abs,u)
    vx,vy=map(abs,v)
    return max(max(ux,uy)*(vx+vy), (ux+uy)*max(vx,vy))/16

def ZF_single(vectors): # Zauberformel, (no penalties)
    return sum(combine(p,p2)
               for p in vectors
               for p2 in vectors)

def ZF_value(vectors):
    return sum_x_flat*sum(abs(x)+abs(y) for x,y in vectors)/8 + ZF_single(vectors)




def forward_dp(Vecs,Corners,keep_tables):
    """The workhorse of this algorithm.
    Vecs (x,y) satisfy 0<=y<=x and are sorted by slope."""
    # special treatment of corners
    
    global angle_area_store, area_threshold, sum_fixed_vectors, sum_x_flat

    num = len(Vecs)

    # try all assignments of corners to the four directions, up to symmetry
    Start_configurations = set()
    for refls4 in itertools.product((rotate0,reflect),repeat=4):
        for c4 in itertools.permutations(Corners):
            [v1],[v2],[v3],[v4] = [ref([c]) for ref,c in zip(refls4,c4)]
            (x1,y1)=v1
            (x2,y2)=v2; x2=-x2 # rotate the corner vectors into the 4 quadrants.
            (x3,y3)=v3; x3=-x3; y3=-y3
            (x4,y4)=v4; y4=-y4
            P2 = (x1,y1),(x2,y2),(x3,y3),(x4,y4)
            Q = frozenset(frozenset(r1(r2(P2))) for r1,r2 in ALL_SYMMETRIES)
            Start_configurations.add(Q)
    print(len(Start_configurations), "start configurations of corners:")

    TABLES = dict() # save the dictionaries to allow backtracing the solutions
    Table = TABLES[num] = dict()
    
    for Q in Start_configurations:
        Q0 = list(Q)[0]
        print(*Q0)
        if Q0=={(0,0)}:
            Table[(0,0),(0,0),(0,0),(0,0)] = 0
        else: # split into the four quadrants:
            [(x1,y1)] = [(x,y) for x,y in Q0 if x>0 and y>0]
            [(x2,y2)] = [(x,y) for x,y in Q0 if x<0 and y>0]
            [(x3,y3)] = [(x,y) for x,y in Q0 if x<0 and y<0]
            [(x4,y4)] = [(x,y) for x,y in Q0 if x>0 and y<0]
            Table[(x1,y1),(x2,y2),(x3,y3),(x4,y4)] = 0
        # objective value of the initial solution, consisting only of four
        # single vectors.
    print("start table:", Table)
    angle_area_store = dict()
    def insert(key,val):
        (x1,y1),(x2,y2),(x3,y3),(x4,y4) = key
        sx,sy = x1+x2+x3+x4, y1+y2+y3+y4
        if sx<0:
            (x1,y1),(x2,y2),(x3,y3),(x4,y4) = key = (-x2,y2),(-x1,y1),(-x4,y4),(-x3,y3)
            sx = -sx
        if sy<0:
            (x1,y1),(x2,y2),(x3,y3),(x4,y4) = key = (x4,-y4),(x3,-y3),(x2,-y2),(x1,-y1)
            sy = -sy
        if sy>sx:
            # (x1,y1),(x2,y2),(x3,y3),(x4,y4) =
            key = (y1,x1),(y4,x4),(y3,x3),(y2,x2)
        if key not in newTable or val>newTable[key]:
            newTable[key] = val

    # a very weak form of meet-in-the middle:
    # compute reachable vectors from the end
    NUM_FROM_END = min((num)//2+6, num) # how many vectors are combined, at most
    reachablevectors = {0: {(0,0)}}
    # reachablevectors[i] combines the last i vectors Vecs[0]...Vecs[i-1]
    # in all possible ways.
    # Compute reachablevectors[1] ... reachablevectors[NUM_FROM_END]
    for i in range(NUM_FROM_END):
        (x0,y0) = Vecs[i]
        all_directions = set((z for (x1,y1) in ((x0,y0),(x0,-y0),(-x0,y0),(-x0,-y0))
                              for z in ((x1,y1),(y1,x1))))
        reachablevectors[i+1] = set(((x+u,y+v) for (x,y) in reachablevectors[i]
                                     for (u,v) in all_directions
                                     ))
        print(f"len(reachablevectors[{i+1}])={len(reachablevectors[i+1])}")

    # START THE MAIN LOOP
    for i in reversed(range(num)):
        (x0,y0) = Vecs[i]
        fixed_vectors = [sqrt(x**2+y**2) for x,y in Vecs[:i]] # lengths
        sum_fixed_vectors = sum(fixed_vectors)
        angle_area_store = dict() # clear the cache
        if i==NUM_FROM_END:
            print("reachability information kicks in.")
        if i<=NUM_FROM_END:
            # redefine "insert". from now on: check reachability from the end
            reachable = reachablevectors[i]
            def insert(key,val): 
                (x1,y1),(x2,y2),(x3,y3),(x4,y4) = key
                sx,sy = x1+x2+x3+x4, y1+y2+y3+y4
                if (sx, sy) not in reachable:
                    return # throw away entries that cannot be completed
                if sx<0:
                    (x1,y1),(x2,y2),(x3,y3),(x4,y4) = key = (-x2,y2),(-x1,y1),(-x4,y4),(-x3,y3)
                    sx = -sx
                if sy<0:
                    (x1,y1),(x2,y2),(x3,y3),(x4,y4) = key = (x4,-y4),(x3,-y3),(x2,-y2),(x1,-y1)
                    sy = -sy
                if sy>sx:
                    # (x1,y1),(x2,y2),(x3,y3),(x4,y4) =
                    key = (y1,x1),(y4,x4),(y3,x3),(y2,x2)
                if key not in newTable or val>newTable[key]:
                    newTable[key] = val
        newTable = dict()
        for (((x1,y1),(x2,y2),(x3,y3),(x4,y4)),vol) in Table.items():
            # order of product is ccw from positive x-axis
            # first octant, bundle in the direction to the right:
            (x,y) = (x0,y0)
            insert(((x1+x,y1+y),(x2,y2),(x3,y3),(x4,y4)), vol+det((x,y),(x1,y1)))
            # second octant, bundle in the upward direction:
            if y0!=x0 and y0!=0: 
                # if y0==x0, the vector could eqally go to the second or first octant.
                # avoid unnecessary blow-up in this case by assigning to the first octant only.
                # if y0==0, the vector could eqally go to the second or third octant.
                # avoid unnecessary blow-up in this case by assigning to the third octant only.
                (x,y) = (y0,x0)
                insert(((x1+x,y1+y),(x2,y2),(x3,y3),(x4,y4)), vol+det((x1,y1),(x,y)))
            # second quadrant
            (x,y) = (-y0,x0) # bundle in the upward direction
            insert(((x1,y1),(x2+x,y2+y),(x3,y3),(x4,y4)) ,vol+det((x,y),(x2,y2)))
            if y0!=x0 and y0!=0: 
                (x,y) = (-x0,y0) # bundle in the direction to the left
                insert(((x1,y1),(x2+x,y2+y),(x3,y3),(x4,y4)), vol+det((x2,y2),(x,y)))
            # third quadrant
            (x,y) = (-x0,-y0)
            insert(((x1,y1),(x2,y2),(x3+x,y3+y),(x4,y4)), vol+det((x,y),(x3,y3)))
            if y0!=x0 and y0!=0: 
                (x,y) = (-y0,-x0)
                insert(((x1,y1),(x2,y2),(x3+x,y3+y),(x4,y4)), vol+det((x3,y3),(x,y)))
            # fourth quadrant
            (x,y) = (y0,-x0)
            insert(((x1,y1),(x2,y2),(x3,y3),(x4+x,y4+y)),vol+det((x,y),(x4,y4)))
            # last octant
            if y0!=x0 and y0!=0: 
                (x,y) = (x0,-y0)
                insert(((x1,y1),(x2,y2),(x3,y3),(x4+x,y4+y)), vol+det((x4,y4),(x,y)))
        # do the pruning only AFTER all insertions:
        beforepruning = len(newTable)
        print(f"iteration {i=} piece={(x0,y0)}, len(Table) before pruning: {beforepruning:,}", flush=True)

        # pruning with Minkowski-mean:
        #print(f"iteration {i=} piece={(x0,y0)}, {len(newTable)=:,}, before Minkowski-pruning", flush=True)
        #ZF_flat = ZF_single(Vecs[:i])
        #sum_x_flat = sum(x for x,y in Vecs[:i])
        ##for key in newTable: # checking
        ##    assert ZF_flat + ZF_value(key) == ZF_single(Vecs[:i]+list(key))
        #newTable = {key:val for key,val in newTable.items()
        #            if ZF_flat + ZF_value(key) +  val/2 > area_threshold}
        #print(f"iteration {i=} piece={(x0,y0)}, {len(newTable)=:,}, after Minkowski-pruning", flush=True)

        newTable = {key:val for key,val in newTable.items()
                    if not prune_solution(fixed_vectors,[sqrt(x**2+y**2) for x,y in key], val/2)}
        TABLES [i] = Table = newTable
        if not keep_tables:
            del TABLES[i+1] # release space
        print(f"iteration {i=} piece={(x0,y0)}, {len(Table)=:,} (before pruning: {beforepruning:,})", flush=True)
    return TABLES


def run_forward_dp(Vecs, Corners, threshold, keep_tables):
    """Forward pass of the dynamic program.
    threshold: the value area_threshold = target_area + penalty used for pruning.
    keep_tables: False = only the optimal value is needed (saves space, the tables of
                 the previous stages are thrown away); True = keep all tables
                 for the reconstruction of the solutions.
    Returns (best, TABLES, finals); best is None if there is no solution at all.
    finals = list of (key, value) of the optimal final states."""
    global area_threshold, bin_search_calls, bin_search_decisions
    area_threshold = threshold
    bin_search_calls = bin_search_decisions = 0
    TABLES = forward_dp(Vecs, Corners, keep_tables)
    print(f"Total number of binary search calls = {bin_search_calls:,},",
          f"average number of decisions per call =",
          f"{bin_search_decisions/max(bin_search_calls,1):.5}")
    Table = TABLES[0]
    best = max((vol + det((x1,y1),(x2,y2)) + det((x3,y3),(x4,y4))
               # det(v1+v2,v3+v4)==0 because the two vectors are opposite.
               # Therefore only det(v1,v2)+det(v3,v4) remains to be added
               for (((x1,y1),(x2,y2),(x3,y3),(x4,y4)), vol) in Table.items()
               if (x1+x4+x3+x2, y1+y4+y3+y2)==(0,0)), default=None)
    if best is None:
        return None, TABLES, [] # (the caller reports this)
    numbest = 0
    finals = []
    for ((x1,y1),(x2,y2),(x3,y3),(x4,y4)),vol in Table.items():
        if (x1+x4+x3+x2, y1+y4+y3+y2)==(0,0) and vol + det((x1,y1),(x2,y2)) + det((x3,y3),(x4,y4))==best:
            bestsol = (x1,y1),(x2,y2),(x3,y3),(x4,y4)
            numbest += 1
            if numbest<=30:
                print(f"A best solution ({numbest}): {bestsol}")
            finals.append((bestsol,vol))
    if numbest>30:
        print("...")
    print(f"{numbest}  best solutions in total")
    print(best,"best value", flush=True)
    return best, TABLES, finals


def prepare_instance(nline, vecs, neg_penalty):
    """Turn one entry (vecs, neg_penalty) of the Step 2 file into the input of the
    dynamic program. Returns Vecs, Corners, shift, penalty."""
    penalty = -neg_penalty
    assert penalty>0
    print("\n======================")
    print(f"Vector set #{nline}:",vecs)
    print(f"{penalty=}")
    Vecs = []
    Corners = []
    for (x,y),m in vecs:
        # (Don't scale Vectors by 2 to make them integral)
        y,x=sorted((abs(x),abs(y))) # now 0<=y<=x, first octant
        if x==int(x):
            Vecs += [(y/x,(x,y))]*m
        else:
            Corners += [(x,y)]*m
        continue
    assert len(Corners)==4 # corners need not be sorted by slope.
    # for shift in -2.5,-1.5,-0.5,+0.5,+1.5,+2.5:
        # Perhaps one could also try -1.5,-2.5? Is this legal?
    for shift in -0.5,+0.5,+1.5,+2.5,3.5,4.5:
        # first try to move the reference vectors to the INNER endpoint;
        # only if necessary, move them outside.
        # "inside" seems to lead to stronger bounds.
        Corners_new = [(x+shift,y+shift) for x,y in Corners] #move the reference point
        if (all((x,y)==(0,0) for (x,y) in Corners_new) or
            (all(x>0 and y>0 for (x,y) in Corners_new) and
             max(s for (s,p) in Vecs)
                <= min(y/x for (x,y) in Corners_new))):
            # corners are the steepest
            print(f"shifting reference points by {shift}",
                  "(positive=outwards, negative=inwards).")
            break
    else:
        print("Could not make the corner pieces the steepest pieces.")

    for (x,y) in Corners:
        penalty += (x+shift)*(y+shift)/2 - x*y/2
    for s,(x,y) in Vecs:
        penalty += x*shift # x is the larger coordinate
        
        
    Corners = Corners_new
    print(f"adusted {penalty=}, after moving the reference vectors.")
    Vecs = [p for (s,p) in sorted(Vecs)] # steepest last
    assert all(c==int(c) for p in Vecs+Corners for c in p)
    Vecs    = [(int(x),int(y)) for x,y in Vecs]
    Corners = [(int(x),int(y)) for x,y in Corners]
    
    """ previously: (scaled by 2)
    one exceptional case, where one cannot start with the corner pieces:
    print(Vecs[-7:])
        [.... (8, 4), (5, 3), (6, 4), (6, 4), (5, 5), (5, 5), (5, 5)]
    SOLUTION: take the OUTER point as the reference point.
    Then 5/3 become 6/4 (and no need for scaling)
    """
    
    total_parity_sum = sum(x+y for (x,y) in Vecs+Corners)
    print("total L1-distance must be even:", total_parity_sum)
    assert total_parity_sum%2 == 0 # otherwise solution is impossible.
    return Vecs, Corners, shift, penalty


def solve_vector_set(nline, vecs, neg_penalty, target_area, keep_tables):
    """Prepare one vector set and run the forward pass. Returns a namespace with
    Vecs, Corners, shift, penalty, best, TABLES, finals,
    enclosed_area (None if there is no construction at all) and
    reaches_goal (True if enclosed_area > target_area, i.e. the result is guaranteed optimal)."""
    Vecs, Corners, shift, penalty = prepare_instance(nline, vecs, neg_penalty)
    best, TABLES, finals = run_forward_dp(Vecs, Corners, target_area + penalty, keep_tables)
    enclosed_area = None
    reaches_goal = False
    if best is None:
        print(f"No solution with area at least {target_area} was found.")
    else:
        B = best
        print(f"best solution value {B/2}. **** enclosed area={B/2-penalty=} ************")
        enclosed_area = B/2-penalty
        if enclosed_area<=target_area:
            print(f"This is not guaranteed to be the optimum, because it is < {target_area=}")
        else:
            reaches_goal = True
            print(f"This is guaranteed to be the optimum.")
        print(flush=True)
    return SimpleNamespace(Vecs=Vecs, Corners=Corners, shift=shift, penalty=penalty,
                           best=best, TABLES=TABLES, finals=finals,
                           enclosed_area=enclosed_area, reaches_goal=reaches_goal)


# ---------------------------------------------------------------------------
# Reconstruction of the solutions (Step 4)
# ---------------------------------------------------------------------------

def reconstruct_one(Vecs, TABLES, bestsol, bestvol):
    """One optimal reference polygon (list of vectors, corner vectors last,
    still with the shifted corner vectors)."""
    num = len(Vecs)
    # Now the optimum value has been determined.
    # We identify the solution by tracking backward:
    SOL = [] # list of vectors used in the solution (in addition to the corner vectors)
    MAT = IDENTITY
    def find(key,vol):
        nonlocal MAT # also: x,y, bestsol, bestvol
        if vol!=bestvol:
            return False
        # normalize the sum into the first octant
        # remember: start table = {((7, 7), (-7, 7), (-7, -7), (7, -7)): 196}
        (x1,y1),(x2,y2),(x3,y3),(x4,y4) = key
        sx,sy = x1+x2+x3+x4, y1+y2+y3+y4
        if sx<0:
            (x1,y1),(x2,y2),(x3,y3),(x4,y4) = key = (-x2,y2),(-x1,y1),(-x4,y4),(-x3,y3)
        if sy<0:
            (x1,y1),(x2,y2),(x3,y3),(x4,y4) = key = (x4,-y4),(x3,-y3),(x2,-y2),(x1,-y1)
        if abs(sy)>abs(sx):
            # (x1,y1),(x2,y2),(x3,y3),(x4,y4) =
            key = (y1,x1),(y4,x4),(y3,x3),(y2,x2)
        if key!=bestsol:
            return False
        # All future solution pieces must be subjected to the same operations, in reverse order
        if abs(sy)>abs(sx):
            MAT = matmul(MAT,REFL_diag)
        if sy<0:
            MAT = matmul(MAT,REFL_vert)
        if sx<0:
            MAT = matmul(MAT,REFL_hori)
        piece = apply_mat(MAT,(x,y))
        print(f"piece used {piece}")
        SOL.append(piece)
        return True
    
    for i in range(num):
        (x0,y0) = Vecs[i]
        for oldkey,oldvol in TABLES[i+1].items():
            (x1,y1),(x2,y2),(x3,y3),(x4,y4) = oldkey
            # order of product is ccw from positive x-axis
            # first octant
            (x,y) = (x0,y0)
            if find(((x1+x,y1+y),(x2,y2),(x3,y3),(x4,y4)),oldvol+det((x,y),(x1,y1))): break
            # second octant
            (x,y) = (y0,x0)
            if find(((x1+x,y1+y),(x2,y2),(x3,y3),(x4,y4)),oldvol+det((x1,y1),(x,y))): break
            # second quadrant
            (x,y) = (-y0,x0)
            if find(((x1,y1),(x2+x,y2+y),(x3,y3),(x4,y4)),oldvol+det((x,y),(x2,y2))): break
            (x,y) = (-x0,y0)
            if find(((x1,y1),(x2+x,y2+y),(x3,y3),(x4,y4)),oldvol+det((x2,y2),(x,y))): break
            # third quadrant
            (x,y) = (-x0,-y0)
            if find(((x1,y1),(x2,y2),(x3+x,y3+y),(x4,y4)),oldvol+det((x,y),(x3,y3))): break
            (x,y) = (-y0,-x0)
            if find(((x1,y1),(x2,y2),(x3+x,y3+y),(x4,y4)),oldvol+det((x3,y3),(x,y))): break
            # fourth quadrant
            (x,y) = (y0,-x0)
            if find(((x1,y1),(x2,y2),(x3,y3),(x4+x,y4+y)),oldvol+det((x,y),(x4,y4))): break
            # last octant
            (x,y) = (x0,-y0)
            if find(((x1,y1),(x2,y2),(x3,y3),(x4+x,y4+y)),oldvol+det((x4,y4),(x,y))): break
        else:
            assert 0 # not found
        bestsol, bestvol = oldkey, oldvol
    return SOL + [apply_mat(MAT,v) for v in bestsol] # add the corner vectors


def reconstruct_all_lists(Vecs, TABLES, finals):
    """All optimal reference polygons, by carrying along the list of all partial
    solutions. Needs a lot of memory. Returns a list of vector lists
    (corner vectors last, still with the shifted corner vectors)."""
    num = len(Vecs)
    SOLS_new = defaultdict(list) # for finding all solutions
    # SOLS_new, SOLS_old: dictionary (key,value) -> [list of (partial solution, MAT)]
    for bestsol,vol in finals:
        SOLS_new[bestsol,vol].append(([],IDENTITY))
    def find2(key,vol): # variation for finding all solutions
        # nonlocal oldkey, oldvol, SOLS_old, SOLS_new, x, y
        (x1,y1),(x2,y2),(x3,y3),(x4,y4) = key
        sx,sy = x1+x2+x3+x4, y1+y2+y3+y4
        if sx<0:
            (x1,y1),(x2,y2),(x3,y3),(x4,y4) = key = (-x2,y2),(-x1,y1),(-x4,y4),(-x3,y3)
        if sy<0:
            (x1,y1),(x2,y2),(x3,y3),(x4,y4) = key = (x4,-y4),(x3,-y3),(x2,-y2),(x1,-y1)
        if abs(sy)>abs(sx):
            # (x1,y1),(x2,y2),(x3,y3),(x4,y4) =
            key = (y1,x1),(y4,x4),(y3,x3),(y2,x2)
        if (key,vol) not in SOLS_old:
            return
        for part_SOL,MAT in SOLS_old[key,vol]:
            # All future solution pieces must be subjected to the same operations, in reverse order
            if abs(sy)>abs(sx):
                MAT = matmul(MAT,REFL_diag)
            if sy<0:
                MAT = matmul(MAT,REFL_vert)
            if sx<0:
                MAT = matmul(MAT,REFL_hori)
            SOLS_new[oldkey,oldvol].append((part_SOL+[ apply_mat(MAT,(x,y))],MAT))
        
    for i in range(num):
        (x0,y0) = Vecs[i]
        SOLS_old = SOLS_new
        print("piece",i, "#entries =", len(SOLS_old),
              "#partial solutions =", sum(len(z) for z in SOLS_old.values()), flush=True)
        #print(SOLS_old)
        SOLS_new = defaultdict(list)
        for oldkey,oldvol in TABLES[i+1].items():
            (x1,y1),(x2,y2),(x3,y3),(x4,y4) = oldkey
            (x,y) = (x0,y0)
            find2(((x1+x,y1+y),(x2,y2),(x3,y3),(x4,y4)),oldvol+det((x,y),(x1,y1)))
            (x,y) = (y0,x0)
            find2(((x1+x,y1+y),(x2,y2),(x3,y3),(x4,y4)),oldvol+det((x1,y1),(x,y)))
            (x,y) = (-y0,x0)
            find2(((x1,y1),(x2+x,y2+y),(x3,y3),(x4,y4)),oldvol+det((x,y),(x2,y2)))
            (x,y) = (-x0,y0)
            find2(((x1,y1),(x2+x,y2+y),(x3,y3),(x4,y4)),oldvol+det((x2,y2),(x,y)))
            (x,y) = (-x0,-y0)
            find2(((x1,y1),(x2,y2),(x3+x,y3+y),(x4,y4)),oldvol+det((x,y),(x3,y3)))
            (x,y) = (-y0,-x0)
            find2(((x1,y1),(x2,y2),(x3+x,y3+y),(x4,y4)),oldvol+det((x3,y3),(x,y)))
            (x,y) = (y0,-x0)
            find2(((x1,y1),(x2,y2),(x3,y3),(x4+x,y4+y)),oldvol+det((x,y),(x4,y4)))
            (x,y) = (x0,-y0)
            find2(((x1,y1),(x2,y2),(x3,y3),(x4+x,y4+y)),oldvol+det((x4,y4),(x,y)))
    return [SOL + [apply_mat(MAT,v) for v in bestsol]
            for (bestsol,_),l in SOLS_new.items() for SOL,MAT in l]


def reconstruct_dag(Vecs, TABLES, finals):
    """All optimal reference polygons as a DAG (saves space, allows counting).
    SOLS[i]: dictionary (key,value) -> (number of partial solutions, [list of (successor node, MAT, piece)])
    MAT = only the matrix for transition between this stage and the next, not the accumulated matrix.
    In the the MAT's have to be applied in one sweep from right to left.
    successor node == (key,value) in SOLS[i+1]
    (This "value" could be moved from the key of the dictionary to the values.)
    Returns SOLS; use backtrack_solutions to walk through the paths."""
    num = len(Vecs)
    SOLS = dict()
    SOLS_new = SOLS[0] = dict() # for finding all solutions
    for bestsol,vol in finals:
        SOLS_new[bestsol,vol]=(1, None)
    def find2(key,vol): # variation for finding all solutions
        # nonlocal oldkey, oldvol, SOLS_old, SOLS_new, x, y
        (x1,y1),(x2,y2),(x3,y3),(x4,y4) = key
        sx,sy = x1+x2+x3+x4, y1+y2+y3+y4
        if sx<0:
            (x1,y1),(x2,y2),(x3,y3),(x4,y4) = key = (-x2,y2),(-x1,y1),(-x4,y4),(-x3,y3)
        if sy<0:
            (x1,y1),(x2,y2),(x3,y3),(x4,y4) = key = (x4,-y4),(x3,-y3),(x2,-y2),(x1,-y1)
        if abs(sy)>abs(sx):
            # (x1,y1),(x2,y2),(x3,y3),(x4,y4) =
            key = (y1,x1),(y4,x4),(y3,x3),(y2,x2)
        try:
            succ_num,_ = SOLS_old[key,vol]
        except KeyError:
            return
        num_accumulated,succ_list =  SOLS_new.setdefault((oldkey,oldvol),(0,[]))
        # for succ_node,MAT in SOLS_old[key,vol]:
        # All future solution pieces must be subjected to the same operations, in reverse order
        MAT = IDENTITY
        if abs(sy)>abs(sx):
            MAT = matmul(MAT,REFL_diag)
        if sy<0:
            MAT = matmul(MAT,REFL_vert)
        if sx<0:
            MAT = matmul(MAT,REFL_hori)
        succ_list.append(((key,vol),MAT,(x,y)))
        SOLS_new[oldkey,oldvol] = (num_accumulated+succ_num,succ_list)
        
    for i in range(num):
        (x0,y0) = Vecs[i]
        SOLS_old = SOLS[i]
        print("piece",i, " #entries =", len(SOLS_old),
              f"#partial solutions = {sum(numb for numb,_ in SOLS_old.values()):,}", flush=True)
        #print(SOLS_old)
        SOLS_new = SOLS[i+1] = dict()
        for oldkey,oldvol in TABLES[i+1].items():
            (x1,y1),(x2,y2),(x3,y3),(x4,y4) = oldkey
            (x,y) = (x0,y0)
            find2(((x1+x,y1+y),(x2,y2),(x3,y3),(x4,y4)),oldvol+det((x,y),(x1,y1)))
            if y0!=x0 and y0!=0: # these tests are not strictly necessary. hope to get a more systematic count, taking into account multiplicities
                (x,y) = (y0,x0)
                find2(((x1+x,y1+y),(x2,y2),(x3,y3),(x4,y4)),oldvol+det((x1,y1),(x,y)))
            (x,y) = (-y0,x0)
            find2(((x1,y1),(x2+x,y2+y),(x3,y3),(x4,y4)),oldvol+det((x,y),(x2,y2)))
            if y0!=x0 and y0!=0: 
                (x,y) = (-x0,y0)
                find2(((x1,y1),(x2+x,y2+y),(x3,y3),(x4,y4)),oldvol+det((x2,y2),(x,y)))
            (x,y) = (-x0,-y0)
            find2(((x1,y1),(x2,y2),(x3+x,y3+y),(x4,y4)),oldvol+det((x,y),(x3,y3)))
            if y0!=x0 and y0!=0: 
                (x,y) = (-y0,-x0)
                find2(((x1,y1),(x2,y2),(x3+x,y3+y),(x4,y4)),oldvol+det((x3,y3),(x,y)))
            (x,y) = (y0,-x0)
            find2(((x1,y1),(x2,y2),(x3,y3),(x4+x,y4+y)),oldvol+det((x,y),(x4,y4)))
            if y0!=x0 and y0!=0: 
                (x,y) = (x0,-y0)
                find2(((x1,y1),(x2,y2),(x3,y3),(x4+x,y4+y)),oldvol+det((x4,y4),(x,y)))
    return SOLS

def backtrack_solutions(SOLS,level,succlist,prefix=[]):
    for (succ,MAT,piece) in succlist:
        if level==0:
            #assert number==1
            yield prefix+[MAT,piece]
        else:
            number,succlist2 = SOLS[level][succ]
            yield from backtrack_solutions(SOLS,level-1,succlist2,prefix+[MAT,piece])


def restore_corners(SOL, Corners, shift):
    """undo the shift of the reference points for the four corner vectors at the end of SOL (in place)"""
    # fix the corner vectors:
    if set(Corners)=={(0,0)}: # this is currently not possible
        SOL[-4:0] = [(-shift,-shift),(-shift,shift),(shift,-shift),(shift,shift)]
    else:
        for i in range(-4,0):
            x,y = SOL[i]
            x -= shift if x>0 else -shift
            y -= shift if y>0 else -shift
            SOL[i] = x,y

def canonical(SOL):
    """representative of a solution under rotations and reflections"""
    return tuple(min(sorted(r1(r2(SOL))) for r1,r2 in ALL_SYMMETRIES))

def enumerate_dag_solutions(SOLS, num, shift, max_solutions):
    """Count the solutions in the DAG built by reconstruct_dag, generate them
    (at most max_solutions different ones) and print them."""
    num_sols = sum(numb for numb,_ in SOLS[num].values())
    print(num_sols,"raw solutions (not taking into account symmetries).")
    ALL_SOLS = set()
    ALL_SOLS_LIST = list()
    multi = defaultdict(int)
    #print(SOLS[num])
    #print(SOLS[num-1])
    jj = 0
    for (corn,_),(_,succlist) in SOLS[num].items():
        corn2 = []
        for x,y in corn:                    
            x -= shift if x>0 else -shift
            y -= shift if y>0 else -shift
            corn2.append((x,y))
        for sol_raw in (backtrack_solutions(SOLS,num-1,succlist,[])):
                jj +=1
                MAT = IDENTITY
                #print(jj,corn,sol_raw)
                sol = []
                for kk in range(1,len(sol_raw),2):
                    piece = sol_raw[-kk]
                    MAT2 = sol_raw[-kk-1]
                    MAT = matmul(MAT,MAT2)
                    sol.append(apply_mat(MAT,piece))
                sol += [apply_mat(MAT,piece) for piece in corn2]
                SOL_Q = tuple(min(sorted(r1(r2(sol))) for r1,r2 in ALL_SYMMETRIES))
                SS = [(x,y) for _,x,y in sorted((atan2(y,x),x,y) for x,y in SOL_Q)]
                if SOL_Q not in ALL_SOLS:
                    # print the solutions as they are produced
                    #print(len(ALL_SOLS)+1,",",SS,",", flush=True)
                    ALL_SOLS.add(SOL_Q)
                    ALL_SOLS_LIST.append(SOL_Q)
                ind = ALL_SOLS_LIST.index(SOL_Q)
                multi[ind] += 1
                if jj<=50: print(jj,f"({ind})",sol)
                if jj==51: print("...\n")
                if len(ALL_SOLS)>=max_solutions:
                    print("aborting after", len(ALL_SOLS), "solutions were found.")
                    break
        if len(ALL_SOLS)>=max_solutions:
            break
    assert len(ALL_SOLS)>=max_solutions or jj==num_sols
    print("multiplicities",dict(multi))
    print("Solutions:") # print the solutions in the end
    for i,S in enumerate(ALL_SOLS_LIST,start=0):
        SS = [(x,y) for _,x,y in sorted((atan2(y,x),x,y) for x,y in S)]
        print(i,",",SS,",")
