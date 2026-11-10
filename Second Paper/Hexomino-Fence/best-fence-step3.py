"""
STEP 3. Dynamic programming with fixed vectors.
Take the vector sets that survived Step 2
(File best_fence_candidate_solutions_step_2_n<n>t<goal>.py)
and run the dynamic program on each of them.
The result is the list of FEASIBLE vector sets, i.e. those for which a fence can actually
be constructed, together with the best enclosed area that is attainable with them.
It is written to the file best_fence_candidate_solutions_step_3_n<n>t<goal>.py,
which is the input of Step 4 (reconstruction of the solutions).
=================================
Usage: pypy3 best-fence-step3.py <n> <goal>

  n       size of the polyominoes (6 for hexominoes)
  goal        target area; the program subtracts 0.001 to guard against rounding errors.
              Must be the same as in the Step 1 and Step 2 runs that produced the input file.
              If a solution with at least this value exists, the program will find it.
              Otherwise the program may find no solution or a solution that is not optimal.

Only the optimal value of every vector set is computed, not the solution itself. That is
sufficient to show that no solution larger than the target area exists, and it saves a lot
of space, because the whole sequence of tables need not be remembered.

Output file: a Python list RESULT with one dictionary per feasible vector set
(vector sets without any construction are not listed):
  set_number     1-based number of the vector set in the Step 2 file
                 (the "Vector set #k" of the output and of the README)
  vectors        the vector set, as in the Step 2 file
  neg_penalty    the correction value (negative penalty), as in the Step 2 file
  enclosed_area  best enclosed area attainable with this vector set, if at least the goal area
"""
import os
import argparse
import importlib
import sys
import itertools
from collections import defaultdict
from math import sqrt,atan2
import math
from types import SimpleNamespace
from watch_memory import memory_checkpoint

_parser = argparse.ArgumentParser(
    description="Step 3: dynamic programming on the vector sets that survived Step 2.")
_parser.add_argument("n", type=int, help="size of the polyominoes (6 for hexominoes)")
_parser.add_argument("goal", type=int,
                     help="target area (0.001 is subtracted internally); must match the Step 1/2 runs")
args = _parser.parse_args()
n_max, goal = args.n, args.goal
target_area = goal - 0.001
 # If a solution with at least this value exists, the program will find it.
 # Otherwise the program may find no solution or a solution that is not optimal.

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
# octant to normalize (saves a factor 8 in most cases).
# (Step 4 does not normalize, because it counts and reconstructs all solutions.)

# Only the current table is kept in memory; the earlier tables are discarded.
# The solutions themselves are reconstructed in Step 4.

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
total_bin_search_calls = total_bin_search_decisions = use_semicircle = 0

def step_module_name(step, n_max, goal):
    return f"best_fence_candidate_solutions_step_{step}_n{n_max}t{goal}"

def load_step_result(step, n_max, goal):
    """import the result file of Step 2 or Step 3; returns the module (its list is .RESULT)"""
    name = step_module_name(step, n_max, goal)
    print (f"importing from {name}.py")
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
    # global area_threshold, R_min, sum_fixed_vectors
    global bin_search_calls, bin_search_decisions, use_semicircle
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
            if R_2**2*math.pi/2 + additional_area <=area_threshold:
                # use area of semicircle as a (loose upper bound)
                return True # prune
            else:
                use_semicircle += 1
                return False # remain on the safe side: don't prune
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

def forward_dp(Vecs,Corners):
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

    Table = list(Table.items()) # from now on, work with a list
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
        for (((x1,y1),(x2,y2),(x3,y3),(x4,y4)),vol) in Table:
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
        print(memory_checkpoint(), f"iteration {i=} piece={(x0,y0)}, len(Table) before pruning: {beforepruning:,}", flush=True)

        # pruning with Minkowski-mean (omitted)

        newTable_list = []
        while newTable:
            (key,val) = newTable.popitem()
            if not prune_solution(fixed_vectors,[sqrt(x**2+y**2) for x,y in key], val/2):
                newTable_list.append((key,val))
        TABLES [i] = Table = newTable = newTable_list
        del TABLES[i+1] # release space
        print(memory_checkpoint(), f"iteration {i=}, piece={(x0,y0)}, len(Table) after pruning: {len(Table):,}", flush=True)
    return TABLES # only TABLES[0] contains something


def run_forward_dp(Vecs, Corners, threshold):
    """Forward pass of the dynamic program.
    threshold: the value area_threshold = target_area + penalty used for pruning.
    Returns best,finals; best is None if there is no solution at all.
    finals = list of (key, value) of the optimal final states."""
    # Only the optimal value is returned (saves space, the tables of the previous stages are thrown away)
    global area_threshold, bin_search_calls, bin_search_decisions, total_bin_search_calls, total_bin_search_decisions
    area_threshold = threshold
    bin_search_calls = bin_search_decisions = 0
    TABLES = forward_dp(Vecs, Corners)
    print(f"Total number of binary search calls = {bin_search_calls:,},",
          f"average number of decisions per call =",
          f"{bin_search_decisions/max(bin_search_calls,1):.5}")
    total_bin_search_decisions += bin_search_decisions
    total_bin_search_calls += bin_search_calls
    
    Table = TABLES[0]
    best = max((vol + det((x1,y1),(x2,y2)) + det((x3,y3),(x4,y4))
               # det(v1+v2,v3+v4)==0 because the two vectors are opposite.
               # Therefore only det(v1,v2)+det(v3,v4) remains to be added
               for (((x1,y1),(x2,y2),(x3,y3),(x4,y4)), vol) in Table
               if (x1+x4+x3+x2, y1+y4+y3+y2)==(0,0)), default=None)
    if best is None:
        return None, [] # (the caller reports this)
    numbest = 0
    finals = []
    for ((x1,y1),(x2,y2),(x3,y3),(x4,y4)),vol in Table:
        if (x1+x4+x3+x2, y1+y4+y3+y2)==(0,0) and vol + det((x1,y1),(x2,y2)) + det((x3,y3),(x4,y4))==best:
            bestsol = (x1,y1),(x2,y2),(x3,y3),(x4,y4)
            numbest += 1
            if numbest<=30: # the best 30 are printed
                print(f"A best solution ({numbest}): {bestsol}")
            finals.append((bestsol,vol))
    if numbest>30:
        print("...")
    print(f"{numbest}  best solutions in total")
    print(best,"best value", flush=True)
    return best, finals

def prepare_instance(nline, vecs, neg_penalty):
    """Turn one entry (vecs, neg_penalty) of the Step 2 file into the input of the
    dynamic program. Returns Vecs, Corners, shift, penalty."""
    penalty = -neg_penalty
    assert penalty>0
    print("\n======================")
    print(f"Vector set #{nline} of {len(INSTANCES)}:",vecs)
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
    if len(Corners)!=4:
        print("The number of corner candidates is not 4.")
        print("The program is not prepared for this case. Exit.")
        exit()

    # corners need not be sorted by slope.
    #for shift in -0.5,+0.5,+1.5,+2.5,3.5,4.5 ...:
    for ishift in range(20):
        shift = ishift-0.5
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
        exit(1)

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


def solve_vector_set(nline, vecs, neg_penalty, target_area):
    """Prepare one vector set and run the forward pass. Returns a namespace with
    Vecs, Corners, shift, penalty, best, enclosed_area (None if there is no construction at all)."""
    Vecs, Corners, shift, penalty = prepare_instance(nline, vecs, neg_penalty)
    best, finals = run_forward_dp(Vecs, Corners, target_area + penalty)
    enclosed_area = None
    if best is None:
        print(f"No solution with area at least {target_area} was found.")
    else:
        B = best
        print(f"best solution value {B/2}. **** enclosed area={B/2-penalty=} ************")
        enclosed_area = B/2-penalty
        if enclosed_area<=target_area:
            print(f"This is not guaranteed to be the optimum, because it is < {target_area=}",flush=True)
        else:
            print(f"This is guaranteed to be the optimum.", flush=True)
    return SimpleNamespace(Vecs=Vecs, Corners=Corners, shift=shift, penalty=penalty,
                           best=best, enclosed_area=enclosed_area)


# ---------------------------------------------------------------------------
       
INSTANCES = load_step_result(2, n_max, goal).RESULT

print(f"Only the optimal value is computed, {target_area=}")
print(f"Solve {len(INSTANCES)} vector sets.")
FIND = 0
CANDIDATES = []
INFEASIBLE = []
for nline,(vecs,neg_penalty) in enumerate(INSTANCES,start=1):
    r = solve_vector_set(nline, vecs, neg_penalty, target_area)
    if r.enclosed_area is not None and r.enclosed_area>=target_area: # has reached the goal
        FIND += 1
        CANDIDATES.append({"set_number": nline,
                           "vectors": vecs,
                           "neg_penalty": neg_penalty,
                           "enclosed_area": r.enclosed_area})
    else:
        INFEASIBLE.append(nline)

# write the result file into the current directory (where Step 2 left its file)
outname = step_module_name(3, n_max, goal) + ".py"

with open(outname, "w") as f:
    f.write(f"# generated by the program best-fence-step3.py with target={target_area}\n")
    f.write(f"# Step 2 file: {step_module_name(2, n_max, goal)}.py ({len(INSTANCES)} vector sets)\n")
    f.write(f"# feasible vector sets (a fence can be constructed): {[c['set_number'] for c in CANDIDATES]}\n")
    f.write("RESULT=[\n")
    for c in CANDIDATES:
        f.write(repr(c) + " ,\n")
    f.write("]\n")

print()
print(f"Feasible vector sets: {[c['set_number'] for c in CANDIDATES]}")
print(f"Without any solution at least goal: {INFEASIBLE}")
print(f"Vector sets reaching the goal: {[c['set_number'] for c in CANDIDATES]}")
print(f"Among {len(INSTANCES)} vector sets, a solution larger than {target_area} was found in {FIND} instance(s).")
print(f"Wrote {len(CANDIDATES)} feasible vector set(s) to `{outname}`.")
print(f"Overall Total number of binary search calls = {total_bin_search_calls:,},",
      f"average number of decisions per call =",
      f"{total_bin_search_decisions/max(total_bin_search_calls,1):.5}")
print("# cases that were not pruned when the center of the circumscribed circle does not lie in the polygon:",
      use_semicircle)


print("runtime and memory",memory_checkpoint())
