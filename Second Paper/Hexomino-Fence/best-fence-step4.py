"""
STEP 4

NO normalization of sum vector into the first octant.

Normalization of the solution:
EQUAL vectors are inserted into octants in weakly increasing order.
Every solution is generated once. (Or up to 8 times, depending
on the initial symmetry of the corner vectors.)

GENERATE all solutions.

IMPROVEMENT AND SPEED-UP AND MEMORY-SAVING IDEAS.
a) get rid of the "count" data in the first phase.

b) On the other hand, one can speed up the second phase by using the "last_octant" key
If last_octant has not been clipped to 0 because a different piece start, it tells
how the last previous piece was inserted, and this is enough to reconstruct the
previous partial solution and to figure out where to look in the old table.
(instead of scanning the whole old table.)
(This method would be even more powerful if the old table were kept as a dictionary.)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dynamic programming from some given list of vectors plus penalty

Pruning by area that can be reached with the current four sides (stored as "key")
and the remaining vectors, forming the largest polygon with these edge lengths.
That is: Accept only if the value "val" is decent enough.

No artificial or heuristic pruning of the search space. Yields the optimal
solution with the given reference vectors, provided "target_area"
is not too high.

For improvement of runtime (fewer calls of pruning function):
We check the pruning criterium only AFTER the table is formed,
and throw out the entries that fail.
(Some memory could be saved by checking the pruning criterium with every insertion.)

Possible variation or addition: use the Zauberformel (Minkowski-mean)
This has been tried, but it yields very poor bounds.

A weak form of meet-in-the middle:
Compute reachable vector sums from the end (a quite moderate number at most W**2)
Reduction of reachable set to first octant is not done, because the set is small anyway.

Here, we DO NOT rotate vector sum in solution into first
octant to normalize, because this messes up the counts.

If --storage-fileprefix is given, throw away the tables after every step, and
save it to disk, in order to save main memory.

Another idea: Push the reference vectors even more INSIDE the area.
Have to think carefully whether this is valid. (some arguments are in the paper).
This may kill some more options when the corner pieces are not exactly diagonal.
Empirically, it has no influence on the bound when all four corner vectors are
equal (and equal to (2.5,2.5).)

  3. Symmetry breaking:
 A) If the start configuration is fully symmetric (i.e. four equal diagonal
    corner vectors) and the steepest piece (the first one processed) is unique and
    lies neither on an axis nor on the diagonal, that piece is forced into octant 7
    (without loss of generality), which automatically excludes all other placements of that piece.
    Then every oriented polygon is generated exactly once. This is the case for n=6.
 B) If the start configuration is "almost" fully symmetric (i.e. three equal
    diagonal corner vectors and on different vector), then
    symmetries are removed by the (arbitrary) choice of the starting configuration already.
    (This happens for n=5).
 C) Otherwise the program aborts if start configurations can be formed that have any nontrivial symmetries.
    (Polygons might appear several times, and would have to be brought into canonical form.)

"""

import itertools
from collections import defaultdict, Counter
from math import sqrt,atan2
import math
#from marshal import dump,load
#from pickle import dump,load
# import bz2
import gzip # compression
import argparse
import importlib
import sys
from ast import literal_eval as evaluate # "somewhat safer" than eval()

from watch_memory import memory_checkpoint


_parser = argparse.ArgumentParser(
    description="Step 4: count and generate all optimal reference polygons.")
_parser.add_argument("n", type=int, help="size of the polyominoes (6 for hexominoes)")
_parser.add_argument("goal", type=int,
                     help="target area (0.001 is subtracted internally); must match the Step 1/2/3 runs")
_parser.add_argument("--only", type=int, metavar="SET_NUMBER", default=None,
                     help="process only this vector set (1-based number in the Step 3 file)")
_parser.add_argument("--max-solutions", type=int, metavar="N", default=10000,
                     help="stop after generating N different solutions (default: 10000)")
_parser.add_argument("--storage-fileprefix", type=str, metavar="DIR/xxx",
                     help='directory + file prefix for temporary files for intermediate store of the DP table (default: table is kept in memory)')
_parser.add_argument("--compress", action='store_true',
                     help='use compression for the temporary files (default: off)')
args = _parser.parse_args()
n_max, goal = args.n, args.goal
target_area = goal - 0.001


EXTERNAL_MEMORY = args.storage_fileprefix is not None
LOCAL_FILE_PREFIX = "/srv/data/rote/a"
LOCAL_FILE_PREFIX = "localdir/b"
LOCAL_FILE_PREFIX = args.storage_fileprefix

if args.compress:
    def open_dump_write(fname):
        return gzip.open(fname+".gz", "wt", encoding="utf-8", compresslevel=6)
    def open_dump_read(fname):
        return gzip.open(fname+".gz", "rt", encoding="utf-8")
else:
    def open_dump_write(fname):
        return open(fname+".txt", 'w')
    def open_dump_read(fname):
        return open(fname+".txt", 'r')
    
#target_area = 1596.999

def step_module_name(step, n_max, goal):
    return f"best_fence_candidate_solutions_step_{step}_n{n_max}t{goal}"
# Y for termporary hand-crafted files

def load_step_result(step, n_max, goal):
    """import the result file of Step 2 or Step 3; returns the module (its list is RESULT)"""
    name = step_module_name(step, n_max, goal)
    try:
        return importlib.import_module(name)
    except ModuleNotFoundError as e:
        if e.name != name:
            raise
        sys.exit(f'Input file "{name}.py" not found. Run Steps 1 to {step} first with the same '
                 f'arguments (n={n_max}, goal={goal}) and run this program from the directory '
                 f'that contains the result files and the scripts.')

CANDIDATES = load_step_result(3, n_max, goal).RESULT
SELECTED = CANDIDATES
if args.only is not None:
    SELECTED = [c for c in CANDIDATES if c["set_number"]==args.only]
    if not SELECTED:
        _parser.error(f"--only: vector set #{args.only} is not feasible. Feasible vector sets "
           f"according to Step 3: {[c['set_number'] for c in CANDIDATES]}")

print(f"Canonical-order DAG, {target_area=}")
print(f"Reconstruct {len(SELECTED)} of {len(CANDIDATES)} feasible vector sets: "
      f"{[c['set_number'] for c in SELECTED]}")


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


def best_polygon(Vecs,Corners):
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
    print(len(Start_configurations), "start configuration(s) of corners:")
    
    # symmetry breaking (see module docstring, point 3)
    fully_symmetric = (len(Start_configurations) == 1 and
                       all(len(Q) == 1 for Q in Start_configurations))
    if not fully_symmetric and any(len(Q)!=8 for Q in Start_configurations):
        print("There are start configurations with symmetries.")
        print("Unique generation not guaranteed. Aborting.")
        # this would require canonization in the end. Not implemented
        exit(1)
    if fully_symmetric:
        if len(Vecs)>=1:
            xs, ys = Vecs[-1]
            assert xs != ys and ys != 0
        if len(Vecs)>=2:
            assert Vecs[-1] != Vecs[-2]
        # otherwise symmetry is not broken
        force_last_octant = True
        l_start = 7
    else:
        force_last_octant = False
        l_start = 0
    print(f"fully symmetric start: {fully_symmetric};  "
          f"steepest piece {Vecs[-1] if Vecs else None} forced into octant 7: {force_last_octant}")
    
    
    TABLES = dict() # save the dictionaries to allow backtracing the solutions
    start = dict()
    for Q in Start_configurations:
        Q0 = list(Q)[0]
        print(*Q0)
        if Q0 == {(0,0)}:
            start[((0,0),(0,0),(0,0),(0,0)), l_start] = (0, 1)
        else: # split into the four quadrants:
            [(x1,y1)] = [(x,y) for x,y in Q0 if x>0 and y>0]
            [(x2,y2)] = [(x,y) for x,y in Q0 if x<0 and y>0]
            [(x3,y3)] = [(x,y) for x,y in Q0 if x<0 and y<0]
            [(x4,y4)] = [(x,y) for x,y in Q0 if x>0 and y<0]
            start[((x1,y1),(x2,y2),(x3,y3),(x4,y4)), l_start] = (0, 1)
        # The objective value of the initial solution, consisting only of four
        # single vectors, is 0, and is counted once (1).
    Table = TABLES[num] = [(key,l,vol,cnt) for (key,l),(vol,cnt) in start.items()]
    #Table = TABLES[num] = [a+b for a,b in Table.items()] # convert to list
    print("start table:", Table)
    angle_area_store = dict()
    def insert(octant,key,add_volume):
        # nonlocal vol,count,UU
        new_value = vol + add_volume
        if current_piece<=NUM_FROM_END:
            if (sum(x for x,y in key), sum(y for x,y in key)) not in reachable:
                return # throw away entries that cannot be completed
        l = min(octant,UU)
        try:
            old_val,old_count = newTable[key,l]
            if new_value>old_val:
                newTable[key,l] = new_value,count
            elif new_value==old_val:
                newTable[key,l] = new_value,count+old_count
        except KeyError:
            newTable[key,l] = new_value,count

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
    for current_piece in reversed(range(num)):
        (x0,y0) = Vecs[current_piece]
        if current_piece==0 or (x0,y0) != Vecs[current_piece-1]:
            UU = 0 # don't care in which octant the next piece lies
        else:
            UU = 8
        fixed_vectors = [sqrt(x**2+y**2) for x,y in Vecs[:current_piece]] # lengths
        sum_fixed_vectors = sum(fixed_vectors)
        angle_area_store = dict() # clear the cache
        if current_piece==NUM_FROM_END:
            print("reachability information kicks in.")
            # from now on: check reachability from the end in "insert".
        if current_piece<=NUM_FROM_END:
            reachable = reachablevectors[current_piece]
        newTable = dict()
        for ((x1,y1),(x2,y2),(x3,y3),(x4,y4)),last_octant, vol,count in Table:
            # order of product is ccw from positive x-axis
            # first octant, bundle in the direction to the right:
            if last_octant<=0:
                (x,y) = (x0,y0)
                insert(0,((x1+x,y1+y),(x2,y2),(x3,y3),(x4,y4)), det((x,y),(x1,y1)))
            # second octant, bundle in the upward direction:
            if last_octant<=1 and y0!=x0 and y0!=0: 
                # if y0==x0, the vector could eqally go to the second or first octant.
                # avoid unnecessary blow-up in this case by assigning to the first octant only.
                # if y0==0, the vector could eqally go to the second or third octant.
                # avoid unnecessary blow-up in this case by assigning to the third octant only.
                (x,y) = (y0,x0)
                insert(1,((x1+x,y1+y),(x2,y2),(x3,y3),(x4,y4)), det((x1,y1),(x,y)))
            # second quadrant
            if last_octant<=2:
                (x,y) = (-y0,x0) # bundle in the upward direction
                insert(2, ((x1,y1),(x2+x,y2+y),(x3,y3),(x4,y4)) ,det((x,y),(x2,y2)))
            if last_octant<=3 and y0!=x0 and y0!=0: 
                (x,y) = (-x0,y0) # bundle in the direction to the left
                insert(3, ((x1,y1),(x2+x,y2+y),(x3,y3),(x4,y4)), det((x2,y2),(x,y)))
            # third quadrant
            if last_octant<=4:
                (x,y) = (-x0,-y0)
                insert(4, ((x1,y1),(x2,y2),(x3+x,y3+y),(x4,y4)), det((x,y),(x3,y3)))
            if last_octant<=5 and y0!=x0 and y0!=0: 
                (x,y) = (-y0,-x0)
                insert(5, ((x1,y1),(x2,y2),(x3+x,y3+y),(x4,y4)), det((x3,y3),(x,y)))
            # fourth quadrant
            if last_octant<=6:
                (x,y) = (y0,-x0)
                insert(6, ((x1,y1),(x2,y2),(x3,y3),(x4+x,y4+y)),det((x,y),(x4,y4)))
            # last octant
            if y0!=x0 and y0!=0: 
                (x,y) = (x0,-y0)
                insert(7, ((x1,y1),(x2,y2),(x3,y3),(x4+x,y4+y)), det((x4,y4),(x,y)))
        # do the pruning only AFTER all insertions:
        beforepruning = len(newTable)
        print(memory_checkpoint(), f"iteration {current_piece}, piece={(x0,y0)}, len(Table) before pruning: {beforepruning:,}", flush=True)

        
        # newTable = [(key,l,val,c) for (key,l),(val,c) in newTable.items()
        #             if not prune_solution(fixed_vectors,[sqrt(x**2+y**2) for x,y in key], val/2)]

        # storage is at premium:
        # destroy the old dict as the new one is built. (But the hash table is not freed until the end.)
        newTable_list = []
        while newTable:
            (key,l),(val,c) = newTable.popitem()
            if not prune_solution(fixed_vectors,[sqrt(x**2+y**2) for x,y in key], val/2):
                newTable_list.append((key,l,val,c))
        newTable = newTable_list
            
        
        TABLES[current_piece] = Table = newTable
        if EXTERNAL_MEMORY:
            with open_dump_write(LOCAL_FILE_PREFIX+f"fence{n_max}t{goal}-{current_piece+1:02}") as dumpfile:
                for ((x1,y1),(x2,y2),(x3,y3),(x4,y4)),last_oct,vol,count in TABLES[current_piece+1]:
                    print(",".join(str(x) for x in (x1,y1,-x2,y2,-x3,-y3,x4,-y4,last_oct,vol)), file=dumpfile)
                # the signs of xi,yi are not stored, because they are known
                # count is not used in the restoration phase.
            del TABLES[current_piece+1] # release space
        print(memory_checkpoint(), f"iteration {current_piece}, piece={(x0,y0)}, len(Table) after pruning: {len(Table):,}", flush=True)
    best = max((vol + det((x1,y1),(x2,y2)) + det((x3,y3),(x4,y4))
               # det(v1+v2,v3+v4)==0 because the two vectors are opposite.
               # Therefore only det(v1,v2)+det(v3,v4) remains to be added
               for ((x1,y1),(x2,y2),(x3,y3),(x4,y4)),_, vol,_ in Table
               if (x1+x4+x3+x2, y1+y4+y3+y2)==(0,0)), default=None)
    if best is None:
        print(f"No solution with area at least {target_area} was found.")
        return None,None
    
    print(best,"best value", flush=True)
    numbest = 0
    oldSol = defaultdict(list)
    for ((x1,y1),(x2,y2),(x3,y3),(x4,y4)),l,vol,count in Table:
        if (x1+x4+x3+x2, y1+y4+y3+y2)==(0,0) and vol + det((x1,y1),(x2,y2)) + det((x3,y3),(x4,y4))==best:
            numbest += count
            oldSol[((x1,y1),(x2,y2),(x3,y3),(x4,y4)),l,vol].append([])
    print(f"Number of optimal solutions: {numbest}")
    
    def find(octant,key,add_volume):
        # nonlocal vol,count, x,y, (x1,y1)..., last_octant,
        new_value = vol + add_volume
        l = min(octant,UU)
        for partial_sol in oldSol.get((key,l,new_value),[]):
            newSol[((x1,y1),(x2,y2),(x3,y3),(x4,y4)),last_octant,vol].append(partial_sol+[(x,y)])
        
    # DO THE MAIN LOOP IN REVERSE

    def read_from_dumpfile(line):
        x1,y1,x2,y2,x3,y3,x4,y4,last_oct,vol = evaluate(line)
        count = 0 # dummy entry, count is not used
        x1,y1,x2,y2,x3,y3,x4,y4 = x1,y1,-x2,y2,-x3,-y3,x4,-y4
        return ((x1,y1),(x2,y2),(x3,y3),(x4,y4)),last_oct,vol,count
    
    for current_piece in (range(num)):
        (x0,y0) = Vecs[current_piece]
        if current_piece==0 or (x0,y0) != Vecs[current_piece-1]:
            UU = 0 # don't care in which octant the next piece lies
        else:
            UU = 8
        if EXTERNAL_MEMORY:
            dumpfile = open_dump_read(LOCAL_FILE_PREFIX+f"fence{n_max}t{goal}-{current_piece+1:02}")
            def get_next_line():
                for line in dumpfile: # read directly from memory without storing
                    yield read_from_dumpfile(line)
            Table = get_next_line()
        else:
            Table = TABLES[current_piece+1]
        newSol = defaultdict(list)
        # exactly as before, except that "insert" is replaced by "find"
        for ((x1,y1),(x2,y2),(x3,y3),(x4,y4)),last_octant, vol,count in Table:
            # order of product is ccw from positive x-axis
            # first octant, bundle in the direction to the right:
            if last_octant<=0:
                (x,y) = (x0,y0)
                find(0,((x1+x,y1+y),(x2,y2),(x3,y3),(x4,y4)), det((x,y),(x1,y1)))
            # second octant, bundle in the upward direction:
            if last_octant<=1 and y0!=x0 and y0!=0: 
                # if y0==x0, the vector could eqally go to the second or first octant.
                # avoid unnecessary blow-up in this case by assigning to the first octant only.
                # if y0==0, the vector could eqally go to the second or third octant.
                # avoid unnecessary blow-up in this case by assigning to the third octant only.
                (x,y) = (y0,x0)
                find(1,((x1+x,y1+y),(x2,y2),(x3,y3),(x4,y4)), det((x1,y1),(x,y)))
            # second quadrant
            if last_octant<=2:
                (x,y) = (-y0,x0) # bundle in the upward direction
                find(2, ((x1,y1),(x2+x,y2+y),(x3,y3),(x4,y4)) ,det((x,y),(x2,y2)))
            if last_octant<=3 and y0!=x0 and y0!=0: 
                (x,y) = (-x0,y0) # bundle in the direction to the left
                find(3, ((x1,y1),(x2+x,y2+y),(x3,y3),(x4,y4)), det((x2,y2),(x,y)))
            # third quadrant
            if last_octant<=4:
                (x,y) = (-x0,-y0)
                find(4, ((x1,y1),(x2,y2),(x3+x,y3+y),(x4,y4)), det((x,y),(x3,y3)))
            if last_octant<=5 and y0!=x0 and y0!=0: 
                (x,y) = (-y0,-x0)
                find(5, ((x1,y1),(x2,y2),(x3+x,y3+y),(x4,y4)), det((x3,y3),(x,y)))
            # fourth quadrant
            if last_octant<=6:
                (x,y) = (y0,-x0)
                find(6, ((x1,y1),(x2,y2),(x3,y3),(x4+x,y4+y)),det((x,y),(x4,y4)))
            # last octant
            if y0!=x0 and y0!=0: 
                (x,y) = (x0,-y0)
                find(7, ((x1,y1),(x2,y2),(x3,y3),(x4+x,y4+y)), det((x4,y4),(x,y)))
        oldSol = newSol
        if EXTERNAL_MEMORY:
            dumpfile.close()
        print(memory_checkpoint(), f"reverse iteration {current_piece}, piece={(x0,y0)}, {len(newSol)=:,},",
              f"partial solutions {sum(len(x) for x in oldSol.values()):,}", flush=True)
        
    SOLS = []
    for (corn,_,_),sols in oldSol.items():
        corn2 = []
        for x,y in corn:                    
            x -= shift if x>0 else -shift
            y -= shift if y>0 else -shift
            corn2.append((x,y))
        jj = 0
        for part in sols:
            SOL = part+corn2
            print(jj, SOL)
            jj +=1
            S = min([(x,y) for _,x,y in sorted(((atan2(y,-x),x,y)
                                               for x,y in r1(r2(SOL))),
                                               reverse=True)]
                    for r1,r2 in ALL_SYMMETRIES)
            SOLS.append(tuple(S))
    return best,SOLS
    
print(f"Finding best fences. {target_area=}")
print(f"Constructing all solutions.")
print(f"Solve {len(SELECTED)} vector sets.")
print(f"Using {EXTERNAL_MEMORY=}" +
      (f", saving to files starting with {LOCAL_FILE_PREFIX}, compression={args.compress}"
       if EXTERNAL_MEMORY else ""))

FIND = 0
ALL_SOLS = []
for nline,vec_family in enumerate(SELECTED,start=1):
    vecs = vec_family['vectors']
    penalty = -vec_family['neg_penalty']
    assert penalty>0
    print("\n======================")
    print(f"Vector set #{nline} of {len(SELECTED)}, index {vec_family['set_number']}:",vecs)
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
        assert x>=y
        penalty += x*shift # x is the larger coordinate
        
        
    Corners = Corners_new
    print(f"adusted {penalty=}, after moving the reference vectors.")
    Vecs = [p for (s,p) in sorted(Vecs)] # steepest last
    assert all(c==int(c) for p in Vecs+Corners for c in p)
    Vecs    = [(int(x),int(y)) for x,y in Vecs]
    Corners = [(int(x),int(y)) for x,y in Corners]
    
    total_parity_sum = sum(x+y for (x,y) in Vecs+Corners)
    print("total L1-distance must be even:", total_parity_sum)
    assert total_parity_sum%2 == 0 # otherwise solution is impossible.

    area_threshold = target_area + penalty

    bin_search_calls = bin_search_decisions = 0
    ### call the main routine:
    B, SOLS = best_polygon(Vecs,Corners)
    print(f"Total number of binary search calls = {bin_search_calls:,},",
          f"average number of decisions per call =",
          f"{bin_search_decisions/bin_search_calls:.5}")

    if B:
        print(f"best solution value {B/2}. **** enclosed area={B/2-penalty=} ************")
        enclosedarea = B/2-penalty
        if enclosedarea<=target_area:
            print(f"This is not guaranteed to be the optimum, because < {target_area=}")
        else:
            FIND += 1
            print(f"This is guaranteed to be the optimum.")
        print(flush=True)
        ALL_SOLS += [(sol,vec_family['neg_penalty']) for sol in SOLS]

        num = len(Vecs)
with open(f"all-reference-polygons-{n_max}t{goal}.txt", 'w') as resultfile:       
    if ALL_SOLS:
        print("All solutions together in sorted order:") # print the solutions in the end
        for i,(sol,correction) in enumerate(sorted(ALL_SOLS),start=1):
            print(f"{i:3}, {sol}, {correction=}")
            print(f"{(i, correction, sol)}", file = resultfile)
        assert len(set(ALL_SOLS))==len(ALL_SOLS)
print(f"Among {len(SELECTED)} vector sets, a solution larger than {target_area} was found in {FIND} instance(s).")
print("total runtime:", memory_checkpoint())
