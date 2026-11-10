# a weak form of meet-in-the middle:
# compute reachable vector sums from the end (a quite moderate number at most W**2

# this version: table is already pruned when inserting.
# reduction of reachable set to first octant is undone (because the set is small).

# In addition: rotate sum in solution into first
# octant to normalize (saves a factor 8 in most cases)
# Der Rekonstruktionsvorgang wurde entsprechend angepasst.
# Das geht aber nur, wenn die 4 Eckvektoren gleich sind? NÖ! (Warum eigentlich?)

# In addition: accept only if the four vectors
# (x1,y1),(x2,y2),(x3,y3),(x4,y4) are roughly the same (after normalization)
# tested with threshold CUTOFF3

# Another idea. accept only if the value "val" is decent enough.
# (to be figured out what that means.)

# Another idea. Throw away the tables more than, say, 10 steps back.
# saves space. Aber dann muss man sich die Lösung in mehreren Läufen
# rekonstruieren, und Schritt für Schritt zurückhangeln.


"""HACK of the pieces:
variation with the "T" piece (code "S") instead of the "C" piece (code "K")
"""

import itertools
from collections import defaultdict
from math import sqrt
import math

# Parameters:

# CUTOFF bounds the sum of the edge-vectors componentwise (L_max distance)
# idea: sum should stay close to the origin, alread for partial sums.

CUTOFF = 20*2 # pruning of the search space.
#CUTOFF = 8*2 # pruning of the search space.

#CUTOFF2 = 10 # not used

# CUTOFF3: pairwise L1 distance between the four vector sums.
CUTOFF3 = 35

from sys import stdout
import sys

if len(sys.argv)>1:
    CUTOFF = int(sys.argv[1])
if len(sys.argv)>2:
    CUTOFF3 = int(sys.argv[2])

print(f"{CUTOFF=}")
print(f"{CUTOFF3=}")

Polyominoes = [[],[((0,0),)]]

n_max = 8
n_max = 7
n_max = 6

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
def apply_mat(A,x):
    a11,a12,a21,a22 = A
    x1,x2 = x
    return (a11*x1+a12*x2, a21*x1+a22*x2)

def normalize(P): # lowest coordinates become 0
    xmin = min(x for x,y in P)
    ymin = min(y for x,y in P)
    return [(x-xmin,y-ymin) for x,y in P]
def print_polyo(P,letter="X"):
    xmax = max(x for x,y in P)
    ymax = max(y for x,y in P)
    pattern = [[" " for _ in range(1+xmax)] for _ in range(1+ymax)]
    for x,y in P:
        pattern[ymax-y][x]=letter
    print("\n".join("".join(l) for l in pattern))

def print_grid(GRID, text=""):    
    # print the pattern represented in the dictionary GRID
    xmin = min(x for x,y in GRID.keys())
    xmax = max(x for x,y in GRID.keys())
    ymin = min(y for x,y in GRID.keys())
    ymax = max(y for x,y in GRID.keys())
    pattern = [["." for _ in range(1+xmax-xmin)] for _ in range(1+ymax-ymin)]
    for (x,y),letter in GRID.items():
        pattern[ymax-y][x-xmin]=letter
    print("\n".join("".join(l) for l in pattern) + text + f",  {xmin}<=x<={xmax}, {ymin}<=y<={ymax}")

    return (xmin+xmax)//2, (ymin+ymax)//2, #  midpoint, hopefully an interior point

    # ipe-file for graphics, optional
    with open(f"out{n}.ipe","w") as out:
        print(f"""<group matrix="10 0 0 10 100 100">""", file=out)
        for y in range(ymin-1,ymax+3):           
            print(f"""<path stroke="black">
               {xmin-1.5} {y-0.5} m
               {xmax+1.5} {y-0.5} l
               </path>""", file=out)
        for x in range(xmin-1,xmax+3):           
            print(f"""<path stroke="black">
               {x-0.5} {ymin-1.5} m
               {x-0.5} {ymax+1.5} l
               </path>""", file=out)
        for (x,y),letter in GRID.items():
            num = sequencenumber[letter]
            color = f"0.6 {num/40} {1-num/49}"
            if num%2:
                color = f"{num/40} {1-num/49} 0.8"
            print(f"""<path fill="{color}">
            {x-0.5} {y-0.5} m
            {x+0.5} {y-0.5} l
            {x+0.5} {y+0.5} l
            {x-0.5} {y+0.5} l
            h
            </path>""", file=out)
            for (dx,dy) in neighbors:
                xx,yy = x+dx,y+dy
                if GRID.get((xx,yy)," ") != letter:
                    perpx,perpy = -dy/2,dx/2
                    print(f"""<path stroke="black" pen="heavier">
                    {(x+xx)/2+perpx} {(y+yy)/2+perpy} m
                    {(x+xx)/2-perpx} {(y+yy)/2-perpy} l
                    </path>""", file=out)
        print(f"""</group>""", file=out)
    
    return (xmin+xmax)//2, (ymin+ymax)//2, #  midpoint, hopefully an interior point

# When comparing the options for choosing the reference vectors,
# the DIAMETER is measured for the actual vectors, not the
# distance between cell centers. The first case where this makes a difference
# is between (4,1) and (3,3), i.e. starting from octominoes.
# (In a heptomino, these difference vectors cannot arise in the same piece.)
#   a**2    + b**2   is 17 versus 18
#  (a+1)**2 + b**2   is 26 versus 25.

def diam_adjusted(a,b):
    """computes length of the resulting reference vector,
    assuming the piece is not used as a corner piece.
    the larger component is increased by one."""
    return a**2 + b**2 + 2*max(abs(a),abs(b)) + 1

def diameter(P):
    """returns the diameter and all diametrical pairs of cells"""
    d = max(diam_adjusted(x-x2,y-y2)
            for (x,y) in (P)
            for (x2,y2) in P)
    dpairs = [((x,y),(x2,y2))
            for i,(x,y) in enumerate(P)
            for j,(x2,y2) in enumerate(P)
              if i<j and diam_adjusted(x-x2,y-y2)==d]
    return d,len(dpairs),dpairs

neighbors4 = (-1,0),(1,0),(0,1),(0,-1) # 4-neighborhood
neighbors8 = (-1,0),(1,0),(0,1),(0,-1),(1,1),(1,-1),(-1,1),(-1,-1) # 8-neighborhood

def orient(a,b,c): # orientation predicate
    (x1,y1),(x2,y2),(x3,y3) = a,b,c
    return (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)

def det(a,b): # area of parallelogram spanned by two vectors
    (x1,y1),(x2,y2) = a,b
    return x1*y2-y1*x2

def bestposition(d,_,dpairs):
    # the best orientation among all diametrical pairs
    # taking into account the occupied and free cells.
    xmin = min(x for x,y in P)
    ymin = min(y for x,y in P)
    xmax = max(x for x,y in P)
    ymax = max(y for x,y in P)
    boundary = set([(x,y) for x in (-1,1+xmax) for y in range(-1,2+ymax)] +
                   [(x,y) for y in (-1,1+ymax) for x in range(-1,2+xmax)])
    best_part = None
    for dpair in dpairs:
        (x1,y1),(x2,y2) = dpair
        cell = {(i,j):" " for i in range(-1,2+xmax)
                   for j in range(-1,2+ymax)}
        for x,y in P:
            cell[x,y]="X"
        # mark the two boundary segments:
        for x,y in ((x1,y1),(x2,y2)):
            # adjacent to endpoints of diametrical pair becomes "Q"
            for (dx,dy) in neighbors4:
                xx,yy = x+dx,y+dy
                if (xx,yy) in boundary:
                    cell[xx,yy]="Q"
            # add corner cells if appropriate
            if (x,y)==(0,0):
                cell[-1,-1]="Q"
            if (x,y)==(xmax,0):
                cell[xmax+1,-1]="Q"
            if (x,y)==(0,ymax):
                cell[-1,ymax+1]="Q"
            if (x,y)==(xmax,ymax):
                cell[xmax+1,ymax+1]="Q"
        # look for the two components of the boundary, and extend from
        # there into the interior
        num_comp = 0
        for x0,y0 in boundary:
            if cell[x0,y0]==" ":
                num_comp += 1
                size = 0
                cell[x0,y0]= num_comp
                todo = [(x0,y0)]
                while todo:
                    (x,y) = todo.pop()
                    for (dx,dy) in neighbors8:
                        xx,yy = x+dx,y+dy
                        if (-1<=xx<=xmax+1 and -1<=yy<=ymax+1
                            and cell[xx,yy]==" "):
                            cell[xx,yy]= num_comp
                            todo.append((xx,yy))
                            if (xx,yy) not in boundary:
                                size += 1

                ori0 = orient((x1,y1),(x2,y2),(x0,y0))
                assert ori0 != 0
                # compute 2x the number of cells of the enclosing rectangle
                # on the "ori0" side of the diametrical chord:
                count = sum(2 if orient((x1,y1),(x2,y2),(x,y))*ori0>0
                            else 1 if orient((x1,y1),(x2,y2),(x,y))==0
                            else 0
                            for x in range(1+xmax) for y in range(1+ymax))
                size = size-count/2 # the net contribution
                # remember the best solution(s)
                if best_part is None or size>best_part:
                    best_part = size
                    best_sol = [(dpair,ori0)]
                elif size==best_part:
                    best_sol.append((dpair,ori0))
                #print(f"{num_comp=}, {size=}")
        assert n<=2 or num_comp == 2
        #print(f"{((x1,y1),(x2,y2))=}"),
        #for j in range(1+ymax,-2,-1):
        #    for i in range(-1,2+xmax):
        #        print(end=f"{cell[i,j]:1}")
        #    print()
    if best_part is None:
        return
    return best_part,best_sol

def count_area(GRID, mid):
    """determine size of the enclosed area.
    midpoint should supposedly be surrounded by the GRID cells.
    Cells occupied by the fence are identified by having their
    key in the dictionary "GRID"."""
    assert mid not in GRID
    found = set([mid])
    todo = [mid]
    size = 1
    while todo:
        (x,y) = todo.pop()
        for (dx,dy) in neighbors8:
            xx,yy = x+dx,y+dy
            if (xx,yy) not in GRID and (xx,yy) not in found:
                found.add((xx,yy))
                todo.append((xx,yy))
                size += 1
    return size

ALL_SYMMETRIES = list(itertools.product((rotate0,reflect),
                                        (rotate0,rotateL,rotate2,rotateR)))

def prune_keys(key):
    (x1,y1),(x2,y2),(x3,y3),(x4,y4) = key
    sx,sy = x1+x2+x3+x4, y1+y2+y3+y4
    if max(abs(sx), abs(sy))>CUTOFF:
        return True # throw away some entries to save space
    #xrange = max(x1,-x2,-x3,x4) - min(x1,-x2,-x3,x4)
    #if xrange>CUTOFF2:  return True
    #yrange = max(y1,y2,-y3,-y4) - min(y1,y2,-y3,-y4)
    #if yrange>CUTOFF2:  return True
    normpoints = (x1,y1),(-x2,y2),(-x3,-y3),(x4,-y4) # moved in the first quadrant
    maxL1 = (  max(x+y for (x,y) in normpoints) - min(x+y for (x,y) in normpoints)
             + max(x-y for (x,y) in normpoints) - min(x-y for (x,y) in normpoints))
    # maximum pairwise L1 distance
    # previously: maximum pairwise squared Euclidean distance.
    # For smaller cutoff, this the tables tended to dry out in the end when the longer pieces
    # are added. Example: octominoes, cutoff=150.
    if maxL1>CUTOFF3:  return True
    return False
        

def best_polygon(Vecs):
    """The workhorse of this algorithm.
    Vecs (x,y) satisfy 0<=y<=x and are sorted by slope."""
    Vecs = [(int(2*x),int(2*y)) for (x,y) in Vecs] # integers
    def insert(key,val):
        if prune_keys(key):
            return        
        if CORNERS_EQUAL:
            (x1,y1),(x2,y2),(x3,y3),(x4,y4) = key
            sx,sy = x1+x2+x3+x4, y1+y2+y3+y4
            # normalize the sum into the first octant
            # (remember: start table = {((7, 7), (-7, 7), (-7, -7), (7, -7)): 196})
            if sx<0:
                (x1,y1),(x2,y2),(x3,y3),(x4,y4) = key = (-x2,y2),(-x1,y1),(-x4,y4),(-x3,y3)
                sx = -sx
                ## ... adjust val as well????? !!!!
            if sy<0:
                (x1,y1),(x2,y2),(x3,y3),(x4,y4) = key = (x4,-y4),(x3,-y3),(x2,-y2),(x1,-y1)
                sy = -sy
            if sy>sx:
                # (x1,y1),(x2,y2),(x3,y3),(x4,y4) =
                key = (y1,x1),(y4,x4),(y3,x3),(y2,x2)
        if key not in newTable or val>newTable[key]:
            newTable[key] = val
    num = len(Vecs)

    # By convention, the vectors that are to be used as corners
    # are the last 4 vectors. (Those that are closest to the 45° direction.)
    corners = Vecs[-4:]
    assert all(s%2 for v in corners for s in v)
    expected_corner_length = n - 1 + (n%2)
    if n>=6:
        assert set(corners) == {(expected_corner_length,expected_corner_length)}
    else:
        assert set([v for v in corners[1:]]) == {(5,5)}
    scorners = list(set(corners))
    CORNERS_EQUAL = len(scorners)==1 and scorners[0][0]==scorners[0][1]
    CORNERS_EQUAL = True
        
    v1,v2,v3,v4 = corners
    (x1,y1)=v1
    (x2,y2)=v2; x2=-x2 # rotate the corner vectors into the 4 quadrants.
    (x3,y3)=v3; x3=-x3; y3=-y3
    (x4,y4)=v4; y4=-y4
    # vol = objective value, will accumulate to half the area
    vol = 0 #det(v1,(x2+x3+x4,y2+y3+y4)) + det((x2,y2), (x3+x4,y3+y4)) + det((x3,y3),(x4,y4))
    # objective value of the initial solution, consisting only of the four
    # corner vectors.

    # a weak form of meet-in-the middle:
    # compute reachable vectors from the end
    NUM_FROM_END = min((num-4)//2+6, num-4) # how many vectors are combined, at most
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
                                     if max(abs(x+u),abs(y+v))<=CUTOFF))
        print(f"len(reachablevectors[{i+1}])={len(reachablevectors[i+1])}")
    
    TABLES = dict() # save the dictionaries to allow backtracing the solutions
    Table = TABLES[num-4] = dict()
    #    Table[(Vecs[0],(0,0),(0,0),(0,0))] = 0
    CORNERS = [v1,(x2,y2),(x3,y3),(x4,y4)]
    #print(f"{CORNERS=}")
    Table[v1,(x2,y2),(x3,y3),(x4,y4)] = vol
    print("start table:", Table)
    for i in reversed(range(num-4)):
        (x0,y0) = Vecs[i]
        if i==NUM_FROM_END:
            print("reachability information kicks in.")
        if i<=NUM_FROM_END:
            # redefine "insert". from now on: check reachability from the end
            reachable = reachablevectors[i]
            def insert(key,val): 
                (x1,y1),(x2,y2),(x3,y3),(x4,y4) = key
                sx,sy = x1+x2+x3+x4, y1+y2+y3+y4
                if (sx, sy) not in reachable:
                    # CUTOFF already included in "reachable"
                    return # throw away some entries to save space
                if prune_keys(key):
                    return                       
                if CORNERS_EQUAL:
                    # normalize the sum into the first octant
                    # remember: start table = {((7, 7), (-7, 7), (-7, -7), (7, -7)): 196}
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
        TABLES [i] = Table = newTable
        print(f"iteration {i=} piece={(x0,y0)}, {len(Table)=}", flush=True)
#        print(Table)
    best = max(vol + det((x1,y1),(x2,y2)) + det((x3,y3),(x4,y4))
               # det(v1+v2,v3+v4)==0 because the two vectors are opposite.
               # Therefore only det(v1,v2)+det(v3,v4) remains to be added
               for (((x1,y1),(x2,y2),(x3,y3),(x4,y4)), vol) in Table.items()
               if (x1+x4+x3+x2, y1+y4+y3+y2)==(0,0))
    for (((x1,y1),(x2,y2),(x3,y3),(x4,y4)),vol) in Table.items():
        if (x1+x4+x3+x2, y1+y4+y3+y2)==(0,0) and vol + det((x1,y1),(x2,y2)) + det((x3,y3),(x4,y4))==best:
            bestsol = (x1,y1),(x2,y2),(x3,y3),(x4,y4)
            bestvol = vol
            print(f"A best solution: {bestsol}")

    print(best,"best value", flush=True)

    # Now the optimum value has been determined.
    # We identify the solution by tracking backward:
    global SOL # list of vectors used in the solution (in addition to the corner vectors)
    SOL = []
    MAT = IDENTITY
    def find(key,vol):
        nonlocal MAT # also: x,y, bestsol, bestvol
        if vol!=bestvol:
            return False
        if CORNERS_EQUAL:
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
        if CORNERS_EQUAL:
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
    for i in range(num-4):
        (x0,y0) = Vecs[i]
        for (((x1,y1),(x2,y2),(x3,y3),(x4,y4)),vol) in TABLES[i+1].items():
            # order of product is ccw from positive x-axis
            # first octant
            (x,y) = (x0,y0)
            if find(((x1+x,y1+y),(x2,y2),(x3,y3),(x4,y4)),vol+det((x,y),(x1,y1))): break
            # second octant
            (x,y) = (y0,x0)
            if find(((x1+x,y1+y),(x2,y2),(x3,y3),(x4,y4)),vol+det((x1,y1),(x,y))): break
            # second quadrant
            (x,y) = (-y0,x0)
            if find(((x1,y1),(x2+x,y2+y),(x3,y3),(x4,y4)),vol+det((x,y),(x2,y2))): break
            (x,y) = (-x0,y0)
            if find(((x1,y1),(x2+x,y2+y),(x3,y3),(x4,y4)),vol+det((x2,y2),(x,y))): break
            # third quadrant
            (x,y) = (-x0,-y0)
            if find(((x1,y1),(x2,y2),(x3+x,y3+y),(x4,y4)),vol+det((x,y),(x3,y3))): break
            (x,y) = (-y0,-x0)
            if find(((x1,y1),(x2,y2),(x3+x,y3+y),(x4,y4)),vol+det((x3,y3),(x,y))): break
            # fourth quadrant
            (x,y) = (y0,-x0)
            if find(((x1,y1),(x2,y2),(x3,y3),(x4+x,y4+y)),vol+det((x,y),(x4,y4))): break
            # last octant
            (x,y) = (x0,-y0)
            if find(((x1,y1),(x2,y2),(x3,y3),(x4+x,y4+y)),vol+det((x4,y4),(x,y))): break
        else:
            assert 0 # not found
        bestsol, bestvol = ((x1,y1),(x2,y2),(x3,y3),(x4,y4)),vol
            
    return best, SOL, [apply_mat(MAT,v) for v in CORNERS]
               
def ideal_polygon(Vectors):
    """distribute 1/8-th of each vector on the eight quadrants"""
    assert all(0<=y<=x for x,y in Vectors) # sorted by slope y/x
    Vecs = [(int(2*x),int(2*y)) for (x,y) in Vectors] # integers
    Sx = sum(x for x,y in Vecs)
    Sy = sum(y for x,y in Vecs)
    pointx,pointy = (0,-Sx-Sy) # start on the vertical x-axis
    vol = 0
    for x,y in Vecs:
        vol += det((pointx,pointy),(x,y))
        pointx += x
        pointy += y
    assert pointx==Sx and pointy==-Sx # end on the SE diagonal y=-x
    # vol is the volume (= 2*area) of the 1/8-th central slice of
    # the ideal polygon with all edges scaled by 16
    
    Zauberformel = ( sum (x1*(x2+y2)
                        for i,(x1,y1) in enumerate(Vecs)
                        for j,(x2,y2) in enumerate(Vecs) if i<j) +
                     sum (x*(x+y) for x,y in Vecs)/2 )/32
    # from the formula for centrally symmetric polygons,
    # taking into account the 8 copies of each vector.
    assert Zauberformel == vol/64    
    return vol/64


def construct_fence(SOL, CORNERS, Vectors):
    """construct the fence from the list of vectors in SOL+CORNERS.
    Info about the available polyomino pieces is stored in "Vectors"."""
    global sequencenumber # used for assigning varying colors
    sequencenumber = dict()
    
    VECS = SOL+CORNERS
    # sorted ccw from downward direction
    vecsN = sorted((y/x,x,y) for (x,y) in VECS if x!=0)
    vecs_down = [(x,y) for (x,y) in VECS if x==0 and y<0]
    vecs_up   = [(x,y) for (x,y) in VECS if x==0 and y>0]

    vecs = vecs_down + [(x,y) for _,x,y in vecsN if x>0] +  vecs_up + [(x,y) for _,x,y in vecsN if x<0]

    vector_store = defaultdict(list) # the list of pieces that can be used for each vector.
    #print(Vectors[:])
    for xy,d,o,l,P in Vectors:
        vector_store[tuple(xy)].append((d,o,l,P))
        print(d,o,l,P)
    print("vectors:",vecs)

    GRID = dict()
    # Centers of grid squares are integers.
    # Grid points are half-integer.
    (x,y) = 0,0.5 # start position; first move is in vertical downward direction
    for cc,(u,v) in enumerate(vecs):
        # find the (normalized) distance vector (r,s) for the diameter
        if u%2: # corners
            assert v%2
            r,s = sorted([abs(u)-1,abs(v)-1])
        else:
            assert v%2==0
            r,s = sorted([abs(u),abs(v)])
            s -= 2 # the larger coordinate is reduced by 1
        r//=2; s//=2
        if r==s:
            uv0 = r,
        else:
            uv0 = r,s
        #print((u,v),(r,s),uv0)
        # get next vector from the correct list:
        _,(((x1,y1),(x2,y2)),o),l,P = vector_store[uv0].pop()

        sequencenumber[l]=cc
        print(f"F{cc}. {(u,v)=} {((((x1,y1),(x2,y2)),o),l,P)=}")
        if o<0:
            # "good" side should be on the left, turn by 180° by swapping endpoints
            (x2,y2),(x1,y1) = (x1,y1),(x2,y2)
        # normalize polyomino so that reference vector points into the first quadrant    
        for rot in (rotate0,rotateL,rotate2,rotateR):
            dirx,diry = rot([(x2-x1,y2-y1)])[0]
            if dirx>0 and diry>=0:
                (x1,y1),(x2,y2) = rot([(x1,y1),(x2,y2)])
                P = rot(P)
                break
        else:
            error

        # find the touching-points the normalized polyomino
        #print(f" {(x1,y1)=} {(x2,y2)=} {(u,v)=}, {P=}")
        if u%2: # corners
            xstart,ystart = x1-0.5,y1
            xend,yend     = x2,y2+0.5
        elif x2-x1 <= y2-y1:
            xstart,ystart = x1,y1-0.5
            xend,yend     = x2,y2+0.5
        else:
            xstart,ystart = x1-0.5,y1
            xend,yend     = x2+0.5,y2
           
        dx,dy = xend-xstart, yend-ystart
        #print(f"before: {((xstart,ystart))=}, {(u,v)=} {(dx,dy)=}")
        for (refl,ori2),rot in itertools.product(((rotate0, False),(reflect, True)),
                                         (rotate0,rotateL,rotate2,rotateR)):
            # try all eight symmetries.
            if rot(refl([(dx,dy)]))[0] == (u/2,v/2):
                break
        else:
            assert 0 # error
            
        if ori2:
            P = rotate2(refl(P)) # reflect at line x+y=0
            (xend,yend),(xstart,ystart) = rotate2(refl([(xstart,ystart),(xend,yend)])) # swap endpoints
            assert xstart<=xend and ystart<=yend
            # determine "rot" again:
            dx,dy = xend-xstart, yend-ystart
            #print(f"before2: {((xstart,ystart))=}, {(u,v)=} {(dx,dy)=}")
            for rot in (rotate0,rotateL,rotate2,rotateR):
                if rot([(dx,dy)])[0] == (u/2,v/2):
                    break
            else:
                assert 0 # error


        # carry out the transformation:
        P_transformed = rot(P)
        #print(f" vor {((xstart,ystart),(xend,yend))=}, {(u,v)=}")
        (xstart,ystart),(xend,yend) = rot([(xstart,ystart),(xend,yend)])
        #print(f"FIND {(xend-xstart,yend-ystart)=} {(u,v)=}")
        #print(f" {P_transformed=} {l=}")
        #print(f" nach {((xstart,ystart),(xend,yend))=}, {(u,v)=} {(x,y)=}")
        assert xend-xstart == u/2
        assert yend-ystart == v/2
        assert xstart-x == int(xstart-x)
        assert ystart-y == int(ystart-y)
        for (i,j) in P_transformed:
            kx = int(i-xstart+x)
            ky = int(j-ystart+y)
            assert (kx,ky) not in GRID            
            GRID[kx,ky]=l # mark the grid cell.
        x += u/2
        y += v/2
        #print_grid(GRID, text=f"  {x=} {y=}")
    assert (x,y) == (0,0.5) # returned to start position
    xmid,ymid = print_grid(GRID)
    print("Enclosed area",count_area(GRID, (xmid,ymid)))    

## Incrementally generate the set of n-ominoes
    
LETTERS = {5: "VPULTFZYWNIX", 4: "oltni", # assigned a posteriori after looking at the output
           2: "A ", 3: "LI"}
Defaultletters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"+"0123456789" # arbitrary
assert len(Defaultletters) == 2*26+10

def Letters(n,i):
    if n in LETTERS:
        return  LETTERS[n][i]
    return Defaultletters[i % len(Defaultletters)]

Vectors = []
for n in range(2,n_max+1):
    Vectors = []
    New = set() # find the set of polyominoes of size n
    for P in Polyominoes[n-1]:
        for x,y in P:
            for q in ((x+1,y),(x,y+1),(x-1,y),(x,y-1)):
                if q not in P:
                    P2 = P+(q,)
                    # normalize under rotations and reflections.
                    Q = min(tuple(sorted(normalize(P4))) for
                            P3 in (P2,rotateL(P2),rotateR(P2),rotate2(P2))
                            for P4 in (P3,reflect(P3)))
                    New.add(Q)
    
    New = sorted(tuple(sorted(reflect(Q))) for Q in New) # (reflected is more wide than high.)
    Polyominoes.append(New)
    print(f"\n{n=}. These are the {len(New)} {n}-ominoes:")
    for i,P in enumerate(New): # compute the diametrical pairs
        name = Letters(n,i)
        diam = diameter(P)
        if n==5 and name=="U": # hack
            print(name, i, "hack U-shape ******")
            assert diam ==(10, 2, [((0, 0), (2, 1)), ((0, 1), (2, 0))])
            diam = (9,1,[((0,1),(2,1))])
        if False and n==6 and name=="L": # hack
            # L       # not chosen
            # LLLLL
            print(name, i, "hack L-shape ******")
            assert diam == (26, 1, [((0, 1), (4, 0))])
            diam =         (25, 1, [((0, 0), (4, 0))])
        if False and n==6 and name=="K": # hack
            # K  K
            # KKKK
            print(name, i, "hack U-shape ******")
            assert diam ==(17, 2, [((0, 0), (3, 1)), ((0, 1), (3, 0))])
            diam =        (16, 1, [((0, 1), (3, 1))])
        if n==6 and name=="S": # hack 2
            #  S  
            #  S  
            # SSSS
            print(name, i, "hack T-shape ******")
            assert diam == (16, 1, [((0, 0), (3, 0))])
            diam =         (13, 1, [((1, 2), (3, 0))])
        if n==7 and name=="e" and i+1==31:
            # e     
            # eeeeee
            print(name, i, "hack L-shape ******")
            assert diam == (37, 1, [((0, 1), (5, 0))]) 
            diam = (36, 1, [((0, 0), (5, 0))])
        if False and n==7 and name=="d" and i+1==30:
            # d   d
            # ddddd
            print(name, i, "hack U-shape ******")
            assert diam == (26, 2, [((0, 0), (4, 1)), ((0, 1), (4, 0))])
            diam = (25, 1, [((0, 1), (4, 1))])
        if n==8 and name=="l" and i+1==100:
            # l     
            # lllllll
            print(name, i, "hack L-shape ******")
            assert diam == (50, 1, [((0, 1), (6, 0))])
            diam = (49, 1, [((0, 0), (6, 0))])
        if False and n==8 and name=="k" and i+1==99:
            # k    k
            # kkkkkk
            print(name, i, "hack U-shape ******")
            assert diam == (37, 2, [((0, 0), (5, 1)), ((0, 1), (5, 0))])
            diam = (36, 1, [((0, 1), (5, 1))])
        print(f"#{i+1}, {name=}: diameter {diam}") #, f"{P=}")
        print(f"polyolist[{n},'{name}']={P}")
        if n>2:
            best_d,best_ori = bestposition(*diam)
            print(f"best orientation: deficit={best_d}, best orientation(s)={best_ori}")
        print_polyo(P,letter=name)
        if n>=4:
            vec = set([frozenset([abs(x1-x2),abs(y1-y2)])
                                 for ((x1,y1),(x2,y2)),o in best_ori])
            assert len(vec)==1
            Vectors.append((sorted(list(vec)[0]),best_d,best_ori[0],name,P))
    print()

    if n<=4:
        continue
    Vectors.sort()
    for s in (Vectors):
        print(s)
    deficit =  sum(b for (v,b,o,l,P) in Vectors)
    print(f"Total deficit from all pieces = {deficit}")

    Vecs = []
    if n==5:
        exceptions = 1 # number of pieces that require exceptional treatment
    if n==4:
        exceptions = 3
    if n==6:
        exceptions = 0
    num_ex = 0
    gain = 0 # gain from treating pieces as corners instead of "straight" pieces

    # identify the four corner pieces
    corners = []
    for xy,d,o,l,P in Vectors:
        if (n==6 and l=="S" and xy==[2]): # hack 2
            xy = [2,2]
        if len(xy)==1 and ( # this should be done in preprocessing
                n <7 or (xy[0]==3 and len(corners)<4)):
            x=xy[0]; y=x
            if n==7:
                assert x==3
            if n==5:
                assert x==2
            if n==4:
                assert x==1
            num_ex += 1
            gain += ((x+1)**2-(x+1/2)**2)/2
            x += 1/2
            y += 1/2
            corners.append((x,y))
        else:
            if len(xy)==1:
                xy=xy[0],xy[0]
            x,y=xy
            if ((x,y)==(1,2) or x==0 and n==4) and exceptions:
                exceptions -= 1
                num_ex += 1
                gain += ((x+1)*(y+1)-(x+1/2)*(y+1/2))/2
                x += 1/2
                y += 1/2
                corners.append((x,y))
            else:
                assert x<=y
                y+=1
        Vecs.append((x/y,(x,y),l))
    print(f"{gain=} {num_ex=}")
    assert num_ex == 4
    assert exceptions == 0
    Vecs.sort()
    print (Vecs)
    print(f"{corners=}")
    # the four "diagonal" pieces should cancel out in terms of even/odd
    total_parity_sum = sum(x+y for _,(x,y),_ in Vecs)
    print("total L1-distance must be even:", total_parity_sum)
    assert total_parity_sum%2 == 0 # otherwise solution is impossible.

    VecsR = [(y,x) for _,(x,y),_ in Vecs] # reflect the vectors again
    ideal = ideal_polygon(VecsR)
    print ("""'Ideal' area:""",ideal)
    perimeter = sum(sqrt(x**2+y**2) for _,(x,y),_ in Vecs)
    rad = perimeter/(2*math.pi)
    area = rad**2*math.pi
    print(f"Circle area {area}, +deficit+gain={area+deficit+gain}")
    perimeter += 1
    rad = perimeter/(2*math.pi)
    area2 = rad**2*math.pi
    print(f"With length+1: Circle area {area2}, delta={area2-area}")

    if n==77:
        continue    
    if 0 and n==6:
        # solution found in a different run: (runtime about 20 hours)
        SOL = [ (0, -8), (0, 8), (0, -8), (0, 8), (-8, 0), (10, 0), (-10,0),
                (12, 0), (2, -10), (-10, 2), (10, -2), (-2, 8), (2, -8), (-2, 8),
                (-2,-8), (2, 8), (-2, -8), (2, 8), (-8, -2), (8, 2), (-8, -2), (-2, -6),
                (6, 2), (-6, 2), (-2, 6), (4, -8),
                (4, 8), (-8, -4), (8, 4), (-8, 4), (8, -4)]
        CORNERS=[(5, 5), (-5, 5), (-5, -5), (5, -5)]
        construct_fence(SOL, CORNERS, Vectors)
        print (f"""'Ideal' area, an upper bound (given the edge vectors): {ideal}, """
               f"+deficit+gain={ideal+deficit+gain}")
        continue
        exit() #

    B, SOL, CORNERS = best_polygon(VecsR)
    print(f"best solution value {B/8}. **** {B/8+deficit+gain=} ************")
    print (f"""'Ideal' area, likely an upper bound: {ideal}, """
           f"+deficit+gain={ideal+deficit+gain}")        
    print()
    print(f"{SOL=}; {CORNERS=}")
    print(f"{sum(x for x,y in SOL)=}; {sum(y for x,y in SOL)=}")
    construct_fence(SOL, CORNERS, Vectors)
    stdout.flush()
