"""
Generate all reference vectors for n-ominoes.
Usage: python best_fence_step0.py <n> [--progress] [--verbose]

With --progress, a progress counter is shown while the boundary vectors are computed.
With --verbose, the pieces are printed as a pattern of X's

The result are two lists, each containing a dictionary of reference vectors for each piece.

Typical dictionary entries looks like this:

  ((1.0, 4.0),   0, 1): 3.5,
  ((2.5, 1.5),   1, 1): -0.125, 
  ((0.0, -3.0), -2, 0): 6.0,

which means: There is a vector (0,-3), with rotation D=-2, parity 0, and penalty 6.
Only the smallest penalty is kept for each entry

The list "vectors":
   vectors are folded by reflection into the octant 0<=x<=y
The list "vectors_precise":
   vectors as they are generated. This is needed for the reconstruction
   of the solution in the end.
"""

import math
import sys
from collections import defaultdict

_positional = [a for a in sys.argv[1:] if not a.startswith("--")]
n_max = int(_positional[0]) if _positional else 7
PROGRESS = "--progress" in sys.argv  # off by default
VERBOSE = "--verbose" in sys.argv  # off by default
DIRECTIONS = (0,-1),(1,0),(0,1),(-1,0) # 0=South, 1=East, 2=North, 3=West
# counterclockwise directions

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

def print_piece(piece, text="", letter="X"):    
    # print the pattern represented in the list piece
    xmin = min(x for x,y in piece)
    xmax = max(x for x,y in piece)
    ymin = min(y for x,y in piece)
    ymax = max(y for x,y in piece)
    pattern = [["." for _ in range(1+xmax-xmin)] for _ in range(1+ymax-ymin)]
    for (x,y) in piece:
        pattern[ymax-y][x-xmin]=letter
    print("\n".join("".join(l) for l in pattern) + text)
    print()


## COMPUTE POLYOMINOES   
## Incrementally generate the set of n-ominoes, adding one cell at a time
    
import string as _s
# Generate enough labels for any number of pieces (covers 369 octominoes and beyond):
_base = _s.ascii_uppercase + _s.ascii_lowercase + _s.digits + _s.punctuation
Letters6 = list(_base) + [f"#{i}" for i in range(len(_base), 500)]

Polyominoes = [[],[((0,0),)]] # start data, size 0 and size 1
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
    print(f"n={n}: {len(New)} {n}-ominoes.", flush=True)
    
Pieces = New

if VERBOSE:
    for i,P in enumerate(Pieces, start=1):
        print_piece(P, text=f" piece {i} of {len(Pieces)}")


## COMPUTE BOUNDARIES. (with n>=7 cells, the boundary is not always a simple polygon.)
vectors = []
vectors_precise = []
vmax = []
Kmax = []
_nPieces = len(Pieces)
print(f"Computing boundary vectors for {_nPieces} pieces...", flush=True)
for iP,P in enumerate(Pieces):
    if PROGRESS and iP % max(1, _nPieces//100) == 0:
        print(f"\r  piece {iP+1}/{_nPieces} ({100*(iP+1)//_nPieces}%)", end="", flush=True)
    boundary_successor=defaultdict(list)
    for (x,y) in P:
        for (dx,dy) in DIRECTIONS:
            xx,yy = x+dx,y+dy
            if (xx,yy) not in P:
                perpx,perpy = -dy/2,dx/2
                start,end=(((x+xx)/2+perpx, (y+yy)/2+perpy),
                           ((x+xx)/2-perpx, (y+yy)/2-perpy))
                # clockwise around the interior
                boundary_successor[start].append(end)
    start = current = min(boundary_successor) # start on the boundary
    path = [start]
    assert len(boundary_successor[current])==1
    current = boundary_successor[current][0]
    while current!= start:
        path.append(current)
        if len(boundary_successor[current])==1:
            current = boundary_successor[current][0]
        else:
            # non-simple polygon:
            # Choose the right turn. (non-full holes remain outside)
            dx,dy = current[0]-path[-2][0], current[1]-path[-2][1]
            dx,dy = dy,-dx # right by 90°
            next = current[0]+dx, current[1]+dy
            assert next in boundary_successor[current]
            current = next
    #assert len(path) == len(boundary_successor)
    L = len(path)
    ## BOUNDARY HAS BEEN COMPUTED.
    
    path0 = path
    #print(L,P,path) # CLOCKWISE boundary traversal
    ref_vectors = dict()
    ref_vectors_precise = dict()
    for path in (path0,reflect(reversed(path0))):
        path2 = path+path
        area_under = []
        directions = []
        for (x1,y1),(x2,y2) in zip(path,path2[1:]):
            area_under.append( (x2-x1)*(y1+y2) ) # = 2 x the signed area on the RIGHT side
            directions.append(DIRECTIONS.index(((x2-x1),(y2-y1))))
        area_under += area_under # double the sequence (two rounds)
        sum_area_under = [0]
        s = 0
        for a in area_under:
            s += a
            sum_area_under.append(s)
        # sum_area_under[i] = area_under[0] + area_under[1] + ... + area_under[i-1]
        assert sum_area_under[L] >= 2*n_max # can have nonsimple polygons
        
        directions *= 2
        absolute_direction = [0]
        dir = 0
        for dir1,dir2 in zip(directions,directions[1:]):
            dir += {3:-1, 0:0, 1:1}[(dir2-dir1) % 4]
            absolute_direction.append(dir)
        # absolute_direction[i]= direction of edge Pi,P(i+1)
        # counterclockwise from edge P0,P1
        assert absolute_direction[L]==-4 # one full clockwise traversal
        assert len(absolute_direction)==2*L
        #print(iP,absolute_direction)
        for i in range(L):
            for j in range(i+1,i+L):           
                # first edge i,i+1
                # last edge j,j+1
                (x1,y1),(x2,y2) = path2[i:i+2]
                xS,yS = (x1+x2)/2, (y1+y2)/2
                dirS = DIRECTIONS.index(((x1-x2),(y1-y2))) # from i+1 to i
                (x1,y1),(x2,y2) = path2[j:j+2]
                xE,yE = (x1+x2)/2, (y1+y2)/2
                dirE = DIRECTIONS.index(((x2-x1),(y2-y1))) # from j to j+1

                (x2,y2) = path2[i+1]
                parity_x = x1-x2
                parity_y = y1-y2
                assert parity_x==int(parity_x)
                assert parity_y==int(parity_y)
                parity = int(parity_x+parity_y) % 2
                # parity of x+y for the "inner" reference vector
                
                along_boundary = (area_under[i]+area_under[j])/2 + sum_area_under[j]- sum_area_under[i+1]
                a_under_ref_vector = (xE-xS)*(yE+yS) # 2 * area under reference vector
                penalty = (along_boundary-a_under_ref_vector)/2 # divide by 2, come back to true area
                rot_old = (dirE-dirS) % 4
                # We let "rot" be an absolute signed quantity,
                # not just a change of directions like rot_old.
                # The total rotation of all selected pieces must be +4.
                # (This makes the bounds tighter.)
                rot = absolute_direction[j]-absolute_direction[i] + 2
                assert (rot-rot_old)%4 == 0
                [v] = rotate([(xE-xS, yE-yS)], -dirS)
                  # reference vector with start-edge = DOWN
                key = (v,rot,parity)
                if key not in ref_vectors_precise or penalty<ref_vectors_precise[key]:
                    ref_vectors_precise[key] = penalty
                  
                # AT THIS LEVEL, ORIENTATION DOES NOT MATTER:
                # We are just computing bounds based on the length and slope.
                # (It matters only when actually putting together the vectors.)
                v = tuple(sorted(map(abs,v)))
                
                key = (v,rot,parity)
                if key not in ref_vectors or penalty<ref_vectors[key]:
                    ref_vectors[key] = penalty
            
    # Null option (v, rot, penalty) = (0, 0, 0): lets the piece be unused.
    ref_vectors[((0,0),0,0)] = 0
    ref_vectors_precise[((0,0),0,0)] = 0

    vectors.append(ref_vectors)
    vectors_precise.append(ref_vectors_precise)
    minPen = min(ref_vectors.values())
    maxLen= max(math.sqrt(x**2+y**2) for (x,y),_,_ in ref_vectors)
    vmax.append(maxLen)
    Kmax.append(-minPen)
    #print("piece","%2d"%iP, Letters6[iP], P, len(ref_vectors),"reference vectors",
          #f"({len(ref_vectors_precise):3} precise)",
          #f"min.pen. {minPen:<6.3g}, max.length {maxLen:.7}")
if PROGRESS:
    print(end="\r  ")          
print(f"Done: {_nPieces} pieces, {sum(len(r) for r in vectors)} total reference vectors.", flush=True)
#print("total number of reference vectors:", sum(len(r) for r in vectors))
