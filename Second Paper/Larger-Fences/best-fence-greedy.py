# Greedy assignment of vectors to the eight arcs.
# in the end, vector sum 0 is reached by local swappings

"""HACK of the pieces:
variation with the T-shaped piece (code "S") or the C-shaped piece (code "K").
Look for "hack" in the program.
"""

import argparse
import itertools
from collections import defaultdict
import sys
import os
import gc
import time

parser = argparse.ArgumentParser(
    description="Greedy construction of a fence from all n-ominoes."
)
parser.add_argument("n", nargs="?", type=int, default=9,
                    help="size of the polyominoes (default: 9)")
parser.add_argument("--no-csv", action="store_true",
                    help="do not write the solution grid to data/sol-*.csv")
parser.add_argument("--no-area-check", action="store_true",
                    help="skip the sweep-line cross-check of the enclosed areas")
parser.add_argument("--no-piece-check", action="store_true",
                    help="skip the check that every n-omino is used exactly once "
                         "and correctly")
parser.add_argument("--no-shoelace", action="store_true",
                    help="skip the boundary/shoelace area computation; the first area is "
                         "then taken from the formula B/8+deficit+gain and the total area "
                         "is that plus total_inner_holes")
parser.add_argument("-v", "--verbose", action="store_true",
                    help="also print every piece with inner holes / ambiguous "
                         "diameters, and the debugging grids")
args = parser.parse_args()
n_max = args.n

Polyominoes = [[], [((0, 0),)]]


def vprint(*a, **k):
    """print only with --verbose"""
    if args.verbose:
        print(*a, **k)


_last_progress = [0.0]


def print_progress(current, total, prefix="Progress:", suffix="Complete", length=50):
    """Carriage-return progress bar (redraws at most ~10 times per second)."""
    now = time.monotonic()
    if current < total and now - _last_progress[0] < 0.1:
        return
    _last_progress[0] = now
    percent = ("{0:.1f}").format(100 * (current / float(total)))
    filled_length = int(length * current // total)
    bar = "█" * filled_length + "-" * (length - filled_length)
    sys.stdout.write(f"\r{prefix} |{bar}| {percent}% {suffix}")
    sys.stdout.flush()
    if current == total:
        sys.stdout.write("\n")
        sys.stdout.flush()


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
    xmin = min(x for x, y in GRID.keys())
    xmax = max(x for x, y in GRID.keys())
    ymin = min(y for x, y in GRID.keys())
    ymax = max(y for x, y in GRID.keys())
    pattern = [["." for _ in range(1 + xmax - xmin)] for _ in range(1 + ymax - ymin)]
    for (x, y), letter in GRID.items():
        pattern[ymax - y][x - xmin] = str(letter)
    print(
        "\n".join("".join(l) for l in pattern)
        + text
        + f",  {xmin}<=x<={xmax}, {ymin}<=y<={ymax}"
    )

    return (
        (xmin + xmax) // 2,
        (ymin + ymax) // 2,
    )  #  midpoint, hopefully an interior point

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
    dpairs = [tuple(sorted(((x,y),(x2,y2))))
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
        # mark from (x1,y1) to the left boundary
        # and mark from (x2,y2) to the right boundary
        # for a diameter, these paths must be free
        bars = [(x,y1) for x in range(-1,x1)] + [(x,y2) for x in range(x2+1,2+xmax)]
        bars2 = [(x,y1) for x in range(-1,x1+1)] + [(x,y2) for x in range(x2,2+xmax)]
        for x,y in bars:
            assert cell[x,y]==" "
            cell[x,y]="Q"
        todo_upper = [(x,y+1) for x,y in bars2 if cell[x,y+1]==" "]
        todo_lower = [(x,y-1) for x,y in bars2 if cell[x,y-1]==" "]
            
        # look for the two components of the boundary, and extend from
        # there into the interior
        below = (x2-x1+1)*(y1+y2+1)/2  # area below diameter inside the rectangle [x1,x2] x [0,ymax]
        above = (x2-x1+1)*(ymax+1) - below
        for num_comp,todo,ori0 in [("1",todo_upper,+1), ("2",todo_lower,-1), ]:
            size = sum(1 for c in todo if c not in boundary)
            for x0,y0 in todo: cell[x0,y0] = num_comp
            while todo:
                (x,y) = todo.pop()
                for (dx,dy) in neighbors8:
                    xx,yy = x+dx,y+dy
                    if (-1<=xx<=xmax+1 and -1<=yy<=ymax+1 and cell[xx,yy]==" "):
                        cell[xx,yy]= num_comp
                        todo.append((xx,yy))
                        if (xx,yy) not in boundary:
                            size += 1

            if 0:
                    print (P)
                    print (f"{dpair=} {above=} {below=}")
                    print(x0,y0,ori0)
                    for yy in reversed(range(-1,ymax+2)):
                        print("".join(str(cell[xx,yy]) for xx in range(-1,xmax+2)))
            if ori0 == -1: # below the diameter
                size = size - below - x1*(y1) - (xmax-x2)*(y2)  # the net contribution
                # the rectangles below the horizontal extensions are extra
            else:
                size = size - above - x1*(ymax-y1) - (xmax-x2)*(ymax-y2)
            # remember the best solution(s)
            if best_part is None or size>best_part:
                best_part = size
                best_sol = [(dpair,ori0)]
            elif size==best_part:
                best_sol.append((dpair,ori0))
                #print(f"{num_comp=}, {size=}")
        assert all(cell[x0,y0]!=" " for x0,y0 in boundary)
        if 0:
                    print (P)
                    print (f"{dpair=}")
                    print(x0,y0)
                    for yy in reversed(range(-1,ymax+2)):
                        print("".join(str(cell[xx,yy]) for xx in range(-1,xmax+2)))
                    exit()
    
    innerholes = sum(1 for i in range(1+xmax)
                     for j in range(1+ymax) if cell[i,j]==" ")
    if innerholes and args.verbose:
        #print(f"{((x1,y1),(x2,y2))=}"),
        for j in range(1+ymax,-2,-1):
            for i in range(-1,2+xmax):
                print(end=f"{cell[i,j]:1}")
            print()
    if best_part is None:
        return
    return best_part, best_sol, innerholes


def poly_canon(P):
    """An optimized, low-allocation version of poly_canon"""
    min_variant = None
    for sym in range(8):
        if sym == 0:
            coords = P
        elif sym == 1:
            coords = [(-y, x) for x, y in P]
        elif sym == 2:
            coords = [(-x, -y) for x, y in P]
        elif sym == 3:
            coords = [(y, -x) for x, y in P]
        elif sym == 4:
            coords = [(y, x) for x, y in P]
        elif sym == 5:
            coords = [(-x, y) for x, y in P]
        elif sym == 6:
            coords = [(-y, -x) for x, y in P]
        else:
            coords = [(x, -y) for x, y in P]

        xmin = min(c[0] for c in coords)
        ymin = min(c[1] for c in coords)
        normalized = tuple(sorted((c[0] - xmin, c[1] - ymin) for c in coords))
        if min_variant is None or normalized < min_variant:
            min_variant = normalized
            
    # Reflect the final absolute minimum
    return tuple(sorted(reflect(min_variant)))


def calculate_enclosed_areas_boundary(occ):
    if not occ:
        return []

    # Helper functions inline for performance
    def pack(x, y):
        return (x << 32) | (y & 0xFFFFFFFF)

    def unpack(val):
        x = val >> 32
        y = val & 0xFFFFFFFF
        if y >= 0x80000000:
            y -= 0x100000000
        return x, y

    # Build adjacency dictionary directly (replaces edges set allocation)
    adj = defaultdict(list)
    for x, y in occ:
        if (x, y - 1) not in occ:
            adj[pack(x, y)].append(pack(x + 1, y))
        if (x + 1, y) not in occ:
            adj[pack(x + 1, y)].append(pack(x + 1, y + 1))
        if (x, y + 1) not in occ:
            adj[pack(x + 1, y + 1)].append(pack(x, y + 1))
        if (x - 1, y) not in occ:
            adj[pack(x, y + 1)].append(pack(x, y))

    areas = []
    total_vertices = len(adj)

    # Destructive, in-place traversal of the adjacency map to conserve memory
    while adj:
        processed_vertices = total_vertices - len(adj)
        if processed_vertices % 5000 == 0 or not adj:
            print_progress(processed_vertices, total_vertices, prefix="Verification Shoelace:", suffix=f"({processed_vertices}/{total_vertices})")

        start_id = next(iter(adj))

        cycle = []
        curr_id = start_id
        p0_id = None

        while True:
            cycle.append(curr_id)

            candidates = adj.get(curr_id)
            if not candidates:
                break

            if len(candidates) == 1:
                next_id = candidates[0]
                del adj[curr_id]
            else:
                # We have a choice, apply preference rule
                cx, cy = unpack(curr_id)

                if p0_id is None:
                    next_id = candidates[0]
                else:
                    px, py = unpack(p0_id)
                    d_in = (cx - px, cy - py)

                    # Preferences based on d_in
                    left = (-d_in[1], d_in[0])
                    straight = d_in
                    right = (d_in[1], -d_in[0])
                    backward = (-d_in[0], -d_in[1])

                    # Map candidate to its direction from curr_id
                    cand_map = {}
                    for c_id in candidates:
                        tx, ty = unpack(c_id)
                        cand_map[(tx - cx, ty - cy)] = c_id

                    next_id = None
                    for d_pref in (left, straight, right, backward):
                        if d_pref in cand_map:
                            next_id = cand_map[d_pref]
                            break
                    if next_id is None:
                        next_id = candidates[0]

                # Remove selected candidate and clean up key if empty
                candidates.remove(next_id)
                if not candidates:
                    del adj[curr_id]

            p0_id = curr_id
            curr_id = next_id

            if curr_id == start_id:
                break

        # Compute area using shoelace formula directly on packed integers in cycle
        if len(cycle) >= 3:
            cycle.append(cycle[0])
            shoelace_sum = 0.0
            x_prev, y_prev = unpack(cycle[0])
            for i in range(1, len(cycle)):
                x_curr, y_curr = unpack(cycle[i])
                shoelace_sum += x_prev * y_curr - x_curr * y_prev
                x_prev, y_prev = x_curr, y_curr

            area = shoelace_sum / 2.0
            areas.append(area)

    print_progress(total_vertices, total_vertices, prefix="Verification Shoelace:", suffix="Complete")

    # Sort enclosed areas (holes have negative areas)
    enclosed_areas = []
    for a in areas:
        if a < 0.0:
            enclosed_areas.append(round(abs(a)))
    enclosed_areas.sort(reverse=True)
    return enclosed_areas


def calculate_enclosed_areas_sweep(occ):
    """Enclosed areas via a scan-line sweep with union-find.

    Free cells are 8-connected (the fence is 4-connected).  No bounding-box
    grid is built.  Each row is described by its
    free *intervals* (the gaps between runs of fence cells).  Sweeping row by
    row, an interval is merged with the intervals of the previous row that it
    touches (diagonal contact counts).  The intervals touching the left/right
    border of the bounding box (and the rows above/below) belong to the
    outside.  Every other union-find component is an enclosed region.

    Time O(M log M), memory O(M) for M fence cells, independent of the
    size of the enclosed regions.
    """
    if not occ:
        return []

    rows = {}
    for x, y in occ:
        rows.setdefault(y, []).append(x)
    lo = min(x for x, _ in occ) - 1  # one column of padding on each side
    hi = max(x for x, _ in occ) + 1
    ymin, ymax = min(rows), max(rows)

    # union-find over free intervals; node 0 = the outside
    parent = [0]
    area = [0]

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]  # path halving
            a = parent[a]
        return a

    def union(a, b):
        a, b = find(a), find(b)
        if a != b:
            if area[a] < area[b]:
                a, b = b, a
            parent[b] = a
            area[a] += area[b]

    def free_intervals(y):
        """sorted, disjoint list of (start, end, node) of the free cells in
        row y.  Consecutive intervals are separated by >= 1 fence cell."""
        xs = rows.get(y)
        if not xs:
            return [(lo, hi, 0)]
        xs.sort()
        out = [(lo, xs[0] - 1, 0)]  # left of the first fence cell: outside
        run_end = xs[0]
        for x in xs[1:]:
            if x != run_end + 1:
                a, b = run_end + 1, x - 1
                parent.append(len(parent))
                area.append(b - a + 1)
                out.append((a, b, len(parent) - 1))
            run_end = x
        out.append((run_end + 1, hi, 0))  # right of the last fence cell: outside
        return out

    def link(prev, cur):
        """union all intervals of two adjacent rows that touch (8-connectivity)"""
        i = j = 0
        while i < len(prev) and j < len(cur):
            a, b, ia = prev[i]
            c, d, ic = cur[j]
            if a <= d + 1 and c <= b + 1:
                union(ia, ic)
            # intervals within a row are >= 1 apart, so the one that ends
            # first cannot touch anything further along in the other row
            if b < d:
                i += 1
            elif b > d:
                j += 1
            else:
                i += 1
                j += 1

    outside_row = [(lo, hi, 0)]
    prev = outside_row  # the row above the fence is entirely outside
    for y in range(ymin, ymax + 1):
        cur = free_intervals(y)
        link(prev, cur)
        prev = cur
    link(prev, outside_row)  # the row below the fence

    root0 = find(0)
    areas = [area[k] for k in range(1, len(parent))
             if parent[k] == k and k != root0]
    areas.sort(reverse=True)
    return areas


def is_connected(occ):
    if not occ:
        return True
    visited = set()
    start = next(iter(occ))
    queue = [start]
    visited.add(start)
    head = 0
    while head < len(queue):
        cx, cy = queue[head]
        head += 1
        for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nb = (cx + dx, cy + dy)
            if nb in occ and nb not in visited:
                visited.add(nb)
                queue.append(nb)
    return len(visited) == len(occ)


def check_pieces(occ, n, expected_pieces):
    """Check that the fence consists of exactly the expected n-ominoes:
    each piece id forms one connected component of n cells, and every
    n-omino shape occurs exactly once.  Returns a list of error messages."""
    exp_set = set(expected_pieces)
    vis = set()
    comps = {}
    errors = []

    for rc, pid in occ.items():
        if rc in vis:
            continue
        comp = []
        queue = [rc]
        vis.add(rc)
        head = 0
        while head < len(queue):
            cxy = queue[head]
            head += 1
            comp.append(cxy)
            for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                nb = (cxy[0] + dx, cxy[1] + dy)
                if occ.get(nb) == pid and nb not in vis:
                    vis.add(nb)
                    queue.append(nb)

        if pid not in comps:
            comps[pid] = comp
        else:
            errors.append(f"Piece {pid}: multiple components")
            comps[pid].extend(comp)

        if len(comp) != n:
            errors.append(f"Piece {pid}: {len(comp)} cells, expected {n}")

    found = {}
    for pid, cells in comps.items():
        if len(cells) != n:
            continue
        pc = poly_canon(cells)
        found[pc] = found.get(pc, 0) + 1

    missing_count = sum(1 for p in exp_set if p not in found)
    if missing_count > 0:
        errors.append(f"{missing_count} n-omino(es) missing")

    dup_count = sum(1 for v in found.values() if v > 1)
    if dup_count > 0:
        errors.append(f"{dup_count} n-omino(es) duplicated")

    invalid_count = sum(1 for k in found if k not in exp_set)
    if invalid_count > 0:
        errors.append(f"{invalid_count} invalid shapes")

    return errors


def area_freq(areas):
    """compact summary of a list of areas: {area: multiplicity}, largest area first"""
    freq = {}
    for a in areas:
        freq[a] = freq.get(a, 0) + 1
    return dict(sorted(freq.items(), reverse=True))


def verify_solution(occ, n, N, expected_pieces,
                    do_pieces=True, do_areas=True, do_shoelace=True,
                    expected_first=None, expected_total=None):
    """Returns (passed, first_area, total_area, message)."""
    errors = []

    if do_pieces:
        errors += check_pieces(occ, n, expected_pieces)
    else:
        print("  [Validation Verification] SKIPPED: piece check (--no-piece-check)")

    if not is_connected(occ):
        errors.append("Fence not 4-connected")

    if do_shoelace:
        areas = calculate_enclosed_areas_boundary(occ)
        first_a = areas[0] if areas else 0
        total_a = sum(areas)
    else:
        areas = None
        first_a, total_a = expected_first, expected_total
        print("  [Validation Verification] SKIPPED: shoelace (--no-shoelace); "
              f"using formula: first={first_a} total={total_a}")

    if do_areas:
        areas_sweep = calculate_enclosed_areas_sweep(occ)
        if areas is not None:
            if areas == areas_sweep:
                print(
                    f"  [Validation Verification] SUCCESS: Boundary Shoelace and sweep-line results MATCH: {area_freq(areas)}"
                )
            else:
                print("  [Validation Verification] ERROR: Area MISMATCH!")
                print(f"    Boundary Shoelace result: {area_freq(areas)}")
                print(f"    Sweep-line result:        {area_freq(areas_sweep)}")
                errors.append("Area discrepancy detected between boundary and sweep-line checks")
        else:
            sw_first = areas_sweep[0] if areas_sweep else 0
            sw_total = sum(areas_sweep)
            if (sw_first, sw_total) == (first_a, total_a):
                print("  [Validation Verification] SUCCESS: sweep-line agrees with formula: "
                      f"first={sw_first} total={sw_total}")
            else:
                print("  [Validation Verification] ERROR: sweep-line differs from formula!")
                print(f"    formula:    first={first_a} total={total_a}")
                print(f"    sweep-line: first={sw_first} total={sw_total}")
                errors.append("Area discrepancy detected between formula and sweep-line check")
    else:
        print("  [Validation Verification] SKIPPED: area cross-check (--no-area-check)")

    num_pieces = len(set(occ.values()))

    if not errors:
        return (
            True, first_a, total_a,
            f"Verification PASSED  first={first_a}  total={total_a}  pieces={num_pieces}",
        )
    return (
        False, first_a, total_a,
        "Verification FAILED:\n" + "\n".join("  ERROR: " + err for err in errors),
    )


def write_csv(path, occ, N):
    if not occ:
        return False
    xs = [c[0] for c in occ.keys()]
    ys = [c[1] for c in occ.keys()]
    xmn = min(xs)
    ymn = min(ys)
    W = max(xs) - xmn + 1
    H = max(ys) - ymn + 1

    if H * W > 10**6: # 2_000_000_000:
        print(
            f"Grid {W}x{H} ({H*W // 1_000_000}M cells) too large — CSV skipped."
        )
        return False

    grid = [[0] * W for _ in range(H)]
    for (x, y), pid in occ.items():
        grid[y - ymn][x - xmn] = pid

    dir_name = os.path.dirname(path)
    if dir_name:
        os.makedirs(dir_name, exist_ok=True)

    w = len(str(N))
    with open(path, "w") as f:
        for r in range(H):
            row_strs = [str(grid[r][c]).zfill(w) for c in range(W)]
            f.write(",".join(row_strs) + "\n")

    print(f"Written {W}x{H} -> {path}")
    return True


def best_polygon(Vecs):
    """The workhorse of this algorithm.
    Vecs (x,y) satisfy 0<=y<=x and are sorted by slope."""
    
    global angle_area_store, area_threshold, sum_fixed_vectors

    Vecs = [(int(2*x),int(2*y)) for (x,y) in Vecs] # integers, doubled

    num = len(Vecs)

    # By convention, the vectors that are to be used as corners
    # are the last 4 vectors. (Those that are closest to the 45° direction.)
    corners = Vecs[-4:]
    assert all(s%2 for v in corners for s in v)
    expected_corner_length = n - 1 + (n%2)
    if n>=6:
        assert set(corners) == {(expected_corner_length,expected_corner_length)}
    else:
        assert set([v for v in corners[1:]]) == {(5, 5)}

    v1, v2, v3, v4 = corners
    (x1, y1) = v1
    (x2, y2) = v2
    x2 = -x2  # rotate the corner vectors into the 4 quadrants.
    (x3, y3) = v3
    x3 = -x3
    y3 = -y3
    (x4, y4) = v4
    y4 = -y4
    vol = 0

    CORNERS = [v1, (x2, y2), (x3, y3), (x4, y4)]

    assert sum(x for x, y in CORNERS) == 0
    assert sum(y for x, y in CORNERS) == 0

    histo = defaultdict(int)
    for x,y in Vecs[:-4]:
        assert x%2 == y%2 == 0
        histo[x//2,y//2] += 1
    ym = max(y for x,y in histo)
    hmax = max(histo.values())
    wid = len(str(hmax))
    print("histogram of vector frequencies")
    for y in range(ym,-1,-1):
        for x in range(0,1+max(x for x,y0 in histo if y0==y)):
            if x<=y and y>0:
                assert (x,y) not in histo
                print(end=" "*(wid+1))
            elif (x,y) not in histo:
                print(end=" "*wid+"0")
            else:
                print(end=f"{histo[x,y]:{wid+1}}")
        print()

    def dirs(v):
        (x0,y0) = v
        for i,v2 in enumerate(((x0,y0),(y0,x0),(-y0,x0),(-x0,y0),(-x0,-y0),(-y0,-x0),(y0,-x0),(x0,-y0))):
            yield v2,i

    vecs = [(x**2+y**2,x,y) for x,y in histo]
    vecs.sort(reverse=True) # from longest to shortest
    #print(vecs)
    sums = [[0,0] for _ in range(8)]
    sx,sy = 0,0
    sol = []
    #from random import shuffle
    #shuffle(vecs)
    for _,x0,y0 in vecs:
        mult = histo[x0,y0]
        for _ in range(mult):
            options = []
            for (x,y),i in dirs((x0,y0)):
                sumsi = sums[i]
                sumsi[0]+=x0
                sumsi[1]+=y0
                balance = (max(x for x,y in sums)-min(x for x,y in sums)+ 
                           max(y for x,y in sums)-min(y for x,y in sums))
                options.append((balance, # primary criterion
                                (sx+x)**2+(sy+y)**2, # sum close to 0, secondary criterion
                                (sx+x,sy+y),i,(x,y)))
                sumsi[0]-=x0
                sumsi[1]-=y0
            
            _,_,(sx,sy),i,(x,y) = min(options)
            sumsi = sums[i]
            sumsi[0] += x0
            sumsi[1] += y0
            sol.append((x, y))
    assert sum(x for x, y in sol) == sx
    assert sum(y for x, y in sol) == sy
    print("the eight arcs", sums)
    print("total", (sx, sy))
    assert (sx + sy) % 2 == 0
    # reduce the total sum vector to zweo
    while abs(sx) >= 2 or abs(sy) >= 2:
        if abs(sx) >= 2:
            if sx > 0:
                sign = +1
            else:
                sign = -1
            def partner(u,v):
                return (-u+sign,v),(u-sign,v),(-u,v)
        elif abs(sy)>=2:
            if sy>0:
                sign = +1
            else:
                sign = -1
            def partner(u,v):
                return (u,-v+sign),(u,v-sign),(u,-v)
        # find two partners that swap places
        for i,(u,v) in enumerate(sol):
            pa,sw1,sw2 = partner(u,v)
            try:
                j = sol.index(pa)
                break
            except ValueError:
                pass
        else:
            print("no partner found")
            exit(1)
        print("swap out", sol[i], sol[j])
        sol[i] = sw1
        sol[j] = sw2
        sx = sum(x for x, y in sol)
        sy = sum(y for x, y in sol)
        print("swapped in:", sol[i], sol[j], ", sum=", (sx, sy))
    if (sx, sy) != (0, 0):
        assert abs(sx) == abs(sy) == 1
        # final step: reflect a piece
        for i, (u, v) in enumerate(sol):
            if abs(abs(u) - abs(v)) != 1:
                continue
            if u > 0:
                unew = abs(v)  # sign stays, absolute values are swapped
            else:
                unew = -abs(v)
            if v > 0:
                vnew = abs(u)
            else:
                vnew = -abs(u)
            if sx + unew - u == 0 == sy + vnew - v:
                break
        else:
            print("no vector found")
            exit(1)
        print("reflect",sol[i])
        sol[i] = unew,vnew
        sx = sum(x for x,y in sol)
        sy = sum(y for x,y in sol)
        print("reflected",sol[i],", sum=",(sx,sy))
                
    global SOL # list of vectors used in the solution (in addition to the corner vectors)
    SOL = [(2*x,2*y) for x,y in sol]
    VECS = SOL+CORNERS
    # sorted ccw from downward direction
    vecsN = sorted((y/x,x,y) for (x,y) in VECS if x!=0)
    vecs_down = [(x,y) for (x,y) in VECS if x==0 and y<0]
    vecs_up   = [(x,y) for (x,y) in VECS if x==0 and y>0]

    vecs = vecs_down + [(x,y) for _,x,y in vecsN if x>0] +  vecs_up + [(x,y) for _,x,y in vecsN if x<0]

    vol = pointx = pointy = 0
    for x,y in vecs:
        vol += det((pointx,pointy),(x,y))
        pointx += x
        pointy += y
    assert pointx==sx and pointy==sy
    # vol is 2*area

    best = vol
    print(f"{n_max=},", best,"best value", flush=True)
    return best, SOL,  CORNERS


def save_fence_visualization(path_vertices, n, filename=None):
    """
    Saves fence path to txt to allow PNG visualization later
    """
    path_dir = "data"
    os.makedirs(path_dir, exist_ok=True)
    out_path = os.path.join(path_dir, f"path-{n}.txt")
    try:
        with open(out_path, "w") as f:
            for x, y in path_vertices:
                f.write(f"{x},{y}\n")
        file_size_kb = os.path.getsize(out_path) // 1024
        print(f"-> Saved scale-invariant path coordinates to {out_path} ({file_size_kb:,} KB)")
        print(f"-> To render the image using your system python, run:")
        print(f"   python3 plot-path.py {out_path} {n}")
    except Exception as e:
        print(f"Could not save path coordinates: {e}")


def construct_fence(SOL, CORNERS, Vectors, expected_first=None, expected_total=None):
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
        #print(d,o,l,P)
    if n<=7: print("vectors:",vecs)

    GRID = dict()
    # Centers of grid squares are integers.
    # Grid points are half-integer.
    (x, y) = 0, 0.5  # start position; first move is in vertical downward direction
    path_vertices = [(x, y)]

    for cc, (u, v) in enumerate(vecs):
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
        _, (((x1, y1), (x2, y2)), o), l, P = vector_store[uv0].pop()

        sequencenumber[l] = cc
        if o < 0:
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
            assert 0  # error

        if ori2:
            P = rotate2(refl(P)) # reflect at line x+y=0
            (xend,yend),(xstart,ystart) = rotate2(refl([(xstart,ystart),(xend,yend)])) # swap endpoints
            assert xstart<=xend and ystart<=yend
            # determine "rot" again:
            dx, dy = xend - xstart, yend - ystart
            for rot in (rotate0, rotateL, rotate2, rotateR):
                if rot([(dx, dy)])[0] == (u / 2, v / 2):
                    break
            else:
                assert 0 # error


        # carry out the transformation:
        P_transformed = rot(P)
        (xstart, ystart), (xend, yend) = rot([(xstart, ystart), (xend, yend)])
        assert xend - xstart == u / 2
        assert yend - ystart == v / 2

        if xstart-x != int(xstart-x):
            print(P)
            print(P_transformed)
            print(x, xstart)

        assert xstart - x == int(xstart - x)
        assert ystart - y == int(ystart - y)
        for i, j in P_transformed:
            kx = int(i - xstart + x)
            ky = int(j - ystart + y)
            assert (kx, ky) not in GRID
            GRID[kx, ky] = cc + 1  # Assign a unique integer ID to each piece
        x += u / 2
        y += v / 2
        path_vertices.append((x, y))
        
    assert (x, y) == (0, 0.5)  # returned to start position

    print("solution constructed, start verification", flush=True)
    passed, first_a, total_a, msg = verify_solution(
        GRID, n_max, len(Polyominoes[n_max]), Polyominoes[n_max],
        do_pieces=not args.no_piece_check, do_areas=not args.no_area_check,
        do_shoelace=not args.no_shoelace,
        expected_first=expected_first, expected_total=expected_total,
    )
    print(msg, flush=True)
    print(f"ENCLOSED: first={first_a} total={total_a}", flush=True)

    # Generate visual output
    save_fence_visualization(path_vertices, n_max)

    if passed:
        out = os.path.join("data", f"sol-{n_max}-{first_a}-{total_a}.csv")
        if args.no_csv:
            print(f"Verification PASSED — first={first_a} total={total_a} (CSV skipped: --no-csv)")
        elif write_csv(out, GRID, len(Polyominoes[n_max])):
            print(f"Verification PASSED — {out}")
        else:
            print(
                f"Verification PASSED — first={first_a} total={total_a} (grid too large for CSV; path would be {out})"
            )
    else:
        print("Verification FAILED")


## Incrementally generate the set of n-ominoes

LETTERS = {
    5: "VPULTFZYWNIX",
    4: "oltni",  # assigned a posteriori after looking at the output
    2: "A ",
    3: "LI",
}
Defaultletters = (
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz" + "0123456789"
)  # arbitrary
assert len(Defaultletters) == 2 * 26 + 10


def Letters(n, i):
    if n in LETTERS:
        return LETTERS[n][i]
    return Defaultletters[i % len(Defaultletters)]

Vectors = []
for n in range(2, n_max + 1):
    print(f"Start generating the {n}-ominoes.", flush=True)
    Vectors = []
    New = set() # find the set of polyominoes of size n
    parents = Polyominoes[n-1]
    for pi, P in enumerate(parents):
        print_progress(pi + 1, len(parents), prefix=f"Generating {n}-ominoes:",
                       suffix=f"({pi+1}/{len(parents)} parents, {len(New)} found)")
        for x,y in P:
            for q in ((x+1,y),(x,y+1),(x-1,y),(x,y-1)):
                if q not in P:
                    P2 = P + (q,)
                    # normalize under rotations and reflections.
                    Q = min(tuple(sorted(normalize(P4))) for
                            P3 in (P2,rotateL(P2),rotateR(P2),rotate2(P2))
                            for P4 in (P3,reflect(P3)))
                    New.add(Q)
    
    New = sorted(tuple(sorted(reflect(Q))) for Q in New) # (reflected is more wide than high.)
    Polyominoes.append(New)
    print(f"\n{n=}. These are the {len(New)} {n}-ominoes:", flush=True)
    num_L_piece = None
    total_inner_holes = 0
    for i, P in enumerate(New):  # compute the diametrical pairs
        if True:
            print_progress(i + 1, len(New), prefix=f"Processing {n}-ominoes:", suffix=f"({i+1}/{len(New)})")
            
        name = Letters(n, i)
        diam = diameter(P)
        if n==5 and name=="U": # hack
            print(name, i, "hack U-pentomino ******")
            print(diam)
            assert diam ==(10, 2, [((0, 0), (2, 1)), ((0, 1), (2, 0))])
            diam = (9,1,[((0,1),(2,1))])
        if False and n==6 and name=="L": # hack
            # L       # not chosen
            # LLLLL
            print(name, i, "hack L-hexomino ******")
            assert diam == (26, 1, [((0, 1), (4, 0))])
            diam =         (25, 1, [((0, 0), (4, 0))])
        if False and n==6 and name=="K": # hack
            # K  K
            # KKKK
            print(name, i, "hack U-hexomino ******")
            assert diam ==(17, 2, [((0, 0), (3, 1)), ((0, 1), (3, 0))])
            diam =        (16, 1, [((0, 1), (3, 1))])
        if n==6 and name=="S": # hack 2
            #  S  
            #  S  
            # SSSS
            print(name, i, "hack T-hexomino ******")
            assert diam == (16, 1, [((0, 0), (3, 0))])
            diam =         (13, 1, [((1, 2), (3, 0))])
        if n==7 and name=="e" and i+1==31:
            # e     
            # eeeeee
            print(name, i, "hack L-heptomino ******")
            assert diam == (37, 1, [((0, 1), (5, 0))]) 
            diam = (36, 1, [((0, 0), (5, 0))])
        if False and n==7 and name=="d" and i+1==30:
            # d   d
            # ddddd
            print(name, i, "hack U-heptomino ******")
            assert diam == (26, 2, [((0, 0), (4, 1)), ((0, 1), (4, 0))])
            diam = (25, 1, [((0, 1), (4, 1))])
        if n==8 and name=="l" and i+1==100:
            # l     
            # lllllll
            print(name, i, "hack L-octomino ******")
            assert diam == (50, 1, [((0, 1), (6, 0))])
            diam = (49, 1, [((0, 0), (6, 0))])
        if False and n==8 and name=="k" and i+1==99:
            # k    k
            # kkkkkk
            print(name, i, "hack U-octomino ******")
            assert diam == (37, 2, [((0, 0), (5, 1)), ((0, 1), (5, 0))])
            diam = (36, 1, [((0, 1), (5, 1))])

        if n > 2:
            best_d, best_ori, innerholes = bestposition(*diam)
            total_inner_holes += innerholes
            if innerholes:
                vprint(innerholes,P)
        if n >= 4:
            vec_list = [
                frozenset([abs(x1 - x2), abs(y1 - y2)])
                for ((x1, y1), (x2, y2)), o in best_ori
            ]
            if len(set(vec_list))>1:# ... would be interesting
                vprint(P)
                vprint(best_ori)
                #exit()
            chosen_vec = vec_list[0]
            Vectors.append((sorted(list(chosen_vec)), best_d, best_ori[0], name, P))

        L_piece = ((0,0),(0,1))+tuple((i,0) for i in range(1,n_max-1)) 
        if n_max>8 and P == L_piece:
            print(diam)
            print(((n_max-1)**2+1, 1, [((0,1),(n_max-2,0))]))
            assert diam == ((n_max-1)**2+1, 1, [((0,1),(n_max-2,0))])
            num_L_piece = i
            diam2 = ((n_max-1)**2, 1, [((0,0),(n_max-2,0))])
            best_d2,best_ori2, innerholes = bestposition(*diam2)
            vec2 = set([frozenset([abs(x1-x2),abs(y1-y2)])
                                 for ((x1,y1),(x2,y2)),o in best_ori2])
            assert len(vec2)==1
            Vectors_alternate = (sorted(list(vec2)[0]),best_d2,best_ori2[0],name,P)
            
    print(f"Done with processing {n}-ominoes.", flush=True)

    if n < n_max:
        print("Proceed to next larger n", flush=True)
        continue
    if n>8: # for n<=7 has been hacked by hand
        assert num_L_piece is not None
        print(Vectors[num_L_piece])
        print(Vectors_alternate)

        total_parity_sum = sum(sum(v)for v,_,_,_,_ in Vectors if len(v)==2)
        print(total_parity_sum)
        if (total_parity_sum + len(Vectors)) %2:
            # +1 is added to all but 4 vectors
            Vectors[num_L_piece]= Vectors_alternate
            # modify usage of the L-piece to achieve correct parity

    print(f"Start sorting {len(Vectors)} vectors.", flush=True)
    Vectors.sort()
    deficit = sum(b for (v, b, o, l, P) in Vectors)
    print(f"Total deficit from all pieces = {deficit}", flush=True)

    Vecs = []
    if n==5:
        exceptions = 1 # number of pieces that require exceptional treatment
    if n==4:
        exceptions = 3
    if n>=6:
        exceptions = 0
    num_ex = 0
    gain = 0 # gain from treating pieces as corners instead of "straight" pieces

    # identify the four corner pieces
    corners = []
    longest_diagonal = max(xy[0] for xy, _, _, _, _ in Vectors if len(xy) == 1)
    expected_corner_length = n - 1 + (n % 2)
    print(n, expected_corner_length, longest_diagonal)
    assert expected_corner_length == 2 * longest_diagonal + 1
    for xy, d, o, l, P in Vectors:
        if n == 6 and l == "S" and xy == [2]:  # hack 2
            xy = [2, 2]
        if len(xy) == 1 and (
            n < 7 or (xy[0] == longest_diagonal and len(corners) < 4)
        ):
            x = xy[0]
            y = x
            if n == 7:
                assert x == 3
            if n == 5:
                assert x == 2
            if n == 4:
                assert x == 1
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
                assert x <= y
                y += 1
        Vecs.append((x / y, (x, y), l))
    print(f"{gain=} {num_ex=}")
    assert num_ex == 4
    assert exceptions == 0
    print("sort again", flush=True)
    Vecs.sort()
    print("done with sorting", flush=True)

    print(f"{corners=}")
    # the four "diagonal" pieces should cancel out in terms of even/odd
    total_parity_sum = sum(x+y for _,(x,y),_ in Vecs)
    print("total L1-distance must be even:", total_parity_sum)
    assert total_parity_sum%2 == 0 # otherwise solution is impossible.

    VecsR = [(y, x) for _, (x, y), _ in Vecs]  # reflect the vectors again
    B, SOL, CORNERS = best_polygon(VecsR)
    print(f"best solution value {B/8}. **** {B/8+deficit+gain=} ************",
          f"{total_inner_holes=}")
    calculated_area = B/8+deficit+gain
    assert calculated_area==int(calculated_area)
    calculated_area=int(calculated_area)
    print()
    if n<=7: print(f"{SOL=}; {CORNERS=}")
    print(f"{sum(x for x,y in SOL)=}; {sum(y for x,y in SOL)=}")
    construct_fence(SOL, CORNERS, Vectors, calculated_area, calculated_area + total_inner_holes)
