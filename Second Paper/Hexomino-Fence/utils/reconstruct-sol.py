"""
Reconstruct polyomino fence solutions from reference polygon solutions.

This script takes solutions in terms of reference vectors (35 for hexominoes,
12 for pentominoes) read from an ``all-reference-polygons-<n>t<area>.txt`` file,
and reconstructs the valid polyomino placements (fences) that produce those
vectors. Each reconstructed fence is verified (single enclosed region of the
target area, every free polyomino used exactly once, no overlaps, no duplicates)
and written to a log file as an ASCII grid.

Ensure that the open source ILP solver HiGHS is available (`pip install highspy`).
Examples:

    # Reconstruct everything in the input file (up to 5 fences per
    # solution), writing to reconstructed-solutions.log
    python3 utils/reconstruct-sol.py data/all-reference-polygons-6t1597.txt

    # Pentominoes (n = 5, area 128), only solution number 1
    python3 utils/reconstruct-sol.py data/all-reference-polygons-5t128.txt --start 1 --end 1

    # Reconstruct up to 50 concrete solutions from vector sets 10 to 20 each
    python3 utils/reconstruct-sol.py data/all-reference-polygons-6t1597.txt --start 10 --end 20 --limit 50
"""

import argparse
import ast
import os
import re
import sys
from collections import defaultdict

import highspy

# ============================================================================
# SECTION 1: Polyomino definitions (same as step0)
# ============================================================================

DIRECTIONS = (0, -1), (1, 0), (0, 1), (-1, 0)  # 0=South, 1=East, 2=North, 3=West

def rotateL(P):
    return [(-y, x) for x, y in P]

def rotateR(P):
    return [(y, -x) for x, y in P]

def rotate2(P):
    return [(-x, -y) for x, y in P]

def rotate0(P):
    return [(x, y) for x, y in P]

def rotate(P, dir):
    """Rotate counterclockwise by 'dir' multiples of 90°"""
    return (rotate0, rotateL, rotate2, rotateR)[dir % 4](P)

def reflect(P):
    return [(y, x) for x, y in P]

def normalize(P):
    xmin = min(x for x, y in P)
    ymin = min(y for x, y in P)
    return [(x - xmin, y - ymin) for x, y in P]

def canonicalize(cells):
    """
    Returns the canonical form of a polyomino (free polyomino).
    The canonical form is the lexicographically smallest among all 8 symmetries,
    normalized to (0,0).
    """
    def bounds(shape):
        xs = [x for x, y in shape]
        ys = [y for x, y in shape]
        return min(xs), min(ys)

    def norm_shape(shape):
        min_x, min_y = bounds(shape)
        return tuple(sorted([(x - min_x, y - min_y) for x, y in shape]))

    def rotate90(shape):
        return [(y, -x) for x, y in shape]

    def reflect_canon(shape):
        return [(-x, y) for x, y in shape]

    variations = []
    current = list(cells)
    # Rotate 4 times
    for _ in range(4):
        variations.append(norm_shape(current))
        variations.append(norm_shape(reflect_canon(current)))
        current = rotate90(current)

    return min(variations)

Letters6 = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijk"

def generate_polyominoes(n_max):
    """Generate all free polyominoes with exactly n_max cells."""
    Polyominoes = [[], [((0, 0),)]]
    for n in range(2, n_max + 1):
        New = set()
        for P in Polyominoes[n - 1]:
            for x, y in P:
                for q in ((x + 1, y), (x, y + 1), (x - 1, y), (x, y - 1)):
                    if q not in P:
                        P2 = P + (q,)
                        Q = min(tuple(sorted(normalize(P4))) for
                                P3 in (P2, rotateL(P2), rotateR(P2), rotate2(P2))
                                for P4 in (P3, reflect(P3)))
                        New.add(Q)
        New = sorted(tuple(sorted(reflect(Q))) for Q in New)
        Polyominoes.append(New)
    return Polyominoes[n_max]


def parse_input_filename(path):
    """
    Extract (n_max, target_area) from the end of an all-reference-polygons-<n>t<area>.txt
    filename. n_max is 5 for pentominoes and 6 for hexominoes.
    """
    name = os.path.basename(path)
    m = re.search(r"(\d+)t(\d+)\.[^.]*$", name)
    if m is None:
        raise ValueError(f"Cannot extract '<n>t<area>' from the end of filename '{name}'")
    n_max, target_area = int(m.group(1)), int(m.group(2))
    if n_max not in (5, 6):
        raise ValueError(f"Only pentominoes (n=5) and hexominoes (n=6) are supported, got n={n_max}")
    return n_max, target_area


# ============================================================================
# SECTION 2: Compute hexomino placements for fence reconstruction
# ============================================================================

DIR_TO_INDEX = {
    (0, -1): 0,  # South
    (1, 0): 1,   # East
    (0, 1): 2,   # North
    (-1, 0): 3,  # West
}

INDEX_TO_DIR = {v: k for k, v in DIR_TO_INDEX.items()}

def get_orientations(P):
    """Get all unique orientations of a polyomino (up to rotation and reflection)"""
    base = [tuple(sorted(normalize(rotate(list(P), r)))) for r in range(4)]
    flipped = [tuple(sorted(normalize(rotate(reflect(list(P)), r)))) for r in range(4)]
    all_orients = base + flipped
    seen = set()
    unique = []
    for o in all_orients:
        if o not in seen:
            seen.add(o)
            unique.append(list(o))
    return unique


def compute_boundary_path(tiles):
    """Compute the clockwise boundary path of a set of tiles"""
    boundary_successor = {}
    edge_to_tile = {}
    
    for (x, y) in tiles:
        for (dx, dy) in DIRECTIONS:
            xx, yy = x + dx, y + dy
            if (xx, yy) not in tiles:
                perpx2, perpy2 = -dy, dx
                midx2, midy2 = x + xx, y + yy
                start_pt = (midx2 + perpx2, midy2 + perpy2)
                end_pt = (midx2 - perpx2, midy2 - perpy2)
                boundary_successor[start_pt] = end_pt
                edge_to_tile[(start_pt, end_pt)] = ((x, y), (xx, yy))
    
    start_pt = min(boundary_successor.keys())
    current = start_pt
    path = [start_pt]
    current = boundary_successor[current]
    while current != start_pt:
        path.append(current)
        current = boundary_successor[current]
    
    return path, edge_to_tile


class Possibility:
    """Represents a possible way to place a hexomino in the fence"""
    def __init__(self, new_orient, area_change_8x, pos_diff, hexo_name, shifted_tiles, next_tile_offset, end_tile_offset, ref_vector, start_orient, hexo_idx, rotation, rotated_vec, penalty):
        self.new_orient = new_orient
        self.area_change_8x = area_change_8x
        self.pos_diff = pos_diff
        self.hexo_name = hexo_name
        self.shifted_tiles = shifted_tiles
        self.next_tile_offset = next_tile_offset
        self.end_tile_offset = end_tile_offset
        self.ref_vector = ref_vector  # World-coordinate reference vector for matching
        self.start_orient = start_orient
        self.hexo_idx = hexo_idx
        self.rotation = rotation  # Rotation value for penalty lookup
        self.rotated_vec = rotated_vec  # Rotated reference vector for penalty lookup
        self.penalty = penalty  # Actual penalty for THIS specific placement


def compute_possibilities(pieces):
    """
    Compute all possibilities for placing each hexomino given a starting orientation.
    Returns dict: (hexo_idx, start_orient) -> list of Possibility
    """
    possibilities = defaultdict(list)
    
    for hexo_idx, P in enumerate(pieces):
        orients = get_orientations(P)
        
        for orient_idx, tiles in enumerate(orients):
            path, edge_to_tile = compute_boundary_path(tiles)
            doubled_path = path + path
            len_path = len(path)
            
            # Compute directions for each edge (step0 style)
            directions = []
            for i in range(len_path):
                p1, p2 = doubled_path[i], doubled_path[i + 1]
                dx, dy = p2[0] - p1[0], p2[1] - p1[1]
                # Our path is in doubled coords, so divide by 2 to get unit direction
                directions.append(DIRECTIONS.index((dx // 2, dy // 2)))
            directions = directions * 2  # Double for wrap-around
            
            # Compute absolute_direction (cumulative direction changes)
            absolute_direction = [0]
            dir_val = 0
            for i in range(len(directions) - 1):
                dir_change = (directions[i + 1] - directions[i]) % 4
                dir_val += {3: -1, 0: 0, 1: 1}.get(dir_change, 0)
                absolute_direction.append(dir_val)
            
            # Compute area_under for penalty calculation
            # Path is in doubled coordinates, so area_under is also doubled (then divided by 4 for half-int coords)
            # Using step0's formula: area_under = (x2-x1)*(y1+y2) where (x1,y1), (x2,y2) are half-integer coords
            # Our path is doubled, so need to scale: path/2 for half-int coords
            area_under = []
            for i in range(len_path):
                x1, y1 = doubled_path[i][0] / 2, doubled_path[i][1] / 2  # Convert to half-int coords
                x2, y2 = doubled_path[i + 1][0] / 2, doubled_path[i + 1][1] / 2
                area_under.append((x2 - x1) * (y1 + y2))
            area_under = area_under * 2  # Double for wrap-around
            
            sum_area_under = [0]
            s = 0
            for a in area_under:
                s += a
                sum_area_under.append(s)
            
            for start_index in range(len_path):
                for end_index in range(start_index + 1, start_index + len_path):
                    p1 = doubled_path[start_index]
                    p2 = doubled_path[start_index + 1]
                    diff_start = (p2[0] - p1[0], p2[1] - p1[1])
                    
                    # Convert to direction index (doubled coordinates)
                    diff_start_norm = (diff_start[0] // 2, diff_start[1] // 2)
                    if diff_start_norm not in DIR_TO_INDEX:
                        continue
                    start_orient = DIR_TO_INDEX[diff_start_norm]
                    
                    p3 = doubled_path[end_index]
                    p4 = doubled_path[end_index + 1]
                    diff_end = (p3[0] - p4[0], p3[1] - p4[1])
                    diff_end_norm = (diff_end[0] // 2, diff_end[1] // 2)
                    if diff_end_norm not in DIR_TO_INDEX:
                        continue
                    new_orient = DIR_TO_INDEX[diff_end_norm]
                    
                    # Compute rotation (step0 style)
                    rotation = absolute_direction[end_index] - absolute_direction[start_index] + 2
                    
                    origin = doubled_path[start_index + 1]
                    inner_path = [
                        (pt[0] - origin[0], pt[1] - origin[1])
                        for pt in doubled_path[(start_index + 1):(end_index + 1)]
                    ]
                    
                    area_change_8x = 0
                    for e in range(len(inner_path) - 1):
                        area_change_8x += (
                            inner_path[e][0] * inner_path[e + 1][1]
                            - inner_path[e + 1][0] * inner_path[e][1]
                        )
                    
                    start_tile_key = (doubled_path[start_index], doubled_path[start_index + 1])
                    start_tile, _ = edge_to_tile[start_tile_key]
                    shifted_tiles = [(t[0] - start_tile[0], t[1] - start_tile[1]) for t in tiles]
                    
                    end_tile_key = (doubled_path[end_index], doubled_path[end_index + 1])
                    end_tile, next_tile = edge_to_tile[end_tile_key]
                    shifted_end_tile = (end_tile[0] - start_tile[0], end_tile[1] - start_tile[1])
                    shifted_next_tile = (next_tile[0] - start_tile[0], next_tile[1] - start_tile[1])
                    
                    # Compute reference vector (world coordinates - NOT rotated)
                    # This is what matches solution vectors
                    # Edge midpoints in regular coords (path uses doubled coords, so / 4)
                    start_mid = ((p1[0] + p2[0]) / 4, (p1[1] + p2[1]) / 4)
                    end_mid = ((p3[0] + p4[0]) / 4, (p3[1] + p4[1]) / 4)
                    ref_vec = (end_mid[0] - start_mid[0], end_mid[1] - start_mid[1])
                    
                    # Compute rotated vector for penalty lookup (step0 style)
                    dirS = DIRECTIONS.index(((p1[0] - p2[0]) // 2, (p1[1] - p2[1]) // 2))
                    rotated_vec = rotate([ref_vec], -dirS)[0]
                    
                    # Compute penalty directly (step0 formula)
                    # start_mid/end_mid are in half-integer coords
                    along_boundary = (area_under[start_index] + area_under[end_index]) / 2 + \
                                     sum_area_under[end_index] - sum_area_under[start_index + 1]
                    a_under_ref_vector = (end_mid[0] - start_mid[0]) * (end_mid[1] + start_mid[1])
                    penalty = (along_boundary - a_under_ref_vector) / 2
                    
                    poss = Possibility(
                        new_orient,
                        area_change_8x,
                        inner_path[-1] if inner_path else (0, 0),
                        Letters6[hexo_idx],
                        shifted_tiles,
                        shifted_next_tile,
                        shifted_end_tile,
                        ref_vec,
                        start_orient,
                        hexo_idx,
                        rotation,
                        rotated_vec,
                        penalty
                    )
                    possibilities[(hexo_idx, start_orient)].append(poss)
    
    return possibilities


# ============================================================================
# SECTION 3: Build mapping from (signed_ref_vec, start_orient) -> possibilities
# ============================================================================

def build_vector_mapping(possibilities):
    """
    Build mapping from (signed_ref_vec, start_orient) -> list of possibilities.
    This matches the verify script's approach.
    """
    mapping = defaultdict(list)
    
    for (hexo_idx, start_orient), poss_list in possibilities.items():
        for poss in poss_list:
            key = (poss.ref_vector, start_orient)
            mapping[key].append(poss)
    
    return mapping


# ============================================================================
# SECTION 4: MIP solver using HiGHS with orientation chaining
# ============================================================================

PENALTY_TOLERANCE = 0.001


def solve_assignment_with_orientation(solution_vectors, possibilities, target_penalty, num_pieces, lim_sol=10):
    """
    Find valid hexomino assignments using MIP with orientation chaining.
    
    1. For each position, find candidates matching the signed vector and orientation
    2. Chain orientations: new_orient of pos i must match start_orient of pos i+1
    3. Each hexomino used exactly once
    4. Penalty constraint (using penalty stored in each Possibility)
    
    Args:
        solution_vectors: List of 35 signed reference vectors from the solution
        possibilities: Precomputed possibilities dict
        target_penalty: Target penalty value (e.g., -66.0)
        num_pieces: Number of polyominoes (12 pentominoes or 35 hexominoes)
        lim_sol: Maximum number of solutions to find
        
    Returns:
        List of solutions, each solution is a list of Possibility objects
    """
    n = len(solution_vectors)
    assert n == num_pieces, f"Expected {num_pieces} vectors, got {n}"
    
    # Build vector mapping
    vec_mapping = build_vector_mapping(possibilities)
    
    # Build candidate lists for each position (considering all possible start orientations)
    # position_candidates[pos] = list of (start_orient, Possibility)
    position_candidates = []
    for i, vec in enumerate(solution_vectors):
        candidates = []
        for orient in range(4):  # 0, 1, 2, 3 for S, E, N, W
            key = (vec, orient)
            if key in vec_mapping:
                for poss in vec_mapping[key]:
                    candidates.append((orient, poss))
        
        if not candidates:
            print(f"  No candidates for position {i} with vector {vec}")
            return []
        position_candidates.append(candidates)
    
    # Map (pos, cand_idx) -> variable info
    var_info = []  # (pos, cand_idx, start_orient, poss)
    for pos in range(n):
        for cand_idx, (start_orient, poss) in enumerate(position_candidates[pos]):
            var_info.append((pos, cand_idx, start_orient, poss))
    
    num_vars = len(var_info)
    
    # Build hexomino to variable mapping for uniqueness constraint
    hexo_to_vars = defaultdict(list)
    for var_idx, (_, _, _, poss) in enumerate(var_info):
        hexo_to_vars[poss.hexo_idx].append(var_idx)
    
    # Create MIP model
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    
    # Create candidate variables
    x = [h.addBinary() for _ in range(num_vars)]
    
    # Create orientation indicator variables
    orient_var = {}
    for pos in range(n):
        for orient in range(4):
            orient_var[(pos, orient)] = h.addBinary()
    
    # Constraint 1: Each position has exactly one candidate selected
    for pos in range(n):
        pos_var_indices = [var_idx for var_idx, (p, _, _, _) in enumerate(var_info) if p == pos]
        if pos_var_indices:
            h.addConstr(sum(x[i] for i in pos_var_indices) == 1)
    
    # Constraint 2: Each hexomino used at most once
    for hexo_idx in range(num_pieces):
        if hexo_idx in hexo_to_vars and hexo_to_vars[hexo_idx]:
            h.addConstr(sum(x[i] for i in hexo_to_vars[hexo_idx]) <= 1)
    
    # Constraint 3a: Link candidate selection to outgoing orientation
    for var_idx, (pos, cand_idx, start_orient, poss) in enumerate(var_info):
        new_orient = poss.new_orient
        h.addConstr(x[var_idx] <= orient_var[(pos, new_orient)])
    
    # Constraint 3b: Exactly one outgoing orientation per position
    for pos in range(n):
        h.addConstr(sum(orient_var[(pos, o)] for o in range(4)) == 1)
    
    # Constraint 3c: Chain orientations - incoming must match outgoing of previous position
    for var_idx, (pos, cand_idx, start_orient, poss) in enumerate(var_info):
        if pos > 0:
            h.addConstr(x[var_idx] <= orient_var[(pos - 1, start_orient)])
    
    # Constraint 3d: Initial orientation - fence starts facing North
    initial_orient = 2  # North
    valid_first_vars = []
    for var_idx, (pos, cand_idx, start_orient, poss) in enumerate(var_info):
        if pos == 0 and start_orient == initial_orient:
            valid_first_vars.append(var_idx)
    
    if not valid_first_vars:
        print(f"  No valid initial candidates with orientation {initial_orient}")
        return []
    
    h.addConstr(sum(x[i] for i in valid_first_vars) >= 1)
    
    # Constraint 3e: Final orientation must equal initial orientation (closed loop)
    h.addConstr(orient_var[(n - 1, initial_orient)] == 1)
    
    # Constraint 4: Penalty sum must equal target
    penalty_coeffs = []
    for var_idx, (pos, cand_idx, start_orient, poss) in enumerate(var_info):
        penalty_coeffs.append(poss.penalty)
    
    penalty_expr = sum(penalty_coeffs[i] * x[i] for i in range(num_vars))  
    h.addConstr(penalty_expr >= target_penalty - PENALTY_TOLERANCE)
    h.addConstr(penalty_expr <= target_penalty + PENALTY_TOLERANCE)
    
    solutions = []
    
    # Repeatedly solve and add exclusion constraint
    while len(solutions) < lim_sol:
        # Show progress
        print(f"\r  Finding solutions: ({len(solutions)}/{lim_sol})", end='', flush=True)
        
        status = h.run()
        
        if h.getModelStatus() != highspy.HighsModelStatus.kOptimal:
            break
        
        # Extract solution
        sol_vals = h.getSolution().col_value
        
        path = [None] * n
        selected_vars = []
        total_penalty = 0
        
        for var_idx in range(num_vars):
            if sol_vals[var_idx] > 0.5:
                pos, cand_idx, start_orient, poss = var_info[var_idx]
                path[pos] = poss
                selected_vars.append(var_idx)
                total_penalty += penalty_coeffs[var_idx]
        
        if any(p is None for p in path):
            print(f"\r  Warning: incomplete solution" + " " * 20)
            break
        
        solutions.append((path, total_penalty))
        
        # Add exclusion constraint
        h.addConstr(sum(x[i] for i in selected_vars) <= n - 1)
    
    # Clear progress line
    print(f"\r  Found {len(solutions)} valid assignments" + " " * 20)
    
    return solutions


# ============================================================================
# SECTION 5: Fence reconstruction from Possibility path
# ============================================================================

def process_possibility_path(path):
    """
    Convert a path of Possibility objects to fence_tiles dict.
    The path is already in order with correct orientation chaining.
    
    Args:
        path: List of Possibility objects
        
    Returns:
        Dict mapping (x, y) -> hexomino letter, or None if overlap detected
    """
    fence_tiles = {}
    curr = (0, 0)
    
    for poss in path:
        for (x, y) in poss.shifted_tiles:
            pos = (x + curr[0], y + curr[1])
            if pos in fence_tiles:
                return None  # Overlap detected
            fence_tiles[pos] = poss.hexo_name
        curr = (curr[0] + poss.next_tile_offset[0], curr[1] + poss.next_tile_offset[1])
    
    return fence_tiles


def piece_contribution(poss):
    """Endpoint-inclusive connector-tile L1 contribution for one piece placement."""
    return abs(poss.end_tile_offset[0]) + abs(poss.end_tile_offset[1]) + 1


def calculate_circumference(path):
    """Calculate circumference as sum of endpoint-inclusive connector-tile L1 lengths."""
    contributions = [piece_contribution(poss) for poss in path]
    return sum(contributions), contributions


def calculate_perimeter_upper_bound(possibilities, num_pieces):
    """
    Compute per-piece maximum contribution over all precomputed placements.

    Returns:
        Tuple (upper_bound_sum, per_piece_max) where per_piece_max is a list
        indexed by hexomino index.
    """
    per_piece_max = [None] * num_pieces

    for (hexo_idx, _start_orient), poss_list in possibilities.items():
        if not poss_list:
            continue
        local_max = max(piece_contribution(poss) for poss in poss_list)
        current = per_piece_max[hexo_idx]
        if current is None or local_max > current:
            per_piece_max[hexo_idx] = local_max

    missing = [idx for idx, val in enumerate(per_piece_max) if val is None]
    if missing:
        raise ValueError(f"Missing possibilities for hexomino indices: {missing}")

    return sum(per_piece_max), per_piece_max


# ============================================================================
# SECTION 6: Area calculation
# ============================================================================

def calculate_enclosed_areas(fence_tiles):
    """Calculate enclosed areas using 8-connectivity flood fill"""
    if not fence_tiles:
        return []
    
    xs = [k[0] for k in fence_tiles.keys()]
    ys = [k[1] for k in fence_tiles.keys()]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    
    areas = []
    visited = set()
    # 8-connectivity (including diagonals)
    dirs8 = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
    
    for x in range(min_x, max_x + 1):
        for y in range(min_y, max_y + 1):
            if (x, y) not in fence_tiles and (x, y) not in visited:
                # BFS for connected component
                visited.add((x, y))
                queue = [(x, y)]
                size = 1
                touches_border = False
                
                head = 0
                while head < len(queue):
                    cx, cy = queue[head]
                    head += 1
                    
                    if cx == min_x or cx == max_x or cy == min_y or cy == max_y:
                        touches_border = True
                    
                    for dx, dy in dirs8:
                        nx, ny = cx + dx, cy + dy
                        if min_x <= nx <= max_x and min_y <= ny <= max_y:
                            if (nx, ny) not in fence_tiles and (nx, ny) not in visited:
                                visited.add((nx, ny))
                                queue.append((nx, ny))
                                size += 1
                
                if not touches_border:
                    areas.append(size)
    
    areas.sort(reverse=True)
    return areas


# ============================================================================
# SECTION 7: Parse solution file
# ============================================================================

def parse_solution_line(line):
    """
    Parse one line of an all-reference-polygons-<n>t<area>.txt file.

    The line is a Python tuple literal of the form
        (sol_num, penalty, vectors)

    Returns (sol_num, target_penalty, vectors), or (None, None, None) for blank,
    comment or unparsable lines.

    target_penalty is the value the MIP must reach, i.e. -penalty from the file
    (the file and the Possibility penalties use opposite sign conventions).
    """
    line = line.strip()
    if not line or line.startswith("#"):
        return None, None, None

    # Wrapping in one more pair of parentheses makes this work whether or not the
    # line is already parenthesised, and tolerates a trailing comma.
    try:
        record = ast.literal_eval(f"({line.rstrip(',')})")
    except (ValueError, SyntaxError):
        return None, None, None

    if not isinstance(record, tuple) or len(record) != 3:
        return None, None, None

    try:
        sol_num = int(record[0])
        file_penalty = float(record[1])
        vectors = tuple(tuple(v) for v in record[2])
    except (TypeError, ValueError):
        return None, None, None

    return sol_num, -file_penalty, vectors  # different definitions of penalties


# ============================================================================
# SECTION 8: Main program
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description="Reconstruct polyomino fence solutions")
    parser.add_argument("input", help="Input file with reference polygons, named all-reference-polygons-<n>t<area>.txt "
                                      "(n=5: pentominoes, n=6: hexominoes)")
    parser.add_argument("--output", default="reconstructed-solutions.log", help="Output log file")
    parser.add_argument("--limit", type=int, default=5, help="Max solutions per vector solution")
    parser.add_argument("--start", type=int, default=1, help="Starting solution number (inclusive)")
    parser.add_argument("--end", type=int, default=None, help="Ending solution number (inclusive)")
    args = parser.parse_args()
    
    try:
        n_max, target_area = parse_input_filename(args.input)
    except ValueError as e:
        parser.error(str(e))
    
    pieces = generate_polyominoes(n_max)
    num_pieces = len(pieces)
    kind = "pentominoes" if n_max == 5 else "hexominoes"
    print(f"Generated {num_pieces} {kind}; target area {target_area}")
    
    print("Computing placement possibilities...")
    possibilities = compute_possibilities(pieces)
    perimeter_upper_bound, _ = calculate_perimeter_upper_bound(possibilities, num_pieces)
    print(f"Perimeter upper bound (sum of per-piece maxima): {perimeter_upper_bound}")
    
    print(f"Reading solutions from {args.input}...")
    solutions_to_process = []
    
    with open(args.input, "r") as f:
        for line in f:
            sol_num, penalty, vectors_sol = parse_solution_line(line)
            if sol_num is not None:
                if args.start <= sol_num and (args.end is None or sol_num <= args.end):
                    solutions_to_process.append((sol_num, penalty, vectors_sol))
    
    print(f"Found {len(solutions_to_process)} solutions to process")
    
    total_reconstructed = 0
    
    # Get canonical forms of all polyominoes for verification
    all_hexominoes_canonical = set()
    for P in pieces:
        all_hexominoes_canonical.add(canonicalize(P))
    
    # Track unique grids to detect duplicates
    seen_grids = set()
    
    with open(args.output, "w") as log:
        for sol_num, penalty, vectors_sol in solutions_to_process:
            print(f"\nProcessing solution {sol_num}...")
            log.write(f"\n{'='*60}\n")
            log.write(f"Solution {sol_num} from vector file\n")
            log.write(f"{'='*60}\n\n")
            log.flush()
            
            if len(vectors_sol) != num_pieces:
                print(f"  Skipping: expected {num_pieces} vectors, got {len(vectors_sol)}")
                log.write(f"SKIPPED: expected {num_pieces} vectors, got {len(vectors_sol)}\n")
                log.flush()
                continue
            
            # Find hexomino assignments with orientation chaining
            assignments = solve_assignment_with_orientation(list(vectors_sol), possibilities, penalty, num_pieces, lim_sol=args.limit)
            
            if not assignments:
                log.write("No valid hexomino assignments found (penalty/orientation constraints not satisfied)\n")
                log.flush()
                continue
            
            # Each assignment is already a complete fence path (list of Possibility objects)
            for assign_idx, (path, total_penalty) in enumerate(assignments):
                # Convert path to fence tiles
                fence_tiles = process_possibility_path(path)
                circumference, piece_contributions = calculate_circumference(path)
                
                if fence_tiles is None:
                    print(f"  ERROR: Overlapping tiles in assignment {assign_idx + 1}")
                    sys.exit(1)
                
                # Check 1: Verify enclosed area is exactly the target area
                enclosed = calculate_enclosed_areas(fence_tiles)
                if len(enclosed) != 1 or enclosed[0] != target_area:
                    print(f"  ERROR: Expected single enclosed area of {target_area}, got {enclosed}")
                    sys.exit(1)
                
                # Check 2: Verify that all polyominoes are distinct (not just rotations/reflections)
                # Group tiles by hexomino letter
                hexomino_tiles = {}
                for (x, y), letter in fence_tiles.items():
                    if letter not in hexomino_tiles:
                        hexomino_tiles[letter] = []
                    hexomino_tiles[letter].append((x, y))
                
                # Canonicalize each hexomino found
                found_canonicals = set()
                for letter, tiles in hexomino_tiles.items():
                    canonical = canonicalize(tiles)
                    found_canonicals.add(canonical)
                
                # Verify: exactly num_pieces distinct polyominoes matching the canonical set
                if found_canonicals != all_hexominoes_canonical:
                    print(f"  ERROR: Not all {num_pieces} unique polyominoes present")
                    sys.exit(1)
                
                # Build grid representation
                xs = [k[0] for k in fence_tiles.keys()]
                ys = [k[1] for k in fence_tiles.keys()]
                min_x, max_x = min(xs), max(xs)
                min_y, max_y = min(ys), max(ys)
                
                grid_rows = []
                for y in range(max_y, min_y - 1, -1):
                    row = ""
                    for x in range(min_x, max_x + 1):
                        if (x, y) in fence_tiles:
                            row += fence_tiles[(x, y)]
                        else:
                            row += "."
                    grid_rows.append(row)
                
                # Check for duplicate grid
                grid_tuple = tuple(grid_rows)
                if grid_tuple in seen_grids:
                    print(f"  ERROR: Duplicate grid detected in assignment {assign_idx + 1}")
                    sys.exit(1)
                seen_grids.add(grid_tuple)
                
                # Valid solution - log it
                total_reconstructed += 1
                log.write(f"Assignment {assign_idx + 1}:\n")
                log.write(f"Circumference: {circumference}\n")
                # Count how often each contribution value appears
                contrib_counts = {}
                for v in piece_contributions:
                    contrib_counts[v] = contrib_counts.get(v, 0) + 1
                #log.write(f"Piece contribution counts: {contrib_counts}\n\n")
                
                # Write grid
                for row in grid_rows:
                    log.write(row + "\n")
                
                log.write("\n")
                log.flush()
    
    print(f"\nTotal reconstructed solutions: {total_reconstructed}")
    print(f"Results written to {args.output}")


if __name__ == "__main__":
    main()