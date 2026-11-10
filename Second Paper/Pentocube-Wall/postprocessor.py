#!/usr/bin/env python3

import sys
import numpy as np
import re
import json
from pathlib import Path
from polycubetools import Coordinate
from polycubetools.utils import load_pentacubes


def debug_field(x, y, z, arr):
    val = arr[x, y, z]
    if val == 0 or val is None:
        return " . "
    s = str(val)
    if s.startswith("#"):
        s = s[1:]
    # Shorten to first 1-2 chars but try to fit in 3 chars
    # Center or align
    return f"{s[:3]:^3}"

def pretty_print(matrix, zero=" . ", sep=" "):
    if matrix is None:
        print("[]")
        return

    arr = np.array(matrix)

    if arr.ndim == 1:
        print(sep.join(debug_field(x,0,0, np.array([x for x in arr]).reshape((len(arr),1,1))).strip() for x in arr)) # minimal hack or just rewrite logic
        # Actually easier to just check 0/None
    elif arr.ndim == 2:
        for row in arr:
             print(sep.join(str(x) if (x != 0 and x is not None) else zero for x in row))
    elif arr.ndim == 3:
        width, height, depth = arr.shape

        for row in range(height):
            row_str = []
            for d in range(depth):
                # Convert row of current layer to string
                row_str.append(
                    " ".join(debug_field(x,row,d, arr) for x in range(width))
                )
            # Join layers with some spacing
            print("   ".join(row_str))
    else:
        for i in range(arr.shape[0]):
            pretty_print(arr[i], zero=zero, sep=sep)

def parse_solution(file_path):
    placements = []

    text = Path(file_path).read_text(encoding="utf-8")
    wall = set()
    if "=======\n" in text:
        instance, text = text.split("=======\n", 1)
        item_line = next(
            line for line in instance.splitlines()
            if line.strip() and not line.lstrip().startswith("|")
        )
        wall = {
            (int(token[0:2]), int(token[2:4]), int(token[4:6]))
            for token in item_line.split()
            if re.fullmatch(r"\d{6}", token)
        }
        text = text.split("====\n", 1)[0]

    for line in text.splitlines():
            line = line.strip()
            if not line:
                continue

            # Extract ID usually starting with #
            # Format example: #15 111 211 212 221 311 (1 of 63)
            # Or: 116 126 136 146 236 #10 (1 of 23)

            parts = line.split()
            cube_id = None
            coords = []

            for part in parts:
                if part.startswith('(') or part.endswith(')') or part == 'of':
                    continue

                # Check for coordinates
                is_coord = False
                if len(part) == 6 and part.isdigit():
                     is_coord = True
                     x = int(part[0:2])
                     y = int(part[2:4])
                     z = int(part[4:6])
                     coords.append(Coordinate(x, y, z))
                elif re.match(r'^\d{3}$', part):
                     is_coord = True
                     x = int(part[0])
                     y = int(part[1])
                     z = int(part[2])
                     coords.append(Coordinate(x, y, z))

                if not is_coord:
                    # It's an ID
                    # Only take the first one?
                    # If we already have one, keep the first one?
                    if cube_id is None:
                        cube_id = part

            if cube_id is not None and coords:
                placements.append({'id': cube_id, 'coords': coords})
            elif "=======" not in line and not line.endswith(":"):
                print(f"Warning: Could not parse line: {line}")

    incomplete = [placement for placement in placements if len(placement['coords']) != 5]
    if incomplete:
        if len(incomplete) != 1 or not wall:
            raise ValueError("Cannot restore incomplete pentacube without its DLX wall definition")

        placement = incomplete[0]
        piece_id = int(placement['id'].removeprefix("P"))
        covered = {
            (coord.x, coord.y, coord.z) for coord in placement['coords']
        }
        cube = next(cube for cube in load_pentacubes() if cube.id == piece_id)
        full_placements = set()
        for rotation in cube.rotations:
            shape = {(coord.x, coord.y, coord.z) for coord in rotation.coords}
            for shape_coord in shape:
                for covered_coord in covered:
                    offset = tuple(covered_coord[axis] - shape_coord[axis] for axis in range(3))
                    full_shape = frozenset(
                        tuple(coord[axis] + offset[axis] for axis in range(3))
                        for coord in shape
                    )
                    if full_shape.intersection(wall) == covered:
                        full_placements.add(full_shape)

        if len(full_placements) != 1:
            raise ValueError(
                f"Expected one full placement for P{piece_id}, found {len(full_placements)}"
            )
        placement['coords'] = [Coordinate(*coord) for coord in sorted(full_placements.pop())]

    return placements

def to_grid_json(placements):
    all_coords = [coord for placement in placements for coord in placement['coords']]
    if not all_coords:
        raise ValueError("No polycube coordinates found in solution")

    dimensions = [max(getattr(coord, axis) for coord in all_coords) + 3
                  for axis in ("x", "y", "z")]
    rows = [[int(placement['id'].removeprefix("P")), coord.x, coord.y, coord.z]
            for placement in placements for coord in placement['coords']]
    return [dimensions, *rows]

def visualize(file_path, output_path=None):
    placements = parse_solution(file_path)
    if not placements:
        print("No placements found.")
        return

    output_path = Path(output_path) if output_path else Path(file_path).with_suffix(".json")
    output_path.write_text(json.dumps(to_grid_json(placements), indent=2) + "\n", encoding="utf-8")

    # Determine bounds
    all_coords = [c for p in placements for c in p['coords']]
    if not all_coords:
        return

    max_x = max(c.x for c in all_coords)
    max_y = max(c.y for c in all_coords)
    max_z = max(c.z for c in all_coords)

    # Create board (1-indexed input effectively, so +2 for clean borders)
    # We'll map coordinates directly. If they are 1-based, index 0 will be empty.
    board = np.zeros((max_x + 2, max_y + 2, max_z + 2), dtype=object)

    for p in placements:
        cid = p['id']
        for c in p['coords']:
            # Check bounds just in case
            if c.x < board.shape[0] and c.y < board.shape[1] and c.z < board.shape[2]:
                board[c.x, c.y, c.z] = cid
            else:
                print(f"Coord out of bounds: {c}")

    # print(f"Board Size: {max_x}x{max_y}x{max_z}")
    # pretty_print(board)

if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "solution.dlx"
    output_path = sys.argv[2] if len(sys.argv) > 2 else None
    visualize(path, output_path)
