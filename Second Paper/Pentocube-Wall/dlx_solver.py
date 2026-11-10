import subprocess
import sys
from pathlib import Path

from polycubetools.utils import load_pentacubes


INNER_CELLS = {
    (x, y, z)
    for x in range(1, 5)
    for y in range(1, 5)
    for z in range(1, 4)
}
EXTENSIONS = (
    [[(x, y, 4) for x in (x0, x0 + 1) for y in (y0, y0 + 1)]
     for x0, y0 in ((1, 1), (1, 2), (2, 2))]
    + [[(5, y, z) for y in (y0, y0 + 1) for z in (1, 2)]
       for y0 in (1, 2)]
)

PROJECT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = PROJECT_DIR / "output"
SOLVER_PATH = PROJECT_DIR / "ssxcc"
POSTPROCESSOR_PATH = PROJECT_DIR / "postprocessor.py"
DLX_PATH = OUTPUT_DIR / "solfound.dlx"
JSON_PATH = OUTPUT_DIR / "solfound.json"


def find_coordinate_bounds(points):
    min_coords = tuple(min(point[axis] for point in points) for axis in range(3))
    max_coords = tuple(max(point[axis] for point in points) for axis in range(3))
    return min_coords, max_coords


def build_wall(inner_cells):
    wall = {
        (x + dx, y + dy, z + dz)
        for x, y, z in inner_cells
        for dx in range(-1, 2)
        for dy in range(-1, 2)
        for dz in range(-1, 2)
    }
    wall.difference_update(inner_cells)
    return wall


def build_dlx_input(wall, pieces, special_piece):
    (min_x, min_y, min_z), (max_x, max_y, max_z) = find_coordinate_bounds(wall)
    wall_widths = (max_x - min_x, max_y - min_y, max_z - min_z)

    def encode_coordinate(point):
        x, y, z = point
        return f"{x - min_x:02d}{y - min_y:02d}{z - min_z:02d}"

    piece_names = [f"P{piece.id}" for piece in pieces]
    placement_rows = set()

    for piece in pieces:
        piece_name = f"P{piece.id}"
        expected_coverage = 4 if piece.id == special_piece.id else 5
        rotations = {
            frozenset((coord.x, coord.y, coord.z) for coord in rotation.coords)
            for rotation in piece.rotations
        }

        for shape in rotations:
            (_, _, _), (shape_max_x, shape_max_y, shape_max_z) = (
                find_coordinate_bounds(shape)
            )
            for translate_x in range(-shape_max_x, wall_widths[0] + 1):
                for translate_y in range(-shape_max_y, wall_widths[1] + 1):
                    for translate_z in range(-shape_max_z, wall_widths[2] + 1):
                        translated_shape = frozenset(
                            (x + translate_x, y + translate_y, z + translate_z)
                            for x, y, z in shape
                        )
                        covered_cells = translated_shape.intersection(wall)
                        if len(covered_cells) == expected_coverage:
                            placement_rows.add(
                                (piece_name, frozenset(covered_cells))
                            )

    wall_items = " ".join(encode_coordinate(point) for point in sorted(wall))
    item_header = f"{wall_items} {' '.join(piece_names)}"
    lines = [
        "| 5x6x6 wall with 2x2 extension",
        f"| Piece P{special_piece.id} uses 4 wall cells",
        item_header,
    ]
    lines.extend(
        f"{piece_name} " + " ".join(encode_coordinate(point) for point in sorted(cells))
        for piece_name, cells in sorted(placement_rows)
    )
    return "\n".join(lines) + "\n", len(placement_rows), wall_widths


def main():
    pieces = list(load_pentacubes())
    if not SOLVER_PATH.is_file():
        raise FileNotFoundError(f"DLX solver executable not found: {SOLVER_PATH}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    for extension in EXTENSIONS[:1]:
        inner_cells = INNER_CELLS.union(extension)
        wall = build_wall(inner_cells)

        print(f"Wall: {len(wall)} cells", flush=True)

        for special_piece in pieces:
            input_text, placement_count, wall_widths = build_dlx_input(
                wall, pieces, special_piece
            )
            print(
                f"  Four-cell piece P{special_piece.id}: {placement_count} rows; "
                f"bounds {wall_widths}",
                flush=True,
            )
            result = subprocess.run(
                [str(SOLVER_PATH), "m1", "t1"],
                input=input_text,
                capture_output=True,
                text=True,
                check=False,
            )
            if result.returncode:
                raise RuntimeError(
                    f"ssxcc failed with exit code {result.returncode}:\n{result.stderr}"
                )

            if result.stdout.strip():
                DLX_PATH.write_text(
                    input_text + "=======\n" + result.stdout + "====\n",
                    encoding="utf-8",
                )
                subprocess.run(
                    [
                        sys.executable,
                        str(POSTPROCESSOR_PATH),
                        str(DLX_PATH),
                        str(JSON_PATH),
                    ],
                    check=True,
                )
                return

            summary_lines = result.stderr.strip().splitlines()
            if summary_lines:
                print(f"    {summary_lines[-1]}", flush=True)

    print("No solution found for the first wall geometry.")


if __name__ == "__main__":
    main()
