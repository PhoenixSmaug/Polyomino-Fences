import logging
from functools import cache
from pathlib import Path
from random import Random
from typing import Sequence

from polycubetools import (
    Coordinate, CuboidHull, HullGrid, ParallelSolver, Placement, Polycube,
    PolycubePlacement, RotatedPolycube, SIX_NEIGHBORHOOD_DIRECTIONS, SNAPSHOT, utils,
)


HULL_MIN = Coordinate(3, 3, 3)
HULL_MAX = Coordinate(15, 15, 15)
GRID_SIZE = Coordinate(19, 19, 19)
TARGET_VOLUME = 11 ** 3

# These values have been determined empirically to find a 11 ** 3 solution
START_ANCHOR = Coordinate(15, 13, 11)
FIRST_PIECE_ID = 14
TIE_BREAK_SEED = 7

FRONTIER_BONUS = 2
INVALID_SCORE = -1000
INVALID_POSITION = Coordinate(-100, -100, -100)
OUTPUT_FILE = Path(__file__).parent / "output" / "hexa_greedy.json"


def scorer(hull_snap: SNAPSHOT, grid_snap: SNAPSHOT, candidate_snap: SNAPSHOT) -> tuple[SNAPSHOT, int]:
    hull = CuboidHull.from_snapshot(hull_snap)
    grid = HullGrid[CuboidHull].from_snapshot(grid_snap, hull)
    candidate = PolycubePlacement.from_snapshot(candidate_snap)
    frontier = grid.frontier

    @cache
    def cell_score(pos: Coordinate) -> int:
        contacts = sum(pos + d in grid for d in SIX_NEIGHBORHOOD_DIRECTIONS)
        return contacts + FRONTIER_BONUS * (pos in frontier)

    def _score(rpc: RotatedPolycube, position: Coordinate) -> int:
        score = 0
        for coord in rpc.coords:
            pos = position + coord
            if pos in grid or not grid.in_bounds(pos) or grid.in_inner_part(pos):
                return INVALID_SCORE
            score += cell_score(pos)
        return score

    rotations = sorted(candidate.polycube.rotations,
                       key=lambda rpc: sorted(c.to_tuple() for c in rpc.coords))
    anchors = sorted(grid.active_frontier, key=Coordinate.to_tuple)

    best = Placement(rotations[0], INVALID_POSITION)
    best_score = INVALID_SCORE
    for rpc in rotations:
        for anchor in anchors:
            for coord in sorted(rpc.coords, key=Coordinate.to_tuple):
                position = anchor - coord
                score = _score(rpc, position)
                if score > best_score:
                    best_score = score
                    best = Placement(rpc, position)

    return best.to_snapshot(), best_score


class GreedySolver(ParallelSolver[CuboidHull]):
    def __init__(self, grid: HullGrid[CuboidHull], polycubes: tuple[Polycube, ...],
                 first: int = FIRST_PIECE_ID) -> None:
        super().__init__(grid, polycubes, scorer)
        self.first: int | None = first
        self.order = sorted(self.remaining)
        Random(TIE_BREAK_SEED).shuffle(self.order)

    def next_candidates(self) -> Sequence[PolycubePlacement]:
        if not self.grid.active_frontier:
            return []

        anchor = min(self.grid.active_frontier, key=Coordinate.to_tuple)
        ids = [self.first] if self.first is not None else self.order
        self.first = None
        return [PolycubePlacement(self.remaining[id_], anchor)
                for id_ in ids if id_ in self.remaining]


def validate_solution(grid: HullGrid[CuboidHull]) -> int:
    if grid.frontier:
        raise RuntimeError(f"Incomplete shell: {len(grid.frontier)} cells uncovered")
    volume = utils.compute_volume(set(grid._cells))
    if volume != TARGET_VOLUME:
        raise RuntimeError(f"Expected volume {TARGET_VOLUME}, found {volume}")
    return volume


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    grid = HullGrid(GRID_SIZE, CuboidHull(HULL_MIN, HULL_MAX))
    grid.active_frontier = {START_ANCHOR}
    solver = GreedySolver(grid, utils.load_hexacubes())
    solver.run(output_file=None)

    volume = validate_solution(grid)
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    grid.to_json(OUTPUT_FILE)
    print(f"Enclosed volume: {volume}; solution saved to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
