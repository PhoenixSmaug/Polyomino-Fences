import logging
from collections import Counter
from pathlib import Path
from random import Random
from typing import Iterable, Sequence

import numpy as np
from polycubetools import (
    Coordinate, CuboidHull, HullGrid, ParallelSolver, Placement, Polycube,
    PolycubePlacement, SIX_NEIGHBORHOOD_DIRECTIONS, SNAPSHOT, utils,
)


INNER_SIZE = Coordinate(28, 28, 33)
HULL_MIN = Coordinate(6, 6, 6)
HULL_MAX = HULL_MIN + INNER_SIZE + Coordinate(1, 1, 1)
GRID_SIZE = HULL_MIN + HULL_MAX + Coordinate(1, 1, 1)
MINIMUM_VOLUME = 25544
TARGET_VOLUME = INNER_SIZE.x * INNER_SIZE.y * INNER_SIZE.z

START_ANCHOR = HULL_MIN
FIRST_PIECE_ID: int | None = None
TIE_BREAK_SEED = 1
REORDER_AFTER = 800
TAIL_TIE_BREAK_SEED = 12

FRONTIER_BONUS = 4
FACE_CAPACITY_PENALTY = 4
INVALID_SCORE = -1000
INVALID_POSITION = Coordinate(-100, -100, -100)
PROGRESS_INTERVAL = 50
OUTPUT_FILE = Path(__file__).parent / "output" / "hepta_greedy.json"


def scorer(_hull_snap: SNAPSHOT, score_snap: SNAPSHOT, candidate_snap: SNAPSHOT) -> tuple[SNAPSHOT, int]:
    candidate = PolycubePlacement.from_snapshot(candidate_snap)
    rotations = sorted(candidate.polycube.rotations, key=lambda rpc: sorted(c.to_tuple() for c in rpc.coords))
    size = np.array(score_snap["size"])
    strides = np.array((size[1] * size[2], size[2], 1))
    anchors = np.array(score_snap["anchors"])
    field = np.frombuffer(score_snap["scores"], dtype=np.int16)
    best = Placement(rotations[0], INVALID_POSITION)
    best_score = INVALID_SCORE

    for rpc in rotations:
        coords = np.array(sorted(c.to_tuple() for c in rpc.coords))
        positions = (anchors[:, None, :] - coords).reshape(-1, 3)
        in_bounds = np.all((positions >= 0) & (positions + coords.max(axis=0) < size), axis=1)
        positions = positions[in_bounds]
        if not len(positions):
            continue
        cells = (positions @ strides)[:, None] + coords @ strides
        cell_scores = field[cells]
        scores = cell_scores.sum(axis=1)
        scores[np.any(cell_scores < 0, axis=1)] = INVALID_SCORE
        index = int(scores.argmax())
        if scores[index] > best_score:
            best_score = int(scores[index])
            best = Placement(rpc, Coordinate(*map(int, positions[index])))

    if best_score != INVALID_SCORE:
        points = [c.to_tuple() for c in candidate.polycube.rotations[0].coords]
        face_capacity = max(count for axis in range(3) for count in Counter(p[axis] for p in points).values())
        best_score -= FACE_CAPACITY_PENALTY * face_capacity
    return best.to_snapshot(), best_score


class GreedySolver(ParallelSolver[CuboidHull]):
    def __init__(self, grid: HullGrid[CuboidHull], polycubes: tuple[Polycube, ...],
                 first: int | None = FIRST_PIECE_ID) -> None:
        super().__init__(grid, polycubes, scorer)
        self.first = first
        self.order = sorted(self.remaining)
        Random(TIE_BREAK_SEED).shuffle(self.order)
        self._base_field = np.zeros(grid.size.to_tuple(), dtype=np.int16)
        frontier = np.array([c.to_tuple() for c in grid.frontier])
        self._base_field[tuple(frontier.T)] = FRONTIER_BONUS
        lo, hi = grid.hull.min_coord, grid.hull.max_coord
        self._base_field[lo.x + 1:hi.x, lo.y + 1:hi.y, lo.z + 1:hi.z] = INVALID_SCORE

    def next_candidates(self) -> Sequence[PolycubePlacement]:
        if not self.grid.active_frontier:
            return []
        anchor = min(self.grid.active_frontier, key=Coordinate.to_tuple)
        ids = [self.first] if self.first is not None else self.order
        self.first = None
        return [PolycubePlacement(self.remaining[id_], anchor) for id_ in ids if id_ in self.remaining]

    def _score_snapshot(self) -> SNAPSHOT:
        field = self._base_field.copy()
        occupied = np.array([c.to_tuple() for c in self.grid._cells], dtype=np.int16).reshape(-1, 3)
        for direction in SIX_NEIGHBORHOOD_DIRECTIONS:
            neighbors = occupied + direction.to_tuple()
            valid = np.all((neighbors >= 0) & (neighbors < self.grid.size.to_tuple()), axis=1)
            field[tuple(neighbors[valid].T)] += 1
        field[tuple(occupied.T)] = INVALID_SCORE
        return {
            "size": self.grid.size.to_tuple(),
            "anchors": sorted(c.to_tuple() for c in self.grid.active_frontier),
            "scores": field.tobytes(),
        }

    def _compute_best(self, candidates: Iterable[PolycubePlacement]) -> tuple[Placement | None, int]:
        snapshot = self._score_snapshot()
        futures = [self.score(({}, snapshot, candidate.to_snapshot())) for candidate in candidates]
        placement, score = max((future.result() for future in futures), key=lambda result: result[1])
        if score == INVALID_SCORE:
            return None, score
        return Placement.from_snapshot(placement), score

    def after_place(self, placement: Placement) -> None:
        super().after_place(placement)
        placed = len(self.polycubes) - len(self.remaining)
        if placed == REORDER_AFTER:
            Random(TAIL_TIE_BREAK_SEED).shuffle(self.order)
        if placed % PROGRESS_INTERVAL == 0:
            self._logger.info("Placed %d pieces; %d shell cells remain", placed, len(self.grid.frontier))


def validate_solution(grid: HullGrid[CuboidHull]) -> int:
    if grid.frontier:
        raise RuntimeError(f"Incomplete shell: {len(grid.frontier)} cells uncovered")
    volume = utils.compute_volume(set(grid._cells))
    if volume != TARGET_VOLUME or volume < MINIMUM_VOLUME:
        raise RuntimeError(f"Expected volume {TARGET_VOLUME}, found {volume}")
    return volume


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    grid = HullGrid(GRID_SIZE, CuboidHull(HULL_MIN, HULL_MAX))
    grid.active_frontier = {START_ANCHOR}
    solver = GreedySolver(grid, utils.load_heptacubes())
    solver.run(output_file=None)

    volume = validate_solution(grid)
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    grid.to_json(OUTPUT_FILE)
    print(f"Enclosed volume: {volume}; solution saved to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
