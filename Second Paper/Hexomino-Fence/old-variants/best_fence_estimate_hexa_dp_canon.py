"""
Canonical-order variant of the dynamic program (port of best-fence-estimate-hexa-step3-count2.py
to the module structure of best_fence_estimate_hexa_dp.py).  Used by
best-fence-estimate-hexa-step4-canon.py.

WHY THIS EXISTS
  best_fence_estimate_hexa_dp.reconstruct_dag counts every ORDERING of identical pieces
  as a separate path (11 copies of (4,1) over 8 directions: millions of paths per polygon).
  Here the order is fixed inside the DP itself, so that
        one path in the DAG  ==  one oriented polygon.

DIFFERENCES TO best_fence_estimate_hexa_dp.forward_dp  (all taken from count2)
  1. State = (key, last_octant).  Consecutive IDENTICAL pieces must be placed into octants
     in weakly increasing order (octants 0..7 = right, up-right, up, up-left, left, ...,
     see 'moves').  last_octant is remembered only between identical pieces
     (afterwards it is clipped to 0 = "don't care").
  2. NO normalisation of the sum vector into the first octant (this would conflict with 1,
     because the octant order is only meaningful in a fixed frame).
  3. Symmetry breaking instead: if the start configuration is fully symmetric (e.g. four equal
     diagonal corner vectors) and the steepest piece (the first one processed) is unique and
     lies neither on an axis nor on the diagonal, that piece is forced into octant 7.
     Then every oriented polygon is generated exactly once.  Otherwise (not fully symmetric
     start, or no suitable piece) every polygon may appear in several images; the final
     identification under symmetries (dp.canonical) takes care of this.
  4. Each table entry carries a COUNT (number of optimal paths), so the number of solutions
     is known right after the forward pass.
The prune/reachability logic and the numbers (area_threshold etc.) are those of dp.py;
the pruning functions of dp.py are reused through the module's globals.
"""
import itertools
from math import sqrt, atan2
from types import SimpleNamespace

import best_fence_estimate_hexa_dp as dp
from best_fence_estimate_hexa_dp import det, rotate0, reflect


def _same_as_previous(Vecs, i):
    """UU of count2: 8 if the next piece to be processed in the forward pass (i-1) is identical to piece i,
    i.e. if last_octant must be remembered; 0 otherwise (clip last_octant to 0)."""
    return 8 if (i > 0 and Vecs[i] == Vecs[i-1]) else 0


def moves(x0, y0, key, last_octant):
    """All ways to add the piece (x0,y0) to the four bundles 'key' (same order and guards as
    in the forward DP).  Yields (octant, new_key, added_volume, vector)."""
    (x1,y1),(x2,y2),(x3,y3),(x4,y4) = key
    gen = (y0 != x0 and y0 != 0)   # the 'diagonal' octants only for generic pieces
    if last_octant <= 0:
        x,y = x0,y0
        yield 0, ((x1+x,y1+y),(x2,y2),(x3,y3),(x4,y4)), det((x,y),(x1,y1)), (x,y)
    if last_octant <= 1 and gen:
        x,y = y0,x0
        yield 1, ((x1+x,y1+y),(x2,y2),(x3,y3),(x4,y4)), det((x1,y1),(x,y)), (x,y)
    if last_octant <= 2:
        x,y = -y0,x0
        yield 2, ((x1,y1),(x2+x,y2+y),(x3,y3),(x4,y4)), det((x,y),(x2,y2)), (x,y)
    if last_octant <= 3 and gen:
        x,y = -x0,y0
        yield 3, ((x1,y1),(x2+x,y2+y),(x3,y3),(x4,y4)), det((x2,y2),(x,y)), (x,y)
    if last_octant <= 4:
        x,y = -x0,-y0
        yield 4, ((x1,y1),(x2,y2),(x3+x,y3+y),(x4,y4)), det((x,y),(x3,y3)), (x,y)
    if last_octant <= 5 and gen:
        x,y = -y0,-x0
        yield 5, ((x1,y1),(x2,y2),(x3+x,y3+y),(x4,y4)), det((x3,y3),(x,y)), (x,y)
    if last_octant <= 6:
        x,y = y0,-x0
        yield 6, ((x1,y1),(x2,y2),(x3,y3),(x4+x,y4+y)), det((x,y),(x4,y4)), (x,y)
    if gen:
        x,y = x0,-y0
        yield 7, ((x1,y1),(x2,y2),(x3,y3),(x4+x,y4+y)), det((x4,y4),(x,y)), (x,y)


def forward_dp_canonical(Vecs, Corners, keep_tables):
    """Forward pass.  TABLES[i] = list of (key, last_octant, vol, count)."""
    num = len(Vecs)

    # start configurations of the corners, exactly as in dp.forward_dp
    Start_configurations = set()
    for refls4 in itertools.product((rotate0, reflect), repeat=4):
        for c4 in itertools.permutations(Corners):
            [v1],[v2],[v3],[v4] = [ref([c]) for ref,c in zip(refls4,c4)]
            (x1,y1)=v1
            (x2,y2)=v2; x2=-x2
            (x3,y3)=v3; x3=-x3; y3=-y3
            (x4,y4)=v4; y4=-y4
            P2 = (x1,y1),(x2,y2),(x3,y3),(x4,y4)
            Q = frozenset(frozenset(r1(r2(P2))) for r1,r2 in dp.ALL_SYMMETRIES)
            Start_configurations.add(Q)
    print(len(Start_configurations), "start configurations of corners:")

    # symmetry breaking (see module docstring, point 3)
    fully_symmetric = (len(Start_configurations) == 1 and
                       all(len(Q) == 1 for Q in Start_configurations))
    xs, ys = Vecs[-1] if num else (0, 0)
    force_last_octant = (fully_symmetric and num > 0 and xs != ys and ys != 0
                         and (num == 1 or Vecs[-1] != Vecs[-2]))
    l_start = 7 if force_last_octant else 0
    print(f"fully symmetric start: {fully_symmetric};  "
          f"steepest piece {Vecs[-1] if num else None} forced into octant 7: {force_last_octant}")

    TABLES = dict()
    start = dict()
    for Q in Start_configurations:
        Q0 = list(Q)[0]
        print(*Q0)
        if Q0 == {(0,0)}:
            start[((0,0),(0,0),(0,0),(0,0)), l_start] = (0, 1)
        else:
            [(x1,y1)] = [(x,y) for x,y in Q0 if x>0 and y>0]
            [(x2,y2)] = [(x,y) for x,y in Q0 if x<0 and y>0]
            [(x3,y3)] = [(x,y) for x,y in Q0 if x<0 and y<0]
            [(x4,y4)] = [(x,y) for x,y in Q0 if x>0 and y<0]
            start[((x1,y1),(x2,y2),(x3,y3),(x4,y4)), l_start] = (0, 1)
    Table = TABLES[num] = [(key,l,vol,cnt) for (key,l),(vol,cnt) in start.items()]
    print("start table:", Table)
    dp.angle_area_store = dict()

    # weak meet-in-the-middle, as in dp.py
    NUM_FROM_END = min(num//2+6, num)
    reachablevectors = {0: {(0,0)}}
    for i in range(NUM_FROM_END):
        (x0,y0) = Vecs[i]
        all_directions = set((z for (x1,y1) in ((x0,y0),(x0,-y0),(-x0,y0),(-x0,-y0))
                              for z in ((x1,y1),(y1,x1))))
        reachablevectors[i+1] = set(((x+u,y+v) for (x,y) in reachablevectors[i]
                                     for (u,v) in all_directions))
        print(f"len(reachablevectors[{i+1}])={len(reachablevectors[i+1])}")

    for i in reversed(range(num)):
        (x0,y0) = Vecs[i]
        UU = _same_as_previous(Vecs, i)
        fixed_vectors = [sqrt(x**2+y**2) for x,y in Vecs[:i]]
        dp.sum_fixed_vectors = sum(fixed_vectors)
        dp.angle_area_store = dict()
        if i == NUM_FROM_END:
            print("reachability information kicks in.")
        check = (i <= NUM_FROM_END)
        reachable = reachablevectors[i] if check else None
        newTable = dict()

        def insert(octant, key, new_value, count):
            if check and (sum(x for x,y in key), sum(y for x,y in key)) not in reachable:
                return   # cannot be completed
            k2 = (key, octant if octant < UU else UU)      # min(octant,UU)
            old = newTable.get(k2)
            if old is None or new_value > old[0]:
                newTable[k2] = (new_value, count)
            elif new_value == old[0]:
                newTable[k2] = (new_value, count + old[1])

        for ((x1,y1),(x2,y2),(x3,y3),(x4,y4)), last_octant, vol, count in Table:
            key = ((x1,y1),(x2,y2),(x3,y3),(x4,y4))
            for octant, newkey, add, _vec in moves(x0, y0, key, last_octant):
                insert(octant, newkey, vol+add, count)
        beforepruning = len(newTable)
        print(f"iteration {i=} piece={(x0,y0)}, len(Table) before pruning: {beforepruning:,}", flush=True)
        Table = [(key,l,val,c) for (key,l),(val,c) in newTable.items()
                 if not dp.prune_solution(fixed_vectors, [sqrt(x**2+y**2) for x,y in key], val/2)]
        TABLES[i] = Table
        if not keep_tables:
            del TABLES[i+1]
        print(f"iteration {i=} piece={(x0,y0)}, {len(Table)=:,} (before pruning: {beforepruning:,})", flush=True)
    return TABLES


def run_forward_dp_canonical(Vecs, Corners, threshold, keep_tables):
    """Returns (best, TABLES, finals, numbest); finals = list of (key, last_octant, vol);
    numbest = number of optimal solutions according to the counts of the forward pass."""
    dp.area_threshold = threshold
    dp.bin_search_calls = dp.bin_search_decisions = 0
    TABLES = forward_dp_canonical(Vecs, Corners, keep_tables)
    print(f"Total number of binary search calls = {dp.bin_search_calls:,},",
          f"average number of decisions per call =",
          f"{dp.bin_search_decisions/max(dp.bin_search_calls,1):.5}")
    Table = TABLES[0]
    def total(key, vol):
        (x1,y1),(x2,y2),(x3,y3),(x4,y4) = key
        if (x1+x4+x3+x2, y1+y4+y3+y2) != (0,0):
            return None
        return vol + det((x1,y1),(x2,y2)) + det((x3,y3),(x4,y4))
    totals = [total(key,vol) for key,l,vol,c in Table]
    best = max((t for t in totals if t is not None), default=None)
    if best is None:
        return None, TABLES, [], 0
    finals, numbest = [], 0
    for (key,l,vol,count),t in zip(Table, totals):
        if t == best:
            finals.append((key,l,vol))
            numbest += count
    print(f"{len(finals)} best final states, {numbest} optimal solutions counted in the forward pass.")
    print(best, "best value", flush=True)
    return best, TABLES, finals, numbest


def solve_vector_set(nline, vecs, neg_penalty, target_area, keep_tables):
    """Same interface as dp.solve_vector_set, canonical-order version.  Additional field: count."""
    Vecs, Corners, shift, penalty = dp.prepare_instance(nline, vecs, neg_penalty)
    best, TABLES, finals, numbest = run_forward_dp_canonical(Vecs, Corners, target_area + penalty, keep_tables)
    enclosed_area = None
    reaches_goal = False
    if best is None:
        print(f"No solution with area at least {target_area} was found.")
    else:
        enclosed_area = best/2 - penalty
        print(f"best solution value {best/2}. **** enclosed area={enclosed_area} ************")
        if enclosed_area <= target_area:
            print(f"This is not guaranteed to be the optimum, because it is < {target_area=}")
        else:
            reaches_goal = True
            print("This is guaranteed to be the optimum.")
        print(flush=True)
    return SimpleNamespace(Vecs=Vecs, Corners=Corners, shift=shift, penalty=penalty, best=best,
                           TABLES=TABLES, finals=finals, count=numbest,
                           enclosed_area=enclosed_area, reaches_goal=reaches_goal)


def reconstruct_dag_canonical(Vecs, TABLES, finals):
    """DAG of all optimal canonical solutions.
    SOLS[i]: dict node -> [number of paths from the final states to this node, [(successor node, vector)]]
    node = (key, last_octant, vol); successors live in SOLS[i-1]... (SOLS[0] = final states).
    One full path == one oriented polygon, so the counts are the numbers of solutions."""
    num = len(Vecs)
    SOLS = {0: {(key,l,vol): [1, []] for key,l,vol in finals}}
    for i in range(num):
        (x0,y0) = Vecs[i]
        UU = _same_as_previous(Vecs, i)
        SOLS_old = SOLS[i]
        print(f"piece {i}  #entries = {len(SOLS_old)}  "
              f"#partial solutions = {sum(e[0] for e in SOLS_old.values()):,}", flush=True)
        SOLS_new = SOLS[i+1] = dict()
        for oldkey, old_l, oldvol, _cnt in TABLES[i+1]:
            for octant, key, add, vec in moves(x0, y0, oldkey, old_l):
                node = (key, octant if octant < UU else UU, oldvol+add)
                entry = SOLS_old.get(node)
                if entry is None:
                    continue
                e = SOLS_new.setdefault((oldkey,old_l,oldvol), [0, []])
                e[0] += entry[0]
                e[1].append((node, vec))
    print(f"piece {num}  #entries = {len(SOLS[num])}  "
          f"#partial solutions = {sum(e[0] for e in SOLS[num].values()):,}", flush=True)
    return SOLS


def _walk(SOLS, level, succlist, prefix):
    for succ, vec in succlist:
        if level == 0:
            yield prefix + [vec]
        else:
            yield from _walk(SOLS, level-1, SOLS[level][succ][1], prefix + [vec])


def enumerate_canonical_solutions(SOLS, num, shift, max_solutions):
    """Generate all solutions of the DAG (at most max_solutions different ones up to symmetry), print them.
    Returns dict canonical solution -> number of oriented solutions (images) found."""
    num_sols = sum(e[0] for e in SOLS[num].values())
    print(num_sols, "oriented solutions in the DAG (before identifying symmetric images).")
    found = {}
    jj = 0
    stop = False
    for (corn,_l,_vol),(_n,succlist) in SOLS[num].items():
        corn2 = []
        for x,y in corn:
            x -= shift if x>0 else -shift
            y -= shift if y>0 else -shift
            corn2.append((x,y))
        for pieces in _walk(SOLS, num-1, succlist, []):
            jj += 1
            sol = pieces + corn2
            Q = dp.canonical(sol)
            found[Q] = found.get(Q, 0) + 1
            if jj <= 50:
                print(jj, f"({list(found).index(Q)})", sol)
            if jj == 51:
                print("...\n")
            if len(found) >= max_solutions:
                print("aborting after", len(found), "solutions were found.")
                stop = True
                break
        if stop:
            break
    assert stop or jj == num_sols, (jj, num_sols)
    print("multiplicities (images per symmetry class)", dict(enumerate(found.values())))
    print("Solutions:")
    for i,S in enumerate(found):
        SS = [(x,y) for _,x,y in sorted((atan2(y,x),x,y) for x,y in S)]
        print(i, ",", SS, ",")
    return found
