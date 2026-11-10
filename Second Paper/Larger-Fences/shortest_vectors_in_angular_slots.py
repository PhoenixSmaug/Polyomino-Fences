#!/usr/bin/env python3
"""Shortest and second-shortest integer vectors in each angle slot.

The range 0..45 degrees is split into N equal slots (default 225, i.e.
0.2 degrees each). Slot k holds the directions with angle in
[45k/N, 45(k+1)/N); the last slot also holds 45 degrees itself. For every
slot the program reports the shortest and the second-shortest integer
vector (x, y), 0 <= y <= x, whose direction lies in the slot, with its
length, slope and angle.

By default all integer vectors count, so the second-shortest can be twice
the shortest (marked "multiple"). With --primitive only vectors with
gcd(x, y) = 1 count, i.e. only distinct directions. Vectors of equal length
are ordered by angle.

No lattice direction lies exactly on an inner slot boundary (tan of a
rational number of degrees is irrational except at 0 and 45 degrees), so
floating-point angles put every vector in the right slot.

Usage:  python3 shortest_vectors.py [-N 225] [--primitive] [--csv FILE]
"""

import argparse
import csv
import math
import sys


def slot_vectors(N, primitive=False, per_slot=2):
    """Return a list of N lists, each with the `per_slot` shortest vectors
    (x, y) of that slot, sorted by length (then angle)."""
    R = 8
    while True:
        slots = [[] for _ in range(N)]
        R2 = R * R
        for x in range(1, R + 1):
            for y in range(0, x + 1):
                if x * x + y * y > R2:
                    break
                if primitive and math.gcd(x, y) != 1:
                    continue
                a = math.degrees(math.atan2(y, x))
                k = min(int(a / 45.0 * N), N - 1)
                slots[k].append((x * x + y * y, a, x, y))
        # Every vector of length <= R has x <= R and was enumerated, so a
        # slot with at least per_slot vectors inside the radius is final.
        if all(len(s) >= per_slot for s in slots):
            return [[(x, y) for _, _, x, y in sorted(s)[:per_slot]]
                    for s in slots]
        R *= 2


def describe(x, y, first=None):
    g = math.gcd(x, y)
    tag = ""
    if first is not None and (x, y) != first and g > 1:
        tag = "multiple"
    return {
        "x": x, "y": y,
        "length": math.hypot(x, y),
        "slope": f"{y // g}/{x // g}",
        "angle": math.degrees(math.atan2(y, x)),
        "note": tag,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("-N", type=int, default=225, help="number of slots")
    ap.add_argument("--primitive", action="store_true",
                    help="only vectors with gcd(x, y) = 1")
    ap.add_argument("--csv", metavar="FILE",
                    help="also write the table to this CSV file")
    args = ap.parse_args()
    if args.N < 1:
        sys.exit("N must be positive")

    N = args.N
    w = 45.0 / N
    res = slot_vectors(N, args.primitive)
    rows = []
    print(f"{N} slots of {w:g} deg"
          f"{', primitive vectors only' if args.primitive else ''}")
    print(" slot  angle range (deg)    shortest            length   angle"
          "     second              length   angle")
    for k, vs in enumerate(res):
        a = describe(*vs[0])
        b = describe(*vs[1], first=vs[0])
        rows.append((k, k * w, (k + 1) * w, a, b))
        print(f"{k:5d}  {k * w:8.3f}-{(k + 1) * w:7.3f}  "
              f"({a['x']:3d},{a['y']:3d}) {a['slope']:>7}  {a['length']:8.3f} "
              f"{a['angle']:7.3f}   "
              f"({b['x']:3d},{b['y']:3d}) {b['slope']:>7}  {b['length']:8.3f} "
              f"{b['angle']:7.3f}  {b['note']}")

    if args.csv:
        with open(args.csv, "w", newline="") as f:
            wr = csv.writer(f)
            wr.writerow(["slot", "from_deg", "to_deg",
                         "x1", "y1", "length1", "angle1",
                         "x2", "y2", "length2", "angle2", "second_is_multiple"])
            for k, lo, hi, a, b in rows:
                wr.writerow([k, f"{lo:.6g}", f"{hi:.6g}",
                             a["x"], a["y"], f"{a['length']:.6f}",
                             f"{a['angle']:.6f}",
                             b["x"], b["y"], f"{b['length']:.6f}",
                             f"{b['angle']:.6f}", int(b["note"] == "multiple")])


if __name__ == "__main__":
    main()
