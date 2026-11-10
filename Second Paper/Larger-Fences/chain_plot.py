#!/usr/bin/env python3
"""Polygonal chain of polyomino diameter vectors, compared with a 1/8 circle,
and the same comparison for convex polygons.

Vector files (data/vec*.csv):
1. Each line after the header holds n,x,y,d and stands for the vector
   (x/d, y/d); n should be 500.  All files together give one set of vectors.
2. Sort the vectors by increasing slope and concatenate them to a polygonal
   chain from (0,0) to (X,Y).
3. Centre C = (0, X+Y), radius R = X+Y (the distance from C to the start).
4. Seen from C, the chain spans 45 degrees, from the direction of (0,0) to the
   direction of (X,Y).  Split that range into k equal angles and intersect
   each of the k+1 rays with the chain.
5. Draw the curve through the k+1 sample points together with the 1/8 circle
   of radius R around C and two segments on the rays from C through the start
   and end point, without axes, on screen and as a plain TikZ file.
6. Plot (|P - C| - R)/R against the ray angle, 0 to 45 degrees measured from
   the direction of (0,0).
7. Print the area of the shape of all free n-ominoes (8 copies of the chain),
   scaled up from the samples with A'(n) = c*lambda^n / (8n).
--unbiased reads every entry as the vector ((x-1)/d, y/d) instead and adds
"-unbiased" to all output names.  Each run also writes a standalone LaTeX
document figures[-unbiased].tex that shows the two figures.

Polygon files (data/path*txt), one convex polygon per file, one vertex x,y
per line:
- C is the centre of the bounding rectangle; the section runs from the
  downward direction (0 deg) to the south-east direction (45 deg) seen from C.
- With more than k vertices in the section, the k+1 equally spaced rays are
  used; otherwise the vertices in the section, plus the rays at 0 and 45 deg.
- R is the distance from C to the boundary point straight below C, and
  (|P - C| - R)/R is plotted in the same figure as 6., on a second vertical
  scale shared by all polygons.  If the section holds fewer than --max-dots
  true corners (points where the boundary turns), they are marked by dots.
- --vec-scale F and --vec-shift S draw the black vector curve F times taller
  and with its zero where S is on the unshifted scale; its own axis still
  shows its true values, and it gets its own (dotted) zero line.

Usage:  python3 chain_plot.py [data_dir] [-k 1000] [-o out_dir]
                            [--rays IN OUT] [--vec-scale F]
                            [--vec-shift S] [--max-dots M] [--max-vec N]
                            [--unbiased] [--no-show]
"""

import argparse
import glob
import os
import re
import sys
from decimal import Decimal, getcontext

import numpy as np

EXPECTED_N = 500
JENSEN_LAMBDA, JENSEN_C = Decimal("4.0626"), Decimal("0.3169")  # Jensen (2000)
VEC_PATTERN = "vec*.csv"
PATH_PATTERN = "path*txt"


def read_vectors(files, max_samples=None):
    """Read the vec files in order.

    A line n,x,y,d is one vector (x,y) of a sample with weight 1/d: a sample
    whose diameter is attained by d vectors gives d consecutive lines with
    the same d.  With max_samples, reading stops after that many samples
    (sum of the weights 1/d), at a sample boundary.
    """
    rows = []
    have = 0.0
    for f in files:
        if max_samples is not None and have >= max_samples - 1e-9:
            break
        a = np.loadtxt(f, delimiter=",", skiprows=1, dtype=np.int64, ndmin=2)
        if not a.size:
            continue
        if max_samples is not None:
            cw = have + np.cumsum(1.0 / a[:, 3])
            a = a[cw <= max_samples + 1e-9]
        rows.append(a)
        have += np.sum(1.0 / a[:, 3])
    a = np.vstack(rows)
    check_groups(a[:, 3])
    if max_samples is not None:
        print(f"read {len(a)} lines = {have:g} samples "
              f"(limit {max_samples})", file=sys.stderr)
    n, x, y, d = a.T
    bad_n = n != EXPECTED_N
    if bad_n.any():
        print(f"warning: skipping {bad_n.sum()} rows with n != {EXPECTED_N}",
              file=sys.stderr)
    zero = (x == 0) & (y == 0)
    if (zero & ~bad_n).any():
        print(f"warning: skipping {(zero & ~bad_n).sum()} zero vectors",
              file=sys.stderr)
    keep = ~bad_n & ~zero
    x, y, d = x[keep], y[keep], d[keep]
    print(f"{len(files)} vector file(s), {keep.sum()} lines, "
          f"{n_samples(d):g} samples", file=sys.stderr)
    return x, y, d


def n_samples(d):
    """Number of samples: every line counts 1/d."""
    return float(np.sum(1.0 / d))


def check_groups(d):
    """Warn if the lines do not form groups of d consecutive lines with the
    same d."""
    i, bad = 0, 0
    while i < len(d):
        k = int(d[i])
        if k < 1 or i + k > len(d) or np.any(d[i:i + k] != k):
            bad += 1
            i += 1
        else:
            i += k
    if bad:
        print(f"warning: {bad} lines do not fit into groups of d equal "
              "consecutive d values", file=sys.stderr)


def build_chain(x, y, d):
    # Sort by slope y/x (d cancels). With integer x, y below a few thousand,
    # atan2(y, x) separates distinct slopes reliably; ties only reorder
    # collinear edges, which does not change the chain.
    order = np.argsort(np.arctan2(y, x), kind="stable")
    vx = x[order] / d[order]
    vy = y[order] / d[order]
    px = np.concatenate(([0.0], np.cumsum(vx)))
    py = np.concatenate(([0.0], np.cumsum(vy)))
    return px, py


def ray_hits(px, py, C, phi, theta):
    """Intersect the rays from C at angles theta with the chain px, py.

    phi are the polar angles of the chain vertices seen from C, increasing
    along the chain; every theta must lie in [phi[0], phi[-1]].
    """
    cx, cy = C
    if np.any(np.diff(phi) < -1e-12):
        sys.exit("curve is not star-shaped with respect to C; "
                 "ray intersections are not unique")
    phi = np.maximum.accumulate(phi)          # remove round-off wiggles
    # Segment i joins vertex i and i+1 and covers angles [phi_i, phi_{i+1}].
    i = np.searchsorted(phi, theta, side="right") - 1
    i = np.clip(i, 0, len(px) - 2)
    ax, ay = px[i] - cx, py[i] - cy
    bx, by = px[i + 1] - px[i], py[i + 1] - py[i]
    ux, uy = np.cos(theta), np.sin(theta)
    # Solve C + t*u = A + s*b for t (distance from C along the ray).
    den = ux * by - uy * bx
    with np.errstate(divide="ignore", invalid="ignore"):
        t = (ax * by - ay * bx) / den
    # A ray along a segment cannot occur when C is strictly inside the
    # convex hull side, but guard anyway.
    deg = np.abs(den) < 1e-15
    t[deg] = np.hypot(ax[deg], ay[deg])
    return cx + t * ux, cy + t * uy


def sample_chain(px, py, k):
    """Rays from C = (0, X+Y) through the vector chain."""
    X, Y = px[-1], py[-1]
    R = X + Y
    C = (0.0, R)
    phi = np.arctan2(py - C[1], px - C[0])
    theta = np.linspace(phi[0], phi[-1], k + 1)  # -90 deg .. direction of (X,Y)
    sx, sy = ray_hits(px, py, C, phi, theta)
    return theta, sx, sy, R, C


def read_polygon(path):
    P = np.loadtxt(path, delimiter=",", ndmin=2)
    if len(P) > 1 and np.allclose(P[0], P[-1]):
        P = P[:-1]
    return P[:, 0], P[:, 1]


def sample_polygon(qx, qy, k):
    """Deviation of a convex polygon from a circle, seen from the centre C of
    its bounding rectangle, over the 45 deg below C towards the south-east.

    Returns alpha (deg, 0 = downward), relative deviation, the number of
    polygon vertices in the section, and the true corners (points where the
    boundary turns, not points on a straight edge) in the section as
    (alpha, deviation) arrays.
    """
    C = ((qx.min() + qx.max()) / 2, (qy.min() + qy.max()) / 2)
    phi = np.arctan2(qy - C[1], qx - C[0])
    # Orient counter-clockwise, so the angle increases from 0 to 45 deg.
    area2 = np.sum(qx * np.roll(qy, -1) - np.roll(qx, -1) * qy)
    if area2 < 0:
        qx, qy, phi = qx[::-1], qy[::-1], phi[::-1]
    # Angle measured from the downward direction, in [0, 2 pi).
    rel = np.mod(phi + np.pi / 2, 2 * np.pi)
    # Start the chain at the vertex with the smallest angle and prepend the
    # vertex before it (angle shifted by -2 pi), so the chain covers
    # [rel_prev - 2 pi, 2 pi) monotonically and contains the angle 0.
    m = np.argmin(rel)
    order = np.roll(np.arange(len(qx)), -m)
    order = np.concatenate(([order[-1]], order, [order[0]]))
    cx_, cy_ = qx[order], qy[order]
    a = rel[order].copy()
    a[0] -= 2 * np.pi
    a[-1] += 2 * np.pi
    quarter = np.pi / 4
    inside = (a >= 0) & (a <= quarter)
    n_in = int(inside[1:-1].sum())
    if n_in > k:
        theta = np.linspace(0.0, quarter, k + 1)
    else:
        theta = np.unique(np.concatenate(([0.0], a[1:-1][inside[1:-1]],
                                          [quarter])))
    sx, sy = ray_hits(cx_, cy_, C, a - np.pi / 2, theta - np.pi / 2)
    dist = np.hypot(sx - C[0], sy - C[1])
    R = dist[0]                                  # boundary point below C
    # True corners: the cross product of the incoming and outgoing edge is
    # non-zero (relative to the edge lengths).
    ex, ey = np.roll(qx, -1) - qx, np.roll(qy, -1) - qy
    px_, py_ = np.roll(ex, 1), np.roll(ey, 1)
    cross = px_ * ey - py_ * ex
    corner = np.abs(cross) > 1e-12 * np.hypot(px_, py_) * np.hypot(ex, ey)
    ca = np.mod(np.arctan2(qy - C[1], qx - C[0]) + np.pi / 2, 2 * np.pi)
    sel = corner & (ca <= quarter)
    crel = (np.hypot(qx[sel] - C[0], qy[sel] - C[1]) - R) / R
    order_c = np.argsort(ca[sel])
    corners = (np.degrees(ca[sel][order_c]), crel[order_c])
    return np.degrees(theta), (dist - R) / R, n_in, corners


def path_n(path):
    """n from a file name like path-13.txt, or None."""
    m = re.search(r"(\d+)", os.path.basename(path)[len("path"):])
    return int(m.group(1)) if m else None


def path_key(path):
    """Sort path files numerically by n; files without a number go last."""
    n = path_n(path)
    return (n is None, n if n is not None else 0, os.path.basename(path))


def print_area(x, y, d, alpha, rel):
    """Area of the full shape = 8 sectors, each 1/2 R^2 int (1+delta)^2 dalpha.
    The N samples stand for A'(n) free n-ominoes, so lengths scale with
    A'(n)/N and R_total = A'(n) * mean(x+y); delta stays unchanged."""
    getcontext().prec = 30                  # numbers exceed the float range
    N = np.sum(1.0 / d)                     # tied rows (d > 1) weigh 1/d
    mean_xy = np.sum((x + y) / d) / N
    f = (1 + rel) ** 2
    I = float(np.sum((f[1:] + f[:-1]) / 2 * np.diff(np.radians(alpha))))
    A = JENSEN_C * JENSEN_LAMBDA ** EXPECTED_N / (8 * EXPECTED_N)   # A'(n)
    R = A * Decimal(float(mean_xy))
    print(f"A'(n) = {A:.4e}, R = {R:.4e}, area/(pi R^2) = {4 * I / np.pi:.6f}, "
          f"area of all free {EXPECTED_N}-ominoes = {4 * Decimal(I) * R * R:.4e}")


def fmt(v):
    return f"{v:.6g}"


def ray_segments(sx, sy, C, extent):
    """Segments on the rays from C through the start and the end point."""
    cx, cy = C
    segs = []
    for px, py in ((sx[0], sy[0]), (sx[-1], sy[-1])):
        segs.append([(cx + f * (px - cx), cy + f * (py - cy)) for f in extent])
    return segs


def write_tikz_curve(path, sx, sy, R, C, theta, extent):
    # Plain TikZ, no axes.  Coordinates are divided by X, so the chain is
    # one unit wide; the unit is set by x=, y= in the tikzpicture options.
    cx, cy = C
    s = 1.0 / sx[-1]

    def coords(xs, ys):
        return " ".join(f"({a * s:.5f},{b * s:.5f})" for a, b in zip(xs, ys))

    with open(path, "w") as f:
        f.write("% Chain sampled along k+1 rays from C (solid), 1/8 circle "
                "around C (dashed),\n% and segments on the rays from C "
                "through the start and end point.\n% Plain TikZ; change "
                "x= and y= to resize (the chain is 1 unit wide).\n")
        f.write("\\begin{tikzpicture}[x=0.9\\linewidth, y=0.9\\linewidth]\n")
        for seg in ray_segments(sx, sy, C, extent):
            f.write("\\draw[thin] plot coordinates {"
                    f"{coords(*zip(*seg))}}};\n")
        # The 1/8 circle as one arc from the start point (0,0), which lies at
        # -90 deg seen from C.  A radius without unit is measured in the
        # x/y units of the picture.
        f.write(f"\\draw[dashed, blue] ({cx * s:.5f},{(cy - R) * s:.5f}) "
                f"arc[start angle=-90, end angle=-45, radius={R * s:.6f}];\n")
        f.write("\\draw[thick] plot coordinates {\n")
        f.write(coords(sx, sy) + "};\n")
        f.write("\\end{tikzpicture}\n")


def padded(lo, hi, pad=0.05):
    lo, hi = min(lo, 0.0), max(hi, 0.0)
    span = (hi - lo) or 1.0
    return lo - pad * span, hi + pad * span


def aligned_limits(left, right):
    """Widen two (lo, hi) ranges so that 0 sits at the same height in both."""
    (l0, l1), (r0, r1) = left, right
    f = max(-l0 / (l1 - l0), -r0 / (r1 - r0))   # height of 0, from below
    def fit(lo, hi):
        # keep the larger side fixed, extend the other one
        if -lo / (hi - lo) < f:
            lo = -f * hi / (1 - f)
        else:
            hi = -lo * (1 - f) / f
        return lo, hi
    return fit(l0, l1), fit(r0, r1)


# Colours for the polygon curves (the vector curve is black).
PATH_COLORS = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e",
               "#8c564b", "#e377c2", "#17becf", "#bcbd22", "#7f7f7f"]


def tex_escape(s):
    return s.replace("_", "\\_").replace("%", "\\%").replace("#", "\\#")


def tick_opts(ticks):
    """Explicit tick positions and labels, copied from the matplotlib plot."""
    pos = ",".join(f"{v:.10g}" for v, _ in ticks)
    lab = ",".join("{$" + t + "$}" for _, t in ticks)
    return f"ytick={{{pos}}}, yticklabels={{{lab}}}"


def write_tikz_deviation(path, vec, paths, limits, shifted, ticks):
    """vec: (label, alpha, rel) or None; paths: list of
    (label, alpha, rel, dots), dots = (alpha, rel) of marked corners or None.
    shifted: the black curve has its own zero line.
    ticks: (left, right) lists of (value, label) for the vertical axes."""
    # scale only axis: the two overlaid axes get the same plot area, even
    # though only the first one has x tick labels, so the zeros line up.
    common = ("scale only axis, width=0.9\\linewidth, height=0.5\\linewidth, "
              "xmin=0, xmax=45, scaled y ticks=false")
    xaxis = ("xtick={0,5,...,45}, xticklabel={$\\pgfmathprintnumber{\\tick}"
             "^\\circ$}, xlabel={$\\phi$}")
    with open(path, "w") as f:
        f.write("% Relative deviation (|P-C|-R)/R against the ray angle.\n"
                "% Vector chain on the left scale, polygons on the right "
                "scale (zero aligned).\n% Needs \\usepackage{pgfplots} "
                "(\\pgfplotsset{compat=1.17} or later).\n")
        for c_i, c in enumerate(PATH_COLORS):
            f.write(f"\\definecolor{{pathcol{c_i}}}{{HTML}}{{{c[1:]}}}\n")
        f.write("\\begin{tikzpicture}\n")
        first = True
        if vec is not None:
            lo, hi = limits[0]
            right = ", axis y line*=left" if paths else ""
            f.write(f"\\begin{{axis}}[{common}, {xaxis}, ymin={lo:.6e}, "
                    f"ymax={hi:.6e}, {tick_opts(ticks[0])}{right}]\n")
            zstyle = "black, densely dotted" if shifted else "gray, dashed"
            f.write(f"\\addplot[{zstyle}, domain=0:45] {{0}};\n")
            f.write("\\addplot[black, thick] coordinates {\n")
            f.writelines(f"({a:.4f},{v:.6e})\n" for a, v in zip(*vec[1:]))
            f.write("};\n\\end{axis}\n")
            first = False
        if paths:
            lo, hi = limits[1]
            if first:
                opts = f"{common}, {xaxis}"
            else:
                opts = (f"{common}, axis y line*=right, axis x line=none")
            f.write(f"\\begin{{axis}}[{opts}, ymin={lo:.6e}, ymax={hi:.6e}, "
                    f"{tick_opts(ticks[1])}, "
                    "legend pos=north west, legend cell align=left]\n")
            if vec is not None:
                f.write("\\addlegendimage{black, thick}\n"
                        f"\\addlegendentry{{{tex_escape(vec[0])} (left)}}\n")
            if vec is None or shifted:
                f.write("\\addplot[gray, dashed, domain=0:45, forget plot] "
                        "{0};\n")
            for c_i, (label, alpha, rel, dots) in enumerate(paths):
                col = f"pathcol{c_i % len(PATH_COLORS)}"
                f.write(f"\\addplot[{col}, thick] coordinates {{\n")
                f.writelines(f"({a:.4f},{v:.6e})\n" for a, v in zip(alpha, rel))
                side = " (right)" if vec is not None else ""
                f.write(f"}};\n\\addlegendentry{{{tex_escape(label)}{side}}}\n")
                if dots is not None and len(dots[0]):
                    f.write(f"\\addplot[{col}, only marks, mark=*, "
                            "mark size=1.2pt, forget plot] coordinates {\n")
                    f.writelines(f"({a:.4f},{v:.6e})\n" for a, v in zip(*dots))
                    f.write("};\n")
            f.write("\\end{axis}\n")
        f.write("\\end{tikzpicture}\n")


def write_figures_doc(path, chain_files, deviation_file, unbiased):
    """Standalone LaTeX document with the figures of this run; it expects
    the figure files in the same directory."""
    which = (" (vectors $((x-1)/d,\\,y/d)$)" if unbiased else "")
    with open(path, "w") as f:
        f.write("\\documentclass{article}\n\\usepackage[margin=2.5cm]{geometry}\n"
                "\\usepackage{pgfplots}\n\\pgfplotsset{compat=1.17}\n"
                "\\begin{document}\n")
        for c in chain_files:
            f.write("\\begin{figure}[ht]\n\\centering\n"
                    f"\\input{{{c}}}\n"
                    f"\\caption{{Chain and circle{which}.}}\n"
                    "\\label{fig:chain}\n\\end{figure}\n")
        f.write("\\begin{figure}[ht]\n\\centering\n"
                f"\\input{{{deviation_file}}}\n"
                f"\\caption{{Relative deviation from the circle{which}.}}\n"
                "\\label{fig:deviation}\n\\end{figure}\n"
                "\\end{document}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("data_dir", nargs="?", default="data")
    ap.add_argument("-k", type=int, default=1000, help="number of angles")
    ap.add_argument("-o", "--out", default=".", help="output directory")
    ap.add_argument("--rays", type=float, nargs=2, default=(0.75, 1.03),
                    metavar=("IN", "OUT"),
                    help="ray segments run from IN to OUT times the distance "
                         "of the start/end point from C (default 0.75 1.03)")
    ap.add_argument("--vec-scale", type=float, default=1.0, metavar="F",
                    help="draw the black vec curve F times taller than on "
                         "the common scale (default 1)")
    ap.add_argument("--vec-shift", type=float, default=0.0, metavar="S",
                    help="put the zero of the black vec curve where S is on "
                         "the common scale (default 0)")
    ap.add_argument("--max-dots", type=int, default=100,
                    help="mark the true corners of a polygon in the section "
                         "if there are fewer than this (default 100)")
    ap.add_argument("--unbiased", action="store_true",
                    help="use the vectors ((x-1)/d, y/d); output files get "
                         "the suffix -unbiased")
    ap.add_argument("--max-vec", type=int, default=None, metavar="N",
                    help="stop reading the vec files after N samples "
                         "(each line counts 1/d; in file name order)")
    ap.add_argument("--no-show", action="store_true",
                    help="do not open plot windows, only write files")
    args = ap.parse_args()

    import matplotlib
    if args.no_show or not (os.environ.get("DISPLAY")
                            or os.environ.get("WAYLAND_DISPLAY")
                            or sys.platform in ("darwin", "win32")):
        matplotlib.use("Agg")
        args.no_show = True
    import matplotlib.pyplot as plt

    vec_files = sorted(glob.glob(os.path.join(args.data_dir, VEC_PATTERN)))
    path_files = sorted(glob.glob(os.path.join(args.data_dir, PATH_PATTERN)),
                        key=path_key)
    if not vec_files and not path_files:
        sys.exit(f"no {VEC_PATTERN} or {PATH_PATTERN} files in {args.data_dir}")
    os.makedirs(args.out, exist_ok=True)
    sfx = "-unbiased" if args.unbiased else ""

    vec = None
    if vec_files:
        x, y, d = read_vectors(vec_files, args.max_vec)
        if args.unbiased:
            x = x - 1
            zero = (x == 0) & (y == 0)
            if zero.any():
                print(f"warning: --unbiased: skipping {zero.sum()} vectors "
                      "that become zero", file=sys.stderr)
                x, y, d = x[~zero], y[~zero], d[~zero]
            steep = y > x
            if steep.any():
                print(f"warning: --unbiased: {steep.sum()} vectors have "
                      "slope > 1", file=sys.stderr)
        px, py = build_chain(x, y, d)
        theta, sx, sy, R, C = sample_chain(px, py, args.k)
        rel = (np.hypot(sx - C[0], sy - C[1]) - R) / R
        alpha = np.degrees(theta - theta[0])  # 0 deg at (0,0), 45 at (X,Y)
        print(f"vectors: end (X,Y) = ({px[-1]:.6g}, {py[-1]:.6g}), "
              f"R = X+Y = {R:.6g}, relative deviation "
              f"{rel.min():.3g} .. {rel.max():.3g}", file=sys.stderr)
        print_area(x, y, d, alpha, rel)
        label = f"{n_samples(d):g} random {EXPECTED_N}-ominoes"
        vec = (label, alpha, rel)

        write_tikz_curve(os.path.join(args.out, f"chain{sfx}.tex"), sx, sy, R, C,
                         theta, args.rays)
        arc_t = np.linspace(-np.pi / 2, -np.pi / 4, 400)
        fig1, ax1 = plt.subplots(figsize=(7, 4))
        for seg in ray_segments(sx, sy, C, args.rays):
            ax1.plot(*zip(*seg), color="black", lw=0.6)
        ax1.plot(C[0] + R * np.cos(arc_t), C[1] + R * np.sin(arc_t),
                 color="gray", ls="--", lw=1)
        ax1.plot(sx, sy, color="black", lw=1.5)
        ax1.set_aspect("equal")
        ax1.axis("off")
        fig1.tight_layout()
        fig1.savefig(os.path.join(args.out, f"chain{sfx}.png"), dpi=150)

    paths = []
    for p in path_files:
        qx, qy = read_polygon(p)
        alpha, rel, n_in, corners = sample_polygon(qx, qy, args.k)
        n_c = len(corners[0])
        dots = corners if n_c < args.max_dots else None
        n = path_n(p)
        label = (f"{n}-omino fence" if n is not None
                 else os.path.splitext(os.path.basename(p))[0])
        how = (f"{args.k + 1} equal rays" if n_in > args.k
               else f"{len(alpha)} rays through vertices and ends")
        print(f"{label}: {len(qx)} vertices, {n_in} in the section "
              f"({n_c} true corners{', marked' if dots else ''}), {how}, "
              f"relative deviation {rel.min():.3g} .. {rel.max():.3g}",
              file=sys.stderr)
        paths.append((label, alpha, rel, dots))

    # Vertical ranges: one for the vectors, one shared by all polygons,
    # with 0 at the same height when both are present.  The unscaled layout
    # fixes how many right-scale units one left-scale unit takes on screen
    # (rho).  The black curve is then drawn at screen height rho*(F*rel + S)
    # in right-scale units, the screen range is widened to hold everything,
    # and the left axis gets the limits that undo F and S, so its tick
    # labels still show the black curve's true values.
    F, S = args.vec_scale, args.vec_shift
    shifted = S != 0
    lim_v = lim_p = None
    if paths:
        lim_p = padded(min(r.min() for _, _, r, _ in paths),
                       max(r.max() for _, _, r, _ in paths))
    if vec:
        lim_v = padded(vec[2].min(), vec[2].max())
        if paths:
            lim_v, lim_p = aligned_limits(lim_v, lim_p)
            rho = (lim_p[1] - lim_p[0]) / (lim_v[1] - lim_v[0])
            screen = lim_p
        else:
            rho, screen = 1.0, lim_v
        b_lo, b_hi = padded(rho * (F * vec[2].min() + S),
                            rho * (F * vec[2].max() + S))
        screen = (min(screen[0], b_lo, rho * S), max(screen[1], b_hi, rho * S))
        if paths:
            lim_p = screen
        lim_v = ((screen[0] / rho - S) / F, (screen[1] / rho - S) / F)

    fig2, ax2 = plt.subplots(figsize=(8, 4))
    handles = []
    if vec and shifted:
        ax2.axhline(0, color="black", ls=":", lw=0.8)
    else:
        ax2.axhline(0, color="gray", ls="--", lw=0.8)
    if vec:
        handles += ax2.plot(vec[1], vec[2], color="black", lw=1.2,
                            label=vec[0] + (" (left)" if paths else ""))
        ax2.set_ylim(*lim_v)
    if paths:
        axp = ax2.twinx() if vec else ax2
        if vec and shifted:
            axp.axhline(0, color="gray", ls="--", lw=0.8)
        for c_i, (label, alpha, rel, dots) in enumerate(paths):
            col = PATH_COLORS[c_i % len(PATH_COLORS)]
            handles += axp.plot(alpha, rel, lw=1.2, color=col,
                                label=label + (" (right)" if vec else ""))
            if dots is not None:
                axp.plot(*dots, "o", ms=3, color=col)
        axp.set_ylim(*lim_p)
    ax2.legend(handles=handles, loc="upper left", frameon=False)
    for a_ in {ax2, axp} if paths else {ax2}:
        a_.ticklabel_format(axis="y", style="plain", useOffset=False)
    ax2.set_xlim(0, 45)
    ax2.set_xticks(range(0, 46, 5))
    ax2.xaxis.set_major_formatter("{x:.0f}°")
    ax2.set_xlabel("φ")
    fig2.tight_layout()
    fig2.savefig(os.path.join(args.out, f"deviation{sfx}.png"), dpi=150)

    def shown_ticks(a_):
        lo, hi = a_.get_ylim()
        return [(v, t.get_text().replace("\u2212", "-"))
                for v, t in zip(a_.get_yticks(), a_.get_yticklabels())
                if lo - 1e-12 <= v <= hi + 1e-12]
    ticks = (shown_ticks(ax2) if vec else None,
             shown_ticks(axp) if paths else None)
    write_tikz_deviation(os.path.join(args.out, f"deviation{sfx}.tex"), vec, paths,
                         (lim_v, lim_p), shifted, ticks)

    write_figures_doc(os.path.join(args.out, f"figures{sfx}.tex"),
                      [f"chain{sfx}.tex"] if vec else [],
                      f"deviation{sfx}.tex", args.unbiased)

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
