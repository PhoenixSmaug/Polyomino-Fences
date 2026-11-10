#!/usr/bin/env python3
"""Histograms of the polyomino diameter vectors by angle.

Reads the same data/vec*.csv files as chain_plot.py.  A line n,x,y,d is a
diameter vector (x,y), 0 <= y <= x, of a sample with weight 1/d: a sample
with d diameter vectors gives d consecutive lines, which share the sample
equally.  The vectors go into L equal buckets by their angle to the x-axis
over 0..45 degrees, and the program draws
  a) the relative frequency of each bucket (sum of the weights 1/d divided
     by the number of samples),
  b) the weighted average length |(x,y)| of the vectors in each bucket,
in one histogram, on screen and as TikZ (pgfplots) file hist.tex: the
average length in blue in the background with its scale, the frequency in
green in front, without a scale, scaled so that its maximum is --fg-scale
(default 0.85) times the height of the length maximum.
Angles of exactly 45 degrees go into the last bucket.  It also prints the
distribution of the d values, the empty and the smallest buckets, and the
slope report for p/q with q <= --max-q, all with weights 1/d.

Usage:  python3 angle_histogram.py [data_dir] [-L 90] [-o out_dir]
                                 [--max-vec N] [--fg-scale F] [--max-q Q]
                                 [--no-show]
"""

import argparse
import bisect
import glob
import os
import sys
from fractions import Fraction
from math import gcd

import numpy as np

from chain_plot import VEC_PATTERN, n_samples, read_vectors


def bucket_stats(x, y, d, L):
    """Weighted bucket statistics; every line has weight 1/d."""
    w = 1.0 / d
    angle = np.degrees(np.arctan2(y, x))
    length = np.hypot(x, y)
    b = np.clip((angle / 45.0 * L).astype(int), 0, L - 1)
    count = np.bincount(b, weights=w, minlength=L)   # samples per bucket
    total = np.bincount(b, weights=w * length, minlength=L)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = total / count                     # NaN for empty buckets
    return count / w.sum(), mean, count


def report_d(d):
    """Distribution of the d values: a sample with d diameter vectors gives
    d lines, so it is counted once per d lines."""
    vals, lines = np.unique(d, return_counts=True)
    samples = lines / vals
    total = samples.sum()
    print(f"\ndistribution of d (number of diameter vectors per sample), "
          f"{total:g} samples:")
    print("       d    samples   rel. freq     lines")
    for v, sm, ln in zip(vals, samples, lines):
        print(f"  {v:6d}  {sm:9g}  {sm / total:10.6f}  {ln:8d}")


def report_buckets(edges, count, freq, mean, n_smallest=5):
    """Print all empty buckets and the n_smallest smallest non-empty ones
    (ties in bucket order)."""
    empty = np.flatnonzero(count == 0)   # count: samples (weights 1/d)
    nonempty = np.flatnonzero(count > 0)
    smallest = nonempty[np.argsort(count[nonempty], kind="stable")][:n_smallest]
    print(f"\nempty buckets: {len(empty)}"
          + ("" if len(empty) else " (none)"))
    if len(empty):
        print("  bucket   angle range")
        for i in empty:
            print(f"  {i:6d}   {edges[i]:6.2f}-{edges[i + 1]:5.2f} deg")
    print(f"\nthe {len(smallest)} smallest non-empty buckets:")
    print("  bucket   angle range        samples  rel. freq  avg length")
    for i in smallest:
        print(f"  {i:6d}   {edges[i]:6.2f}-{edges[i + 1]:5.2f} deg  "
              f"{count[i]:8.2f}  {freq[i]:9.5f}  {mean[i]:10.2f}")


def fmt_slope(s):
    """A slope as reduced fraction, decimal value and arctangent in degrees."""
    if s is None:
        return "(none)"
    deg = np.degrees(np.arctan2(s.numerator, s.denominator))
    return f"{str(s):>7} = {float(s):.6f}, {deg:8.4f} deg"


def report_slopes(x, y, d, max_q):
    """For each reduced p/q with 0 <= p <= q <= max_q: weighted count
    (samples, every line 1/d) and relative frequency of vectors with exactly
    this slope y/x, and the nearest occurring slopes above and below it.
    Slopes are compared exactly (integer arithmetic)."""
    g = np.gcd(x, y)
    pairs, inv = np.unique(np.stack([y // g, x // g], axis=1), axis=0,
                           return_inverse=True)
    cnt = np.bincount(inv.ravel(), weights=1.0 / d)
    occ = sorted(zip((Fraction(int(p), int(q)) for p, q in pairs),
                     cnt.tolist()))
    slopes = [s for s, _ in occ]
    counts = dict(occ)
    total = n_samples(d)
    targets = sorted({Fraction(p, q) for q in range(1, max_q + 1)
                      for p in range(0, q + 1) if gcd(p, q) == 1})
    print(f"\nslopes p/q with 0 <= p <= q <= {max_q} "
          f"({len(slopes)} distinct slopes occur among {total:g} samples):")

    for t in targets:
        c = counts.get(t, 0)
        i = bisect.bisect_right(slopes, t)
        above = slopes[i] if i < len(slopes) else None
        j = bisect.bisect_left(slopes, t) - 1
        below = slopes[j] if j >= 0 else None
        print(f"  {fmt_slope(t)}: samples {c:8.2f}, rel. freq {c / total:.5f}")
        print(f"      next larger:  {fmt_slope(above)}")
        print(f"      next smaller: {fmt_slope(below)}")


def write_tikz_hist(path, edges, length, freq_shown, ymax):
    """Two overlaid bar histograms with the same fixed vertical
    range [0, ymax] in length units: average length (blue, background, with
    its scale) and the rescaled relative frequency (green, foreground, no
    scale).  Empty buckets (NaN) are drawn as gaps."""
    # scale only axis: both axes get the same plot area although only the
    # first one has tick labels.
    common = ("scale only axis, width=0.85\\linewidth, "
              "height=0.4\\linewidth, ybar, bar shift=0pt, "
              f"bar width={edges[1] - edges[0]:.6g}, "
              f"xmin=0, xmax=45, ymin=0, ymax={ymax:.6g}")

    def bars(f, vals, style):
        f.write(f"\\addplot[{style}] coordinates {{\n")
        # bars at the bucket centres; the unitless bar width is in axis
        # units (pgfplots compat 1.7 or later)
        for e0, e1, v in zip(edges[:-1], edges[1:], vals):
            f.write(f"({(e0 + e1) / 2:.4f},{0.0 if np.isnan(v) else v:.6g})\n")
        f.write("};\n")

    with open(path, "w") as f:
        f.write("% Average vector length (blue, left scale) and relative "
                "frequency (green,\n% rescaled, no scale) per angle bucket, "
                "0 to 45 degrees.\n% Needs \\usepackage{pgfplots} "
                "(\\pgfplotsset{compat=1.17} or later).\n")
        f.write("\\begin{tikzpicture}\n")
        f.write(f"\\begin{{axis}}[{common}, xtick={{0,5,...,45}}, "
                "xticklabel={$\\pgfmathprintnumber{\\tick}^\\circ$}, "
                "scaled y ticks=false, "
                "xlabel={$\\alpha$}, ylabel={average length}]\n")
        bars(f, length, "fill=blue!35, draw=blue!70!black, thin")
        f.write("\\end{axis}\n")
        f.write(f"\\begin{{axis}}[{common}, hide axis]\n")
        bars(f, freq_shown, "fill=green!55!black, draw=green!30!black, thin")
        f.write("\\end{axis}\n\\end{tikzpicture}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("data_dir", nargs="?", default="data")
    ap.add_argument("-L", type=int, default=90, help="number of buckets")
    ap.add_argument("-o", "--out", default=".", help="output directory")
    ap.add_argument("--max-vec", type=int, default=None, metavar="N",
                    help="stop reading the vec files after N samples "
                         "(each line counts 1/d)")
    ap.add_argument("--fg-scale", type=float, default=0.85, metavar="F",
                    help="height of the frequency maximum relative to the "
                         "length maximum (default 0.85)")
    ap.add_argument("--max-q", type=int, default=6, metavar="Q",
                    help="report the slopes p/q with 0 <= p <= q <= Q "
                         "(default 6)")
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

    files = sorted(glob.glob(os.path.join(args.data_dir, VEC_PATTERN)))
    if not files:
        sys.exit(f"no {VEC_PATTERN} files in {args.data_dir}")
    x, y, d = read_vectors(files, args.max_vec)
    L = args.L
    freq, mean, count = bucket_stats(x, y, d, L)
    edges = np.linspace(0, 45, L + 1)
    print(f"{n_samples(d):g} samples in {L} buckets of {45 / L:g} deg; "
          f"{(count == 0).sum()} empty; samples per bucket "
          f"{count.min():g}..{count.max():g}", file=sys.stderr)
    report_d(d)
    report_buckets(edges, count, freq, mean)
    report_slopes(x, y, d, args.max_q)

    # Common vertical range in length units, fixed explicitly so that no
    # automatic scaling can change the ratio.  The frequency is mapped so
    # that its maximum would match the length maximum, then times --fg-scale.
    lmax = np.nanmax(mean)
    freq_shown = freq * (lmax / freq.max()) * args.fg_scale
    ymax = 1.05 * max(lmax, freq_shown.max())

    os.makedirs(args.out, exist_ok=True)
    write_tikz_hist(os.path.join(args.out, "hist.tex"), edges, mean,
                    freq_shown, ymax)

    fig, ax = plt.subplots(figsize=(8, 4))
    w = 45 / L
    ax.bar(edges[:-1], np.nan_to_num(mean), width=w, align="edge",
           color="#9ab3ea", edgecolor="#1f3f8f", linewidth=0.4, zorder=1,
           label="average length")
    # Frequency on the same axis (no second scale), drawn in front.
    ax.bar(edges[:-1], freq_shown, width=w, align="edge",
           color="#2e8b3e", edgecolor="#174d20", linewidth=0.4, zorder=2,
           label="relative frequency")
    ax.set_ylim(0, ymax)
    ax.set_xlim(0, 45)
    ax.set_xticks(range(0, 46, 5))
    ax.xaxis.set_major_formatter("{x:.0f}°")
    ax.set_xlabel("α")
    ax.set_ylabel("average length")
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "hist.png"), dpi=150)

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
