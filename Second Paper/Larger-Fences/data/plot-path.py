# plot-path.py
import sys
import os

def plot_path(filepath, n):
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Polygon
    except ImportError:
        print("Error: Matplotlib is not installed")
        return

    print(f"Reading path coordinates from {filepath}...")
    xs, ys = [], []
    with open(filepath, "r") as f:
        for line in f:
            parts = line.strip().split(",")
            if len(parts) == 2:
                xs.append(float(parts[0]))
                ys.append(float(parts[1]))

    num_points = len(xs)
    print(f"Loaded {num_points:,} coordinates.")

    # Downsample path if it exceeds limits to prevent slow render times
    max_render_points = 5000000
    if num_points > max_render_points:
        step = num_points // max_render_points
        xs_plot = xs[::step] + [xs[-1]]
        ys_plot = ys[::step] + [ys[-1]]
    else:
        xs_plot = xs
        ys_plot = ys

    fig, ax = plt.subplots(figsize=(10, 10), dpi=300)
    
    # Fill interior
    polygon = Polygon(list(zip(xs_plot, ys_plot)), closed=True, 
                      facecolor='#f4f7f6', edgecolor='none', zorder=1)
    ax.add_patch(polygon)

    # Calculate exact segment index boundaries to prevent any gaps from rounding
    num_segments = min(256, len(xs_plot) - 1)
    indices = [int(i * (len(xs_plot) - 1) / num_segments) for i in range(num_segments + 1)]
    
    cmap = plt.get_cmap('viridis')
    lw = 1.5 if num_points < 10000 else 0.8

    for i in range(num_segments):
        start_idx = indices[i]
        end_idx = indices[i + 1] + 1  # +1 ensures inclusive slicing for continuous path overlap
        ax.plot(xs_plot[start_idx:end_idx], ys_plot[start_idx:end_idx], 
                color=cmap(i / num_segments), linewidth=lw, zorder=2)

    ax.set_aspect('equal', adjustable='box')
    ax.axis('off')

    plt.title(f"{n}-omino Fence ({num_points-1:,} vectors)", 
              fontsize=14, color='#2c3e50', pad=10)

    out_img = f"data/fence-{n}.png"
    plt.tight_layout()
    plt.savefig(out_img, bbox_inches='tight', dpi=300, facecolor='white')
    plt.close()
    print(f"Saved visualization to {out_img}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 plot-path.py <path_file_path> <n>")
        sys.exit(1)
    plot_path(sys.argv[1], int(sys.argv[2]))
