"""Method for grid plotting"""

from matplotlib import patches, pyplot as plt


def plot_grid(grid, width=1, height=1, primitive=True):
    """
    Method for plotting grids.
    Simple numpy arrays can also be plotted of primitive is True
    """
    nb_rows = len(grid)
    nb_cols = len(grid[0])
    __, axs = plt.subplots(1, figsize=(40, 40))
    for row in range(nb_rows):
        for col in range(nb_cols):
            facecolor = 'none'
            if not primitive:
                facecolor = grid[row, col].color
            r_x = col * width
            r_y = row * height
            axs.add_artist(
                patches.Rectangle(
                    (r_x, r_y),
                    width,
                    height,
                    facecolor=facecolor,
                    edgecolor='black'
                )
            )
            axs.annotate(
                str(grid[row, col]),
                (r_x + width / 2.0, r_y + height / 2.0),
                ha='center',
                va='center'
            )
            # Debug
            axs.annotate(
                str(len(grid[row, col].rect_list)),
                (r_x + width / 2.0, r_y)
            )
    axs.set_xlim((0, nb_cols * width))
    axs.set_ylim((0, nb_rows * height))
    plt.tick_params(
        axis='both',
        which='both',
        bottom=False,
        top=False,
        left=False,
        labelbottom=False,
        labelleft=False
    )
    plt.tight_layout(pad=0)
    plt.savefig("grid_full_python.png", dpi=600)
