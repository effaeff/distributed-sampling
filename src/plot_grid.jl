module PlotGrid
export plot_grid
"""Method for grid plotting"""

using PyCall

plt = pyimport("matplotlib.pyplot")
patches = pyimport("matplotlib.patches")


function plot_grid(grid, width, height, primitive)
    """
    Method for plotting grids.
    Simple numpy arrays can also be plotted of primitive is True
    """
    nb_rows = size(grid, 1)
    nb_cols = size(grid, 2)
    __, axs = plt.subplots(1, figsize=(40, 40))
    for row in 1:nb_rows
        for col in 1:nb_cols
            facecolor = "none"
            if !primitive
                facecolor = grid[row, col].color
            end
            r_x = col * width
            r_y = row * height
            axs.add_artist(
                patches.Rectangle(
                    (r_x, r_y),
                    width,
                    height,
                    facecolor=facecolor,
                    edgecolor="black"
                )
            )
            axs.annotate(
                string(grid[row, col]),
                (r_x + width / 2.0, r_y + height / 2.0),
                ha="center",
                va="center"
            )
            # Debug
            axs.annotate(
                string(length(grid[row, col].rect_list)),
                (r_x + width / 2.0, r_y)
            )
        end
    end
    axs.set_xlim((0, nb_cols * width))
    axs.set_ylim((0, nb_rows * height))
    plt.tick_params(
        axis="both",
        which="both",
        bottom=false,
        top=false,
        left=false,
        labelbottom=false,
        labelleft=false
    )
    plt.tight_layout(pad=0)
    plt.savefig("grid_full.png", dpi=600)
end

end
