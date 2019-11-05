"""Methods for grid plotting"""
module PlotGrid
export plot_grid, plot_rasterized_samples

using PyCall

plt = pyimport("matplotlib.pyplot")
patches = pyimport("matplotlib.patches")


function plot_rasterized_samples(data, grid, cellsize)
    """
    Method for plotting grid ontop of samples
    """
    __, axs = plt.subplots(1)
    plt.scatter(data[:, 1], data[:, 2], s=0.2)
    for grid_y in 1:size(grid, 1)
        for grid_x in 1:size(grid, 2)
            origin_x = (grid_x - 1) * cellsize
            origin_y = (grid_y - 1) * cellsize
            isempty(grid[grid_y, grid_x].point_list) ? color = "g" : color = "r"
                rect = patches.Rectangle(
                (origin_x, origin_y),
                cellsize,
                cellsize,
                linewidth=1,
                edgecolor=color,
                facecolor="none"
            )
            axs.add_patch(rect)
            if grid_x == 1
                axs.annotate(
                    string(grid_y),
                    (origin_x - cellsize / 2.0, origin_y + cellsize / 2.0),
                    ha="center",
                    va="center",
                    fontsize=5
                )
            end
            if grid_y == 1
                axs.annotate(
                    string(grid_x),
                    (origin_x + cellsize / 2.0, origin_y - cellsize / 2.0),
                    ha="center",
                    va="center",
                    fontsize=5
                )
            end
        end
    end
    plt.savefig("rasterized_samples.png", dpi=600)
    plt.close()
end

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
    plt.savefig("grid.png", dpi=600)
    plt.close()
end

end
