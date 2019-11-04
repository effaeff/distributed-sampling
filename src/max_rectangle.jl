module MaxRectangle
export max_rectangle

using Parameters
using StatsBase
using ProgressMeter
using PyCall

plt = pyimport("matplotlib.pyplot")
patches = pyimport("matplotlib.patches")
cmap = pyimport("matplotlib.cm")
colors = pyimport("matplotlib.colors")


@with_kw mutable struct Cell
    dim::Int64
    rect_list::Array{Vector{Int64}} = []
    dist_to_zero::Vector{Int64} = [0 for __ in 1:dim]
    max_rect::Vector{Int64} = [0 for __ in 1:dim]
    color::NTuple{4, Float64} = (1.0, 1.0, 1.0, 0.0)
    point_list::Array{Vector{Float64}} = []
    point_list_rects::Array{Array{Vector{Float64}}} = []
    point_list_max_rect::Array{Vector{Float64}} = []
end

function Base.show(io::IO, object::Cell)
    println(io, "$(object.dist_to_zero)\n$(object.max_rect)")
end

function find_candidates(neighbors, grid, grid_idx, cellsize)
    rect_candidates = Array{Int64}(
        undef,
        sum(neighbor->length(neighbor.rect_list), neighbors),
        length(neighbors)
    )
    rect_candidates_points = []
    candidate_idx = 1
    for (dim, neighbor) in enumerate(neighbors)
        for (rect_idx, rect) in enumerate(neighbor.rect_list)
            local_rect::Vector{Int64} = [
                idx == dim ? rect[idx] + 1 : min(grid[grid_idx].dist_to_zero[idx], rect[idx])
                for idx in 1:length(neighbors)
            ]

            local_points = neighbor.point_list[rect_idx]
            tb_start = CartesianIndex(
                [
                    idx == dim ? grid_idx[idx] : grid_idx[idx] - min(
                        grid[grid_idx].dist_to_zero[idx], rect[idx] + 1
                    )
                    for idx in 1:length(neighbors)
                ]...
            )
            local_points = vcat(
                local_points,
                [grid[tb_idx].point_list for tb_idx in tb_start:grid_idx]
            )

            push!(rect_candidates_points, local_points)
            rect_candidates[candidate_idx] = local_rect
            candidate_idx += 1
    return rect_candidates, rect_candidates_points

function max_rectangle(data::Array{Float64}, cellsize::Float64, dim::Int64)
    nb_cells = convert(Int64, ceil(1 / cellsize))

    # Build grid based on nb_cells and dim
    grid = [
        Cell(dim=dim) for __ in CartesianIndices(Base.OneTo.(Tuple([nb_cells for __ in 1:dim])))
    ]

    # Store points in corresponding grid cells
    for point_idx in 1:size(data, 1)
        grid_idx = CartesianIndex(
            # List comprehension of grid coordinates
            [
                convert(Int64, floor(data[point_idx, dim_idx] / cellsize)) + 1
                for dim_idx in reverse(1:dim)
            ]... # Unpack list components for usage as arguments for constructing CartesianIndex
        )
        push!(grid[grid_idx].point_list, data[point_idx, :])
    end

    global_max = 0
    progress = Progress(nb_cells ^ dim, 1)
    grid_indices = CartesianIndices(grid)
    for grid_idx in grid_indices
        if !isempty(grid[grid_idx].point_list)
            # Admittedly brainfuck-ish
            grid[grid_idx].dist_to_zero = [
                # Construct each element of dist_to_zero independently
                # by using the dist_to_zero of previous cells in each dimension.
                # Note the reversed order of storage of the coordinates in grid,
                # e.g., grid[z, y, x] in three dimensions.
                grid[
                    CartesianIndex(
                        [
                            idx_dim == dim_idx ? max(idx_value - 1, 1) : idx_value
                            for (idx_dim, idx_value) in enumerate(Tuple(grid_idx))
                        ]...
                    )
                ].dist_to_zero[dim - dim_idx + 1] + 1
                for dim_idx in reverse(1:dim)
            ]
            row = grid_idx[1]
            col = grid_idx[2]
            # Check if all of the dimensional components of the index are 1 except for one.
            # Then, max_rect and point_list just accumulate along the dimension which is not 1.
            if count(x->x == 1, Tuple(grid_idx)) >= dim - 1
                # dist_to_zero is already calculated
                grid[grid_idx].max_rect = grid[grid_idx].dist_to_zero
                push!(grid[grid_idx].rect_list, grid[grid_idx].dist_to_zero)
                point_list = grid[grid_idx].point_list
                if grid[grid_idx].dist_to_zero != [1 for __ in 1:dim]
                    # point_list is the accumulation of point_list at current cell and point_list
                    # of previous cell along the dimension whose index component is not 1
                    point_list = vcat(
                        point_list,
                        grid[
                            # max is evaluated for each index component independently.
                            # Due to only being one component differing from 1,
                            # the index is always reduced along the desired dimension.
                            max(grid_idx - oneunit(grid_idx), first(grid_indices))
                        ].point_list_max_rect
                    )
                end
                push!(grid[grid_idx].point_list_rects, point_list)
                grid[grid_idx].point_list_max_rect = point_list
            else
                # If there is more than one dimension differing from 1, the current cell is 
                # surrounded by other cells.
                # The max_rect of the current cell is estimated as a combination of a rect from
                # the rect_list of non-overlapping rects of a neighboring cell
                # and the possibilities, which the inclusion of the current rect offers.
                neighboring_cells = Array{Cell, 1}(undef, dim)
                for neighbor_idx in 1:dim
                    neighboring_cells[neighbor_idx] = grid[
                        # Reduce the dimensional components of grid_idx by 1 one ofter the other
                        # to get neighboring cells. 
                        grid_idx - CartesianIndex(
                            (neighbor_idx == dim_idx ? 1 : 0 for dim_idx in reverse(1:dim))...
                        )
                    ]
                end
                if all(neighbor->isempty(neighbor.rect_list), neighboring_cells)
                    # If rect_list of all neighbors are empty, the current cell is the start
                    # of a new rect.
                    grid[grid_idx].max_rect = [1 for __ in 1:dim]
                    push!(grid[grid_idx].rect_list, [1 for __ in 1:dim])
                    push!(grid[grid_idx].point_list_rects, grid[grid_idx].point_list)
                    grid[grid_idx].point_list_max_rect = grid[grid_idx].point_list
                else
                    max_rect = [0 for __ in 1:dim]
                    point_list_max_rect = []
                    rect_candidates = Array{Int64}(
                        undef,
                        sum(neighbor->length(neighbor.rect_list), neighboring_cells),
                        dim
                    )
                    candidates, candidate_points = find_candidates(
                        neighboring_cells,
                        grid,
                        grid_idx,
                        cellsize
                    )


                    for candidate_idx in 1:size(rect_candidates, 1)
                        candidate = rect_candidates[candidate_idx, :]
                        if prod(candidate) > prod(max_rect)
                            max_rect = candidate
                            point_list_max_rect = rect_candidates_points[candidate_idx]
                        end
                        (
                            isempty(grid[row, col].rect_list) ?
                            append_rect = true :
                            append_rect = false
                        )
                        for rect_idx in 1:size(grid[row, col].rect_list, 1)
                            rect = grid[row, col].rect_list[rect_idx]
                            if rect[1] >= candidate[1] && rect[2] >= candidate[2]
                                append_rect = false
                                break
                            elseif rect[1] <= candidate[1] && rect[2] <= candidate[2]
                                deleteat!(grid[row, col].rect_list, rect_idx)
                                deleteat!(grid[row, col].point_list_rects, rect_idx)
                            end
                            append_rect = true
                        end
                        if append_rect
                            push!(
                                grid[row, col].point_list_rects,
                                rect_candidates_points[candidate_idx]
                            )
                            push!(grid[row, col].rect_list, candidate)
                        end
                    end
                    grid[row, col].max_rect = max_rect
                    grid[row, col].point_list_max_rect = point_list_max_rect
                end
            end
            rect_size = prod(grid[row, col].max_rect)
            if rect_size > global_max
                global_max = rect_size
            end
        end
        next!(progress)
    end
    colormap = cmap.get_cmap("viridis")
    normalize = colors.Normalize(vmin=0, vmax=global_max)
    for grid_idx in CartesianIndices(grid)
        rect_size = prod(grid[grid_idx].max_rect)
        if rect_size != 0
            grid[grid_idx].color = colormap(normalize(rect_size))
        end
    end

    return grid
end

end
