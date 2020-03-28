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
    """
    Routine to find candidates for the maximum rectangle for a cell
    based on rectangles of neighboring cells.
    """
    rect_candidates = Array{Int64, 2}(
        undef,
        sum(neighbor->length(neighbor.rect_list), neighbors),
        length(neighbors)
    )
    rect_candidates_points = []
    candidate_idx = 1
    # Neighbors correspond to neighboring dimensions in decreasing order,
    # e.g. z, y, x in the three dimensional case.
    for (i_neighbor, neighbor) in enumerate(neighbors)
        # i_neighbor iterates grid indices
        # dim iterates dimensions
        dim = length(neighbors) - i_neighbor + 1
        for (rect_idx, rect) in enumerate(neighbor.rect_list)
            # Construct candidate based on neighbor and current cell
            local_rect::Vector{Int64} = [
                (
                    i_dim == dim ?
                    # At the dimension of the neighboring rectangle, the size just increments by 1
                    rect[i_dim] + 1 :
                    # Otherwise, we have to check how far the candidate can be spanned
                    # in the direction of the corresponding dimension
                    min(grid[grid_idx].dist_to_zero[i_dim], rect[i_dim])
                )
                # Dimensional iteration
                for i_dim in 1:length(neighbors)
            ]
            # Points of candidate include all points of the regarded neighboring rectangle
            local_points = neighbor.point_list_rects[rect_idx]
            # Starting index of traceback procedure to collect all additional points
            tb_start = CartesianIndex(
                [
                    (
                        i_index == i_neighbor ?
                        # At the neighboring index, the traceback index equals the grid index
                        grid_idx[i_index] :
                        # The other dimensions have to be iterated.
                        # Starting from the grid index, we go back as far as possible,
                        # limited by dist_to_zero and the neighboring rectangle
                        grid_idx[i_index] - min(
                            grid[grid_idx].dist_to_zero[length(neighbors) - i_index + 1],
                            rect[length(neighbors) - i_index + 1]
                        ) + 1
                    )
                    # Grid index iteration
                    for i_index in 1:length(neighbors)
                ]...
            )
            local_points = vcat(
                local_points,
                [grid[tb_idx].point_list for tb_idx in tb_start:grid_idx]...
            )

            push!(rect_candidates_points, local_points)
            rect_candidates[candidate_idx, :] = local_rect
            candidate_idx += 1
        end
    end
    return rect_candidates, rect_candidates_points
end

function max_rectangle(data::Array{Float64}, cellsize::Float64, dim::Int64)
    """
    Multidimensional implementation of algorithm to find the maximum rectangle in a grid.
    The data is expected to be normalized.
    """
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
            ]... # Splat list components for usage as arguments for constructing CartesianIndex
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
                        max(
                            grid_idx - CartesianIndex(
                                (neighbor_idx == dim_idx ? 1 : 0 for dim_idx in 1:dim)...
                            ),
                            CartesianIndex((1 for __ in 1:dim)...)
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
                    # Find all rectangle candidates using all neighboring cells
                    candidates, candidate_points = find_candidates(
                        neighboring_cells,
                        grid,
                        grid_idx,
                        cellsize
                    )

                    for candidate_idx in 1:size(candidates, 1)
                        candidate = candidates[candidate_idx, :]
                        if prod(candidate) > prod(max_rect)
                            max_rect = candidate
                            point_list_max_rect = candidate_points[candidate_idx]
                        end
                        # Only rectangles, which are not completely overlapped by other
                        # rectangles should be remembered.
                        # So the overlapping between the candidate and the rectangles of
                        # the current cell has to be checked.
                        (
                            # If rect_list is empty, the candidate can be appended to it,
                            # since it will be the only member.
                            isempty(grid[grid_idx].rect_list) ?
                            append_rect = true :
                            append_rect = false
                        )
                        not_overlapped = ones(length(grid[grid_idx].rect_list))
                        for rect_idx in 1:size(grid[grid_idx].rect_list, 1)
                            rect = grid[grid_idx].rect_list[rect_idx]
                            # If a lager rectangle if found,
                            # no other members of rect_list have to be checked
                            if all(rect .>= candidate)
                                append_rect = false
                                break
                            # If the candidate is larger than a member,
                            # the member has to be removed from rect_list
                            elseif all(rect .<= candidate)
                                not_overlapped[rect_idx] = 0
                            end
                            append_rect = true
                        end
                        if any(not_overlapped .!= 1)
                            grid[grid_idx].rect_list = [
                                grid[grid_idx].rect_list[rect_idx]
                                for rect_idx in 1:size(not_overlapped, 1)
                                if not_overlapped[rect_idx] == 1
                            ]
                            grid[grid_idx].point_list_rects = [
                                grid[grid_idx].point_list_rects[rect_idx]
                                for rect_idx in 1:size(not_overlapped, 1)
                                if not_overlapped[rect_idx] == 1
                            ]
                        end
                        if append_rect
                            push!(
                                grid[grid_idx].point_list_rects,
                                candidate_points[candidate_idx]
                            )
                            push!(grid[grid_idx].rect_list, candidate)
                        end
                    end
                    grid[grid_idx].max_rect = max_rect
                    grid[grid_idx].point_list_max_rect = point_list_max_rect
                end
            end
            # Adjust the global maximal product of rectangle components for the colorscale limit
            rect_size = prod(grid[grid_idx].max_rect)
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
