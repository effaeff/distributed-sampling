module MaxRectangle
export max_rectangle

import Base.show

using Parameters
using ProgressMeter
using PyCall

plt = pyimport("matplotlib.pyplot")
patches = pyimport("matplotlib.patches")
cmap = pyimport("matplotlib.cm")
colors = pyimport("matplotlib.colors")


@with_kw mutable struct Cell
    rect_list::Array{Vector{Int64}} = []
    dist_to_zero::Vector{Int64} = [0, 0]
    max_rect::Vector{Int64} = [0, 0]
    color::NTuple{4, Float64} = (1.0, 1.0, 1.0, 0.0)
    point_list::Array{Vector{Float64}} = []
    point_list_rects::Array{Array{Vector{Float64}}} = []
    point_list_max_rect::Array{Vector{Float64}} = []
end

function show(io::IO, object::Cell)
    println(io, "$(object.dist_to_zero)\n$(object.max_rect)")
end

function find_lower_candidates(lower, lower_idx, grid, row, col, cellsize)
    lower_rect = lower.rect_list[lower_idx]
    local_rect = [
        min(grid[row, col].dist_to_zero[1], lower_rect[1]),
        lower_rect[2] + 1
    ]
    local_points = lower.point_list_rects[lower_idx]
    traceback_start = col - min(
        grid[row, col].dist_to_zero[1], lower_rect[1]
    ) + 1
    traceback_end = col
    for traceback_idx in traceback_start:traceback_end
        local_points = vcat(
            local_points,
            grid[row, traceback_idx].point_list
        )
    end
    lowest_x = (col - local_rect[1]) * cellsize
    local_points = [
        local_points[idx] for idx in 1:size(local_points, 1)
        if local_points[idx][1] > lowest_x
    ]
    return local_rect, local_points
end

function find_left_candidates(left, left_idx, grid, row, col, cellsize)
    left_rect = left.rect_list[left_idx]
    local_rect = [
        left_rect[1] + 1,
        min(grid[row, col].dist_to_zero[2], left_rect[2])
    ]
    local_points = left.point_list_rects[left_idx]
    traceback_start = row - min(
        grid[row, col].dist_to_zero[2], left_rect[2]
    ) + 1
    traceback_end = row
    for traceback_idx in traceback_start:traceback_end
        local_points = vcat(
            local_points,
            grid[traceback_idx, col].point_list
        )
    end
    lowest_y = (row - local_rect[2]) * cellsize
    local_points = [
        local_points[idx] for idx in 1:size(local_points, 1)
        if local_points[idx][2] > lowest_y
    ]
    return local_rect, local_points
end

function max_rectangle(data::Array{Float64, 2}, cellsize::Float64)
    nb_cells = convert(Int64, ceil(1 / cellsize))
    nb_rows = nb_cells
    nb_cols = nb_cells
    grid = [Cell() for __ in 1:nb_rows, __ in 1:nb_cols]

    for point_idx in 1:size(data, 1)
        grid_x = convert(Int64, floor(data[point_idx, 1] / cellsize)) + 1
        grid_y = convert(Int64, floor(data[point_idx, 2] / cellsize)) + 1
        push!(grid[grid_y, grid_x].point_list, data[point_idx, :])
    end

    global_max = 0
    progress = Progress(nb_cols * nb_cols, 1)
    for row in 1:nb_rows
        for col in 1:nb_cols
            if !isempty(grid[row, col].point_list)
                grid[row, col].dist_to_zero = [
                    grid[row, col - 1 >= 1 ? col - 1 : 1].dist_to_zero[1] + 1,
                    grid[row - 1 >= 1 ? row - 1 : 1, col].dist_to_zero[2] + 1
                ]
                if row == 1 || col == 1
                    grid[row, col].max_rect = grid[row, col].dist_to_zero
                    push!(grid[row, col].rect_list, grid[row, col].dist_to_zero)
                    dist_to_zero = grid[row, col].dist_to_zero
                    point_list = grid[row, col].point_list
                    if grid[row, col].dist_to_zero != [1, 1]
                        point_list = vcat(
                            point_list,
                            grid[
                                row - 1 > 0 ? row - 1 : row,
                                col - 1 > 0 ? col - 1 : col
                            ].point_list_max_rect
                            
                        )
                    end
                    push!(grid[row, col].point_list_rects, point_list)
                    grid[row, col].point_list_max_rect = point_list
                else
                    lower = grid[row - 1, col]
                    left = grid[row, col - 1]
                    if isempty(left.rect_list) && isempty(lower.rect_list)
                        grid[row, col].max_rect = [1, 1]
                        push!(grid[row, col].rect_list, [1, 1])
                        push!(grid[row, col].point_list_rects, grid[row, col].point_list)
                        grid[row, col].point_list_max_rect = grid[row, col].point_list
                    else
                        max_rect = [0, 0]
                        point_list_max_rect = []
                        rect_candidates = Array{Int64}(
                            undef,
                            length(left.rect_list) + length(lower.rect_list),
                            2
                        )
                        rect_candidates_points = []
                        for lower_idx in 1:size(lower.rect_list, 1)
                            local_rect, local_points = find_lower_candidates(
                                lower, lower_idx, grid, row, col, cellsize
                            )
                            push!(rect_candidates_points, local_points)
                            rect_candidates[lower_idx, :] = local_rect
                        end
                        for left_idx in 1:size(left.rect_list, 1)
                            local_rect, local_points = find_left_candidates(
                                left, left_idx, grid, row, col, cellsize
                            )
                            push!(rect_candidates_points, local_points)
                            rect_candidates[length(lower.rect_list) + left_idx, :] = local_rect
                        end
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
    end
    colormap = cmap.get_cmap("viridis")
    normalize = colors.Normalize(vmin=0, vmax=global_max)
    for row in 1:nb_rows
        for col in 1:nb_cols
            rect_size = prod(grid[row, col].max_rect)
            if rect_size != 0
                grid[row, col].color = colormap(normalize(rect_size))
            end
        end
    end

    return grid
end

end
