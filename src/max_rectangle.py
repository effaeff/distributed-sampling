"""
Methods to determine the maximum rectangle for each cell of a grid,
indluding the set of points lying in the subgrid
"""

import math
import numpy as np
import matplotlib

from tqdm import tqdm


class Cell:
    """Class for one cell of the regarded grid"""
    def __init__(self, max_rect=[0, 0], dist_to_zero=[0, 0]):
        self.rect_list = []
        self.dist_to_zero = dist_to_zero
        self.max_rect = max_rect
        self.color = 'none'
        self.point_list = []
        self.point_list_rects = []
        self.max_rect_point_list = []

    def __str__(self):
        return '{}\n{}'.format(self.dist_to_zero, self.max_rect)

def max_rectangle(data, cellsize):
    """Calculate max rectangle for each cell of grid through multiple new brilliant ideas"""
    nb_cells = int(math.ceil(1 / cellsize))

    nb_rows = nb_cells
    nb_cols = nb_cells

    grid = np.reshape(
        np.array([Cell() for __ in range(nb_rows * nb_cols)]),
        (nb_rows, nb_cols)
    )

    for point in data:
        grid_x = int(math.floor(point[0] / cellsize))
        grid_y = int(math.floor(point[1] / cellsize))
        grid[grid_y, grid_x].point_list.append(point)

    global_max = 0
    progress = tqdm(total=nb_cols * nb_rows)
    for row in range(nb_rows):
        for col in range(nb_cols):
            if grid[row, col].point_list:
                grid[row, col].dist_to_zero = [
                    grid[row, col - 1 if col - 1 >= 0 else 0].dist_to_zero[0] + 1,
                    grid[row - 1 if row - 1 >= 0 else 0, col].dist_to_zero[1] + 1
                ]
                if row == 0 or col == 0:
                    grid[row, col].max_rect = grid[row, col].dist_to_zero.copy()
                    grid[row, col].rect_list.append(grid[row, col].dist_to_zero.copy())
                    dist_to_zero = grid[row, col].dist_to_zero
                    point_list = grid[row, col].point_list.copy()
                    if grid[row, col].dist_to_zero != [1, 1]:
                        point_list += grid[
                            row - 1 if dist_to_zero[1] > 1 else row,
                            col - 1 if dist_to_zero[0] > 1 else col
                        ].max_rect_point_list.copy()
                    grid[row, col].point_list_rects.append(point_list)
                    grid[row, col].max_rect_point_list = point_list.copy()
                else:
                    lower = grid[row - 1, col]
                    left = grid[row, col - 1]
                    if not left.rect_list and not lower.rect_list:
                        grid[row, col].max_rect = [1, 1]
                        grid[row, col].rect_list.append([1, 1])
                        grid[row, col].point_list_rects.append(grid[row, col].point_list.copy())
                        grid[row, col].max_rect_point_list = grid[row, col].point_list.copy()
                    else:
                        max_rect = [0, 0]
                        max_rect_point_list = []
                        rect_candidates = np.empty((len(left.rect_list) + len(lower.rect_list), 2))
                        rect_candidates_points = []
                        for lower_idx, lower_rect in enumerate(lower.rect_list):
                            local_rect = [
                                min(grid[row, col].dist_to_zero[0], lower_rect[0]),
                                lower_rect[1] + 1
                            ]
                            local_points = lower.point_list_rects[lower_idx].copy()
                            traceback_start = col - min(
                                grid[row, col].dist_to_zero[0], lower_rect[0]
                            ) + 1
                            traceback_end = col + 1
                            for traceback_idx in range(traceback_start, traceback_end):
                                local_points += grid[row, traceback_idx].point_list.copy()
                            lowest_x = (col - local_rect[0] + 1) * cellsize
                            local_points = np.array(local_points)
                            local_points = local_points[
                                np.where(local_points[:, 0] > lowest_x)
                            ].tolist()
                            rect_candidates_points.append(local_points)
                            rect_candidates[lower_idx] = np.array(local_rect)
                        for left_idx, left_rect in enumerate(left.rect_list):
                            local_rect = [
                                left_rect[0] + 1,
                                min(grid[row, col].dist_to_zero[1], left_rect[1])
                            ]
                            local_points = left.point_list_rects[left_idx].copy()
                            traceback_start = row - min(
                                grid[row, col].dist_to_zero[1], left_rect[1]
                            ) + 1
                            traceback_end = row + 1
                            for traceback_idx in range(traceback_start, traceback_end):
                                local_points += grid[traceback_idx, col].point_list.copy()
                            lowest_y = (row - local_rect[1] + 1) * cellsize
                            local_points = np.array(local_points)
                            local_points = local_points[
                                np.where(local_points[:, 1] > lowest_y)
                            ].tolist()
                            rect_candidates_points.append(local_points)
                            rect_candidates[len(lower.rect_list) + left_idx] = np.array(local_rect)
                        for candidate_idx, candidate in enumerate(rect_candidates):
                            if np.prod(candidate) > np.prod(max_rect):
                                max_rect = candidate
                                max_rect_point_list = rect_candidates_points[candidate_idx]
                            append_rect = True if not grid[row, col].rect_list else False
                            for rect_idx, rect in enumerate(grid[row, col].rect_list):
                                if rect[0] >= candidate[0] and rect[1] >= candidate[1]:
                                    append_rect = False
                                    break
                                elif rect[0] <= candidate[0] and rect[1] <= candidate[1]:
                                    del grid[row, col].rect_list[rect_idx]
                                    del grid[row, col].point_list_rects[rect_idx]
                                append_rect = True
                            if append_rect:
                                grid[row, col].point_list_rects.append(
                                    rect_candidates_points[candidate_idx]
                                )
                                grid[row, col].rect_list.append(candidate.tolist())
                        grid[row, col].max_rect = max_rect
                        grid[row, col].max_rect_point_list = max_rect_point_list
                rect_size = grid[row, col].max_rect[0] * grid[row, col].max_rect[1]
                if rect_size > global_max:
                    global_max = rect_size
            progress.update(1)
    cmap = matplotlib.cm.get_cmap('viridis')
    normalize = matplotlib.colors.Normalize(vmin=0, vmax=global_max)
    for row in range(nb_rows):
        for col in range(nb_cols):
            rect_size = grid[row, col].max_rect[0] * grid[row, col].max_rect[1]
            if rect_size != 0:
                grid[row, col].color = cmap(normalize(rect_size))
    return grid
