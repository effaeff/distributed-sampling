"""Sample data under the restriction of a certain distribution for the subsets"""

import os
import random
import math
import misc
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from sklearn.neighbors import NearestNeighbors
from max_rectangle import max_rectangle, plot_grid


def distributed_sampling(data):
    """Method for the actual clustering algorithm"""
    print("length of considered remaining data: {}".format(len(data)))
    nbrs = NearestNeighbors(n_neighbors=2, algorithm='kd_tree', n_jobs=-1).fit(data)
    distances, __ = nbrs.kneighbors(data)
    # Use the maximum distance between a point to his nearest neighbors as cell radius
    cell_radius = np.max(distances[:, 1])

    cellsize = cell_radius / math.sqrt(2) * 2
    nb_cells = int(math.ceil(1 / cellsize))
    grid = np.zeros((nb_cells, nb_cells))

    # Generate cost matrix of potential rectangles through dynamic programming
    grid_awesome = max_rectangle(data, cellsize)
    plot_grid(grid_awesome, primitive=False)

    # Draw samples
    __, axs = plt.subplots(1)
    axs.scatter(data[:, 0], data[:, 1])

    # Count samples in each grid cell
    for point in data:
        grid_x, grid_y = int(math.floor(point[0] / cellsize)), int(math.floor(point[1] / cellsize))
        grid[grid_y, grid_x] += 1

    # Draw grid
    for grid_y in range(nb_cells):
        for grid_x in range(nb_cells):
            origin_x = grid_x * cellsize
            origin_y = grid_y * cellsize
            if grid[grid_y, grid_x] > 0:
                rect = patches.Rectangle(
                    (origin_x, origin_y),
                    cellsize,
                    cellsize,
                    linewidth=1,
                    edgecolor='r',
                    facecolor='none'
                )
                axs.add_patch(rect)

    # plt.figure()
    # plt.scatter(
    #     np.array(grid_awesome[1, 39].max_rect_point_list)[:, 0],
    #     np.array(grid_awesome[1, 39].max_rect_point_list)[:, 1]
    # )
    plt.show()

    quit()

def main():
    """Main Method"""
    random_seed = 1234
    random.seed(random_seed)
    test_size = 0.4
    data_dir = 'W:/Projects/SFB876/Publications/Force_Model/Data/4_features'
    # data_dir = '//ls14-vmnas.cs.tu-dortmund.de/home/Projects/SFB876/Publications/Force_Model/Data/4_features'

    filenames = [
        filename for filename in os.listdir(data_dir)
        if filename.endswith('_features.npy')
    ]
    train_files, val_files = train_test_split(
        filenames,
        test_size=test_size,
        random_state=random_seed
    )
    val_files, __ = train_test_split(val_files, test_size=0.5, random_state=random_seed)

    features = []
    for idx in tqdm(range(len(train_files[:10]))):
        x__ = np.load('{}/{}'.format(data_dir, train_files[idx]))
        features.append(np.unique(x__[np.where(x__[:, 0] > 0)], axis=0))
    features = np.vstack(features)

    scaler = MinMaxScaler()
    features[:, :2] = scaler.fit_transform(features[:, :2])
    print("length of complete feature space: {}".format(len(features)))
    data = features[:, :2].copy()
    distributed_sampling(data)


    ###########################
    ## Max rectangle shizzle ##
    ###########################
    # nb_rows = 10
    # nb_cols = 10
    # # nb_zeros = 10
    # data = np.ones((nb_rows, nb_cols))
    # # for __ in range(nb_zeros):
    # #     data[random.randint(0, nb_cells - 1), random.randint(0, nb_cells - 1)] = 0
    # for row in range(nb_rows):
    #     for col in range(nb_cols):
    #         if row == nb_cols - col - 1:
    #             data[row, col] = 0
    # # data[6, 4] = 0
    # grid = max_rectangle_awesome(data)
    # plot_grid(grid, primitive=False)

if __name__ == "__main__":
    misc.to_local_dir(__file__)
    main()
