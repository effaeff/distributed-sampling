"""Sample data under the restriction of a certain distribution for the subsets"""

import time
import os
import random
import math
import misc
from tqdm import tqdm
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from sklearn.neighbors import NearestNeighbors
from max_rectangle import max_rectangle
from plot_grid import plot_grid

from matplotlib import pyplot as plt

def distributed_sampling(data):
    """Method for the actual clustering algorithm"""
    print("length of considered remaining data: {}".format(len(data)))
    nbrs = NearestNeighbors(n_neighbors=2, algorithm='kd_tree', n_jobs=-1).fit(data)
    distances, __ = nbrs.kneighbors(data)
    # Use the maximum distance between a point to his nearest neighbors as cell radius
    cell_radius = np.max(distances[:, 1])

    cellsize = cell_radius / math.sqrt(2) * 2

    # Generate cost matrix of potential rectangles through dynamic programming
    start = time.time()
    grid_awesome = max_rectangle(data, cellsize)
    end = time.time()
    print("Runtime: {} seconds".format(end - start))

    # plot_grid(grid_awesome, primitive=False)

def main():
    """Main Method"""
    random_seed = 1234
    random.seed(random_seed)
    test_size = 0.4
    # data_dir = 'W:/Projects/SFB876/Publications/Force_Model/Data/4_features'
    data_dir = 'C:/Data/Workspace/distributed_sampling/data'

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
    for idx in tqdm(range(len(train_files))):
        x__ = np.load('{}/{}'.format(data_dir, train_files[idx]))
        features.append(np.unique(x__[np.where(x__[:, 0] > 0)], axis=0))
    features = np.vstack(features)

    scaler = MinMaxScaler()
    features[:, :2] = scaler.fit_transform(features[:, :2])
    data = features[:, :2].copy()

    plt.scatter(data[:, 0], data[:, 1])
    plt.savefig("test.png")

    # distributed_sampling(data)

    ###################################
    ## Debug (currently not working) ##
    ###################################
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
