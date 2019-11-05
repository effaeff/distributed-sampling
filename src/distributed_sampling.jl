include("max_rectangle.jl")
include("plot_grid.jl")
using .MaxRectangle
using .PlotGrid

using MLDataUtils
using ProgressBars
using PyCall

np = pyimport("numpy")
model_selection = pyimport("sklearn.model_selection")
neighbors = pyimport("sklearn.neighbors")
preprocessing = pyimport("sklearn.preprocessing")

plt = pyimport("matplotlib.pyplot")


function distributed_sampling(data::Array{Float64}, dim::Int64)
    println("Length of considered remaining data: $(size(data, 1))")

    nbrs = neighbors.NearestNeighbors(n_neighbors=2, algorithm="kd_tree", n_jobs=-1).fit(data)
    distances, __ = nbrs.kneighbors(data)
    # Use the maximum distance between a point to his nearest neighbors as cell radius
    cell_radius = maximum(distances[:, 2])
    cellsize = cell_radius / sqrt(2) * 2

    @time grid = MaxRectangle.max_rectangle(data, cellsize, dim);
    PlotGrid.plot_grid(grid, 1, 1, false)
    PlotGrid.plot_rasterized_samples(data, grid, cellsize)
end

function main()
    random_seed::Int64 = 1234
    test_size::Float64 = 0.4
    data_dir::String = "C:/Data/Workspace/distributed_sampling/data"
    # data_dir::String = "W:/Projects/SFB876/Publications/Force_Model/Data/4_features"

    filenames = filter(filename->occursin("_features.npy", filename), readdir(data_dir))
    train_files, val_files = model_selection.train_test_split(
        filenames,
        test_size=test_size,
        random_state=random_seed
    )
    val_files, __ = model_selection.train_test_split(
        val_files,
        test_size=0.5,
        random_state=random_seed
    )

    features::Array{Array{Float64, 2}} = []
    for (idx, __) in ProgressBar(enumerate(train_files))
        x__ = np.load("$data_dir/$(train_files[idx])")[:, 1:2]
        push!(features, unique(x__[x__[:, 1] .> 0, :], dims=1))
    end
    data = vcat(features...)

    scaler = preprocessing.MinMaxScaler()
    data = scaler.fit_transform(data)
    dim::Int64 = size(data, 2)

    distributed_sampling(data, dim)
end

main()
