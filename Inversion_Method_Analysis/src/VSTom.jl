module VSTom

using PyPlot
using Grid

export read_ncdf
export gradEP
export resize_transfer_matrix
export bilinear
export make_stacked_plot,make_ep_plot
export write_to_file


include("netcdf.jl")
include("gradients.jl")
include("utils.jl")
include("plotting.jl")

end
