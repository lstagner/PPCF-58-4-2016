module VSTom

using PyPlot
using Grid

export read_ncdf
export gradEP
export resize_transfer_matrix
export bilinear

include("netcdf.jl")
include("gradients.jl")
include("utils.jl")

end
