module RegularizedSystems
using PyPlot
using JuMP
using Ipopt

export minimize
export TSVDFunctional, TikhonovFunctional, FisherFunctional, EntropyFunctional
export RegularizedSolution
export interpolate, PolyharmonicSpline, CubicSpline, MonotoneCubicSpline
export discrete_lcurve, lcurve
export scale_by_error

type RegularizedSolution{T<:FloatingPoint}
  x::Array{T,1}
  sigma::Array{T,1}
  reg_err::Array{T,1}
  hyperparameter::Union(T,Int)
  rho::T
  xi::T
end

function scale_by_error(A0::Array{Float64,2},b0::Array{Float64,1},err::Array{Float64,1})
  	nr,nc = size(A0)

  	# Set a minimum error level
  	err[err .< 0.02*mean(b0)] = 0.02*mean(b0)

  	# Scale matrices by error
  	A = zeros(nr,nc)
  	b = zeros(nr)

  	for i=1:nr
  		A[i,:] = (1.0/err[i]).*A0[i,:]
  		b[i] = b0[i]/err[i]
  	end
    return A, b
end

include("tsvd.jl")
include("tikhonov.jl")
include("fisher.jl")
include("entropy.jl")
include("cubic_spline.jl")
include("polyharmonic_spline.jl")
include("lcurve.jl")

end
