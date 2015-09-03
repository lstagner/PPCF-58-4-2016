type TikhonovFunctional{T,N}
  A::Array{T,2}
  b::Array{T,1}
  Ascale::T
  L::NTuple{N,Array{T,2}}
end

function TikhonovFunctional{T<:Real,N}(A::Array{T,2},b::Array{T,1},L::NTuple{N,Array{T,2}})
  Amax = maximum(A)
	A = A .* (1.0/Amax)
	Ascale = 1.0*(1.0/Amax)
  return TikhonovFunctional(A, b, Ascale, L)
end

function TikhonovFunctional{T<:Real,N}(A0::Array{T,2},b0::Array{T,1},err::Array{T,1},L::NTuple{N,Array{T,2}})
  A, b = scale_by_error(A0,b0,err)
  return TikhonovFunctional(A, b, L)
end


function minimize{T<:Real,N}(RF::TikhonovFunctional{T,N},alpha::T; true_sol = Void)
  nr,nc = size(RF.A)
	ATA = RF.A'*RF.A
	ATB = RF.A'*RF.b

  H = zeros(nc,nc)
  for L in RF.L
    H = H .+ L'*L
  end

	x = ((ATA + alpha.*H)\ATB)

  rho = norm(RF.A*x .- RF.b)^2
  xi = 0.0
  for L in RF.L
    xi = xi + norm(L*x)^2
  end

  Adag = inv(ATA + alpha*H)*RF.A'
  sigma = sqrt(diag(Adag*Adag'))

  if true_sol != Void
    x_true = vec(true_sol)/RF.Ascale
    reg_err = (x_true .- Adag*(RF.A*x_true))
  else
    reg_err = zeros(length(x))
  end

	return RegularizedSolution(RF.Ascale*x,RF.Ascale*sigma,RF.Ascale*reg_err,alpha,rho,xi)
end
