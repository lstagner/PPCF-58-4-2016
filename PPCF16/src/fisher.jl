type FisherFunctional{T,N}
  A::Array{T,2}
  b::Array{T,1}
  Ascale::T
  L::NTuple{N,Array{T,2}}
end

function FisherFunctional{T<:Real,N}(A::AbstractArray{T,2},
                                     b::AbstractArray{T,1},
                                     L::NTuple{N,AbstractArray{T,2}})
  Amax = maximum(A)
	A = A .* (1.0/Amax)
	Ascale = 1.0*(1.0/Amax)
  return FisherFunctional(A, b, Ascale, L)
end

function FisherFunctional{T<:Real,N}(A0::AbstractArray{T,2},
                                     b0::AbstractArray{T,1},
                                     err::AbstractArray{T,1},
                                     L::NTuple{N,AbstractArray{T,2}})
  A, b = scale_by_error(A0,b0,err)
  return FisherFunctional(A, b, L)
end

function minimize{T<:Real,N}(RF::FisherFunctional{T,N},
                             alpha::T;
                             maxiter=50, tol=1e-2,
                             wmax=1.0, verbose=false)
  nr,nc = size(RF.A)
	ATA = RF.A'*RF.A
	ATB = RF.A'*RF.b

  W = eye(nc)
  H = zeros(nc,nc)
  for L in RF.L
    H = H .+ L'*W*L
  end
  x = ((ATA + alpha.*H)\ATB)
  cnt=0
  iteration_tol=0.0
  for i=1:maxiter
    xold = copy(x)
    W = diagm(float(1.0 ./ (max(x,1.0/wmax))))
    H = H*0
    for L in RF.L
      H = H .+ L'*W*L
    end
    x = ((ATA + alpha.*H)\ATB)
    cnt = cnt + 1
    iteration_tol = mean(abs((x .- xold)./xold))
    iteration_tol < tol && break
  end

  rho = norm(RF.A*x .- RF.b)^2
  xi = 0.0
  for L in RF.L
    xi = xi + norm(sqrt(W)*L*x)^2
  end

  Adag = inv(ATA + alpha*H)*RF.A'
  sigma = sqrt(diag(Adag*Adag'))
  
  if verbose
    if cnt != maxiter
      info("Converged in $cnt iterations: tol = $iteration_tol")
    else
      warn("Did not converge in $maxiter iterations: tol = $iteration_tol")
    end
  end

  return RegularizedSolution(RF.Ascale*x,RF.Ascale*sigma,alpha,rho,xi)
end
