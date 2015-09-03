type FisherFunctional{T,N}
  A::Array{T,2}
  b::Array{T,1}
  Ascale::T
  L::NTuple{N,Array{T,2}}
end

function FisherFunctional{T<:Real,N}(A::Array{T,2},b::Array{T,1},L::NTuple{N,Array{T,2}})
  Amax = maximum(A)
	A = A .* (1.0/Amax)
	Ascale = 1.0*(1.0/Amax)
  return FisherFunctional(A, b, Ascale, L)
end

function FisherFunctional{T<:Real,N}(A0::Array{T,2},b0::Array{T,1},err::Array{T,1},L::NTuple{N,Array{T,2}})
  A, b = scale_by_error(A0,b0,err)
  return FisherFunctional(A, b, L)
end

function minimize{T<:Real,N}(RF::FisherFunctional{T,N},alpha::T; maxiter=50, tol=1e-2, wmax=1.0,
  verbose=false,true_sol = Void)
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
  
  if true_sol != Void
    x_true = vec(true_sol)/RF.Ascale
    reg_err = x_true .- Adag*(RF.A*x_true)
  else
    reg_err = zeros(length(x))
  end

  if verbose
    if cnt != maxiter
      info("Converged in $cnt iterations: tol = $iteration_tol")
    else
      warn("Did not converge in $maxiter iterations: tol = $iteration_tol")
    end
  end

  return RegularizedSolution(RF.Ascale*x,RF.Ascale*sigma,RF.Ascale*reg_err,alpha,rho,xi)
end
