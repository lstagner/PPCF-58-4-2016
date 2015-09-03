type EntropyFunctional{T}
  A::Array{T,2}
  b::Array{T,1}
  Ascale::T
  m::Model
  alpha::Array{T,1}
  d::Array{T,1}
  blur::Array{T,2}
end

function EntropyFunctional{T<:Real}(A::Array{T,2},b::Array{T,1}; Ascale=1e4, alpha_start=1e0,d =
  Void,blur = Void)

  A = Ascale*A
  nr,nc = size(A)
  alpha = [alpha_start]

  if blur == Void
    blur = eye(nc)
  end
  Ab = A*blur

  # Create default image
  if d == Void
    d = fill((1e-6)*mean(b)/mean(Ab),nc)
  end

  m = Model(solver=IpoptSolver(tol=1e-5,max_iter=500,print_level=0))

  @defVar(m, x[1:nc] >= eps())
  @defVar(m, z[1:nr])

  @addConstraint(m, eq[i=1:nr], z[i] == sum{Ab[i,j]*x[j],j=1:nc} - b[i])

  @setNLObjective(m, Min, -alpha[1]*sum{x[i] - d[i] - x[i]*log(x[i]/d[i]),i=1:nc} + 0.5*sum{z[i]^2,i=1:nr})

  for i = 1:nc
      setValue(x[i],d[i])
  end
  for i = 1:nr
      setValue(z[i], dot(vec(Ab[i,:]),d))
  end

  return EntropyFunctional(A,b,Ascale,m,alpha,d,blur)
end

function EntropyFunctional{T<:Real}(A0::Array{T,2},b0::Array{T,1},err::Array{T,1}; kwargs...)
  A, b = scale_by_error(A0, b0, err)
  return EntropyFunctional(A,b; kwargs...)
end

function minimize{T<:Real}(RF::EntropyFunctional{T}, alpha::T; true_sol = Void)
  nr,nc = size(RF.A)
  A = RF.A*RF.blur
  RF.alpha[1] = alpha
  status=solve(RF.m)

  xx = getVar(RF.m,:x)
  sol = getValue(xx)
  x = sol.innerArray

  rho = norm(A*x .- RF.b)^2
  d = RF.d
  xi = max(-sum(x .- d .- x.*log(x./d)),eps())

  cov = inv(A'*A + alpha*diagm(1.0./x))
  sigma = sqrt(diag(cov))

  if true_sol != Void
    x_true = vec(true_sol)
    A = RF.A/RF.Ascale
    tsol = minimize(EntropyFunctional(A,A*x_true,Ascale=RF.Ascale,blur=RF.blur),alpha)
    reg_err = (x_true .- RF.blur*tsol.x)/RF.Ascale
  else
    reg_err = zeros(length(x))
  end

  return RegularizedSolution(RF.Ascale*x, RF.Ascale*sigma,RF.Ascale*reg_err,alpha, rho, xi)
end
#
# function maxEnt(A::Array{Float64,2}, b::Array{Float64}, d::Array{Float64}, c::Float64)
#
#     nr,nc = size(A)
#
#     m = Model(solver=IpoptSolver(tol=1e-5,max_iter=500,print_level=0))
#
#     @defVar(m, f[1:nc] >= eps())
#     @defVar(m, z[1:nr])
#
#     @addConstraint(m, eq[i=1:nr], z[i] == sum{A[i,j]*f[j],j=1:nc} - b[i])
#     chiSq = [c]
#     @addNLConstraint(m, sum{z[i]^2,i=1:nr} <= chiSq[1])
#
#     @setNLObjective(m, Max, sum{f[i] - d[i] - f[i]*log(f[i]/d[i]),i=1:nc})
#
#     for i = 1:nc
#         setValue(f[i],d[i])
#     end
#     for i = 1:nr
#         setValue(z[i], dot(vec(A[i,:]),d))
#     end
#
#     status=solve(m)
#     count = 0
#     while status != :Optimal && count < 5
#         count+=1
#         sol = getValue(f)
#         fp = sol.innerArray
#         chiSq[1] = ceil((norm(A*fp-b)^2.0)/nr,1)*nr
#         info("Increasing reduced chi-squared to "*string(round(chiSq[1]/nr,2)))
#         status = solve(m)
#     end
#     status == :Optimal && info("Optimal Solution Found!")
#     sol = getValue(f)
#     return sol.innerArray, (norm(A*sol.innerArray-b)^2.0)/nr, status
#
# end
