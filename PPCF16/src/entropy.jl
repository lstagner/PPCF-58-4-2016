type EntropyFunctional{T}
  A::Array{T,2}
  b::Array{T,1}
  Ascale::T
  m::Model
  alpha::Array{T,1}
  d::Array{T,1}
end

function EntropyFunctional{T<:Real}(A::AbstractArray{T,2},
                                    b::AbstractArray{T,1};
                                    Ascale=1e4, alpha_start=1e0,
                                    d = Void)

  A = Ascale*A
  nr,nc = size(A)
  alpha = [alpha_start]

  # Create default image
  if d == Void
    d = fill(Ascale*(1e-6)*mean(b)/mean(A),nc)
  end

  m = Model(solver=IpoptSolver(tol=1e-5,max_iter=500,print_level=0))

  @defVar(m, x[1:nc] >= eps())
  @defVar(m, z[1:nr])

  @addConstraint(m, eq[i=1:nr], z[i] == sum{A[i,j]*x[j],j=1:nc} - b[i])

  @setNLObjective(m, Min, -alpha[1]*sum{x[i] - d[i] - x[i]*log(x[i]/d[i]),i=1:nc} + 0.5*sum{z[i]^2,i=1:nr})

  for i = 1:nc
      setValue(x[i],d[i])
  end
  for i = 1:nr
      setValue(z[i], dot(vec(A[i,:]),d))
  end

  return EntropyFunctional(A,b,Ascale,m,alpha,d,blur)
end

function EntropyFunctional{T<:Real}(A0::AbstractArray{T,2},
                                    b0::AbstractArray{T,1},
                                    err::AbstractArray{T,1}; kwargs...)
  A, b = scale_by_error(A0, b0, err)
  return EntropyFunctional(A,b; kwargs...)
end

function minimize{T<:Real}(RF::EntropyFunctional{T}, alpha::T)
  nr,nc = size(RF.A)
  A = RF.A*RF.blur
  RF.alpha[1] = alpha
  status=solve(RF.m)

  xx = getVar(RF.m,:x)
  x = getValue(xx)

  rho = norm(A*x .- RF.b)^2
  d = RF.d
  xi = max(-sum(x .- d .- x.*log(x./d)),eps())

  cov = inv(A'*A + alpha*diagm(1.0./x))
  sigma = sqrt(diag(cov))

  return RegularizedSolution(RF.Ascale*x, RF.Ascale*sigma,alpha, rho, xi)
end

function maxEnt(A0::AbstractArray{Float64,2},
                b0::AbstractArray{Float64},
                err::AbstractArray{Float64};
                Ascale=1.0, rChisq=1.0,
                d=Void,system=:Ipopt)

    nr,nc = size(A0)

    # Scale matrices by error
    A,b = scale_by_error(A0,b0,err)

    # Create default image
    if d == Void
        d = fill((1e-6)*mean(b)/mean(A),nc)
        #d = ones(nc)./map.nelements
    end

    A = Ascale*A
    if system == :Ipopt
        f, rchsq, status = maxEntJuMP(A,b,d,float(rChisq*nr))
    else
        error("Unknown Solver System")
    end

    return f.*Ascale, rchsq, status

end

function maxEntJuMP(A::AbstractArray{Float64,2},
                    b::AbstractArray{Float64},
                    d::AbstractArray{Float64},
                    c::Float64)

    nr,nc = size(A)

    m = Model(solver=IpoptSolver(tol=1e-5,max_iter=500,print_level=0))

    @defVar(m, f[1:nc] >= eps())
    @defVar(m, z[1:nr])

    @addConstraint(m, eq[i=1:nr], z[i] == sum{A[i,j]*f[j],j=1:nc} - b[i])
    chiSq = [c]
    @addNLConstraint(m, sum{z[i]^2,i=1:nr} <= chiSq[1])

    @setNLObjective(m, Max, sum{f[i] - d[i] - f[i]*log(f[i]/d[i]),i=1:nc})

    for i = 1:nc
        setValue(f[i],d[i])
    end
    for i = 1:nr
        setValue(z[i], dot(vec(A[i,:]),d))
    end

    status=solve(m)
    count = 0
    while status != :Optimal && count < 5
        count+=1
        fp = getValue(f)
        chiSq[1] = ceil((norm(A*fp-b)^2.0)/nr,1)*nr
        info("Increasing reduced chi-squared to "*string(round(chiSq[1]/nr,2)))
        status = solve(m)
    end
    status == :Optimal && info("Optimal Solution Found!")
    sol = getValue(f)
    return sol, (norm(A*sol-b)^2.0)/nr, status

end

type ConvexEntropyFunctional{T}
  A::Array{T,2}
  b::Array{T,1}
  Ascale::T
  alpha::Array{T,1}
  d::Array{T,1}
  m::Convex.Problem
  x::Convex.Variable
end

function ConvexEntropyFunctional{T<:Real}(A::AbstractArray{T,2},
                                          b::AbstractArray{T,1};
                                          Ascale = 1e4,
                                          alpha_start = 1e0,
                                          d = Void)
  
  A = Ascale*A
  nr,nc = size(A)
  alpha = [alpha_start]

  # Create default image
  if d == Void
    d = fill(Ascale*(1e-6)*mean(b)/mean(A),nc)
  end

  x = Convex.Variable(nc)
  alpha = [alpha_start]

  problem = Convex.minimize((-alpha[1]*(sum(x .- d) + entropy(x./d)) + 0.5*sumsquares(A*x - b)))

  return ConvexEntropyFunctional(A,b,Ascale,alpha,d,problem,x)
end

function ConvexEntropyFunctional{T<:Real}(A0::AbstractArray{T,2},
                                          b0::AbstractArray{T,1},
                                          err::AbstractArray{T,1}; kwargs...)

  A, b = scale_by_error(A0, b0, err)
  return ConvexEntropyFunctional(A,b; kwargs...)
end

function minimize{T<:Real}(RF::ConvexEntropyFunctional{T}, alpha::T)
  nr,nc = size(RF.A)
  A = RF.A
  RF.alpha[1] = alpha
  ws = typeof(RF.x.value) != Void
  solve!(RF.m,warmstart=ws)

  status = RF.m.status
  x = RF.x.value[:,1]

  rho = norm(A*x .- RF.b)^2
  d = RF.d
  xi = max(-sum(x .- d .- x.*log(x./d)),eps())

  cov = inv(A'*A + alpha*diagm(1.0./x))
  sigma = sqrt(diag(cov))

  return RegularizedSolution(RF.Ascale*x, RF.Ascale*sigma,alpha, rho, xi)
end
