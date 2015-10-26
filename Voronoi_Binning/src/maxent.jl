using JuMP
using Ipopt

function maxEnt(A0::Array{Float64,2}, b0::Array{Float64}, err::Array{Float64}; Ascale=1.0, rChisq=1.0, d=nothing,system=:Ipopt)
    nr,nc = size(A0)

    # Set a minimum error level
    err[err .< 0.02*mean(b0)] = 0.02*mean(b0)

    # Scale matrices by error
    A = zeros(nr,nc)
    b = zeros(nr)

    for i=1:nr
        A[i,:] = (Ascale/err[i]).*A0[i,:]
        b[i] = b0[i]/err[i]
    end

    # Create default image
    if d == nothing
        d = fill((1e-6)*mean(b)/mean(A),nc)
        #d = ones(nc)./map.nelements
    end

    if system == :Ipopt
        f, rchsq, status = maxEntJuMP(A,b,d,float(rChisq*nr))
    # elseif system == :SCS
    #     f, rchsq, status = maxEntConvex(A,b,d,float(rChisq*nr))
    else
        error("Unknown Solver System")
    end

    return f.*Ascale, rchsq, status

end

function maxEnt(A0::Array{Float64,2}, b0::Array{Float64}, err::Array{Float64},map::BinMapping; Ascale=1.0, rChisq=1.0, d=nothing,system=:Ipopt)
   nr,nc = size(A0)

   # Scale matrices by error
   A = zeros(nr,nc)
   b = zeros(nr)

   for i=1:nr
       A[i,:] = (Ascale/err[i]).*A0[i,:]
       b[i] = b0[i]/err[i]
   end

   # Create default image
   if d == nothing
       d = fill((1e-6)*mean(b)/mean(A),nc) #./map.nelements
   end

   if system == :Ipopt
       f, rchsq, status = maxEntJuMP(A,b,d,float(rChisq*nr))
  #  elseif system == :SCS
  #      f, rchsq, status = maxEntConvex(A,b,d,float(rChisq*nr))
   else
       error("Unknown Solver System")
   end

   return BinnedArray(map,f.*Ascale), rchsq, status

end

function maxEntJuMP(A::Array{Float64,2}, b::Array{Float64}, d::Array{Float64}, c::Float64)

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
        sol = getValue(f)
        fp = sol.innerArray
        chiSq[1] = ceil((norm(A*fp-b)^2.0)/nr,1)*nr
        info("Increasing reduced chi-squared to "*string(round(chiSq[1]/nr,2)))
        status = solve(m)
    end
    status == :Optimal && info("Optimal Solution Found!")
    sol = getValue(f)
    return sol.innerArray, (norm(A*sol.innerArray-b)^2.0)/nr, status

end

# Lets keep this commented out until bugs in Convex.jl are fixed
#function maxEntConvex(A::Array{Float64,2}, b::Array{Float64}, d::Array{Float64}, c::Float64)
#
#    nr, nc = size(A)
#
#    f = Convex.Variable(nc)
#
#    problem = maximize(sum( f .- d .+ entropy(f) .+ f.*log(d)))
#
#    problem.constraints += [sum_squares(A*f - b) <= c, f > eps()]
#
#    solve!(problem,SCSSolver())
#
#    return f.value, (norm(A*f.value-b)^2.0)/nr, problem.status
#
#end
