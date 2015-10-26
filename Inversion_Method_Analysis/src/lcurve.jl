function circumcircle_curvature{T<:Real}(x::AbstractArray{T,1},y::AbstractArray{T,1})
  # The curvature at a point P[i] is defined as 1/R where R is the radius of the osculating
  # circle at P[i]. The osculating circle at P[i] can be approximated by the circumcircle
  # of the points P[i-1],P[i],P[i+1]. The inverse of the radius of this circumcircle is then
  # an approximation of the curvature at P[i].
  nx = length(x)
  kappa = zeros(nx)
  # Endpoints are defined to have zero curvature
  for i=2:nx-1
    P = [x[i-1:i+1] y[i-1:i+1]]
    dis = zeros(3)

    for j=1:3
      dis[j] = norm(P[j,:])^2
    end
    a = det([P ones(3)])
    if abs(a) > eps()
      d = -det([[dis P[:,2]] ones(3)])
      e = det([[dis P[:,1]] ones(3)])
      f = -det([dis P])

      x0 = -d/(2*a)
      y0 = -e/(2*a)
      R  = sqrt(max((d)^2 + (e)^2 - 4*a*f,0.0))/(2*abs(a))
      kappa[i] = sign(y0 - y[i])/R # if center is above curve then positive curvature
    else
      kappa[i] = 0.0
    end
  end

  return kappa
end

function discrete_lcurve(RF, K;
                         nseeds = 21,
                         kwargs = Dict{Symbol,Any}(),
                         doplot=false, filename = Void)

  k_step = (K[end]-K[1])/nseeds
  k_values = unique(int(K[1]:k_step:K[end]))
  k_set = IntSet(k_values)
  num = length(k_values)
  rho = zeros(num)
  xi = zeros(num)

  # Calculate rho,xi for k values
  for (i,k) in enumerate(k_values)
    sol = minimize(RF, k; kwargs...)
    rho[i] = sol.rho
    xi[i]  = sol.xi
  end

  # rho needs to be in ascending order
  if !issorted(rho)
    w = sortperm(rho)
    rho = rho[w]
    xi = xi[w]
    k_values = k_values[w]
  end

  lrho = log10(rho)
  lxi = log10(xi)

  while true
    i1 = MonotoneCubicSpline(lrho,lxi)
    i2 = PolyharmonicSpline(1,[lrho lxi],Float64(k_values))

    curvature = circumcircle_curvature(lrho,lxi)
    cmax_ind = indmax(curvature)

    k_new = Int[]
    if cmax_ind != 1
      lrho_left = lrho[cmax_ind]-0.5*abs(lrho[cmax_ind] - lrho[cmax_ind-1])
      lxi_left = interpolate(i1,lrho_left)[1]
      k = interpolate(i2,lrho_left,lxi_left)[1]
      k1, k2 = floor(Int64,k), ceil(Int64,k)
      k1_in_set = issubset(k1,k_set)
      k2_in_set = issubset(k2,k_set)

      if (!k1_in_set $ !k2_in_set)
        !k1_in_set ? k_left = k1 : k_left = k2
        push!(k_new,k_left)
      elseif !k1_in_set & !k2_in_set
        k_left = round(Int64,k)
        push!(k_new,k_left)
      end

    end

    if cmax_ind != length(lrho)
      lrho_right = lrho[cmax_ind] + 0.5*abs(lrho[cmax_ind+1] - lrho[cmax_ind])
      lxi_right = interpolate(i1,lrho_right)[1]
      k = interpolate(i2,lrho_right,lxi_right)[1]
      k1, k2 = floor(Int64,k), ceil(Int64,k)
      k1_in_set = issubset(k1,k_set)
      k2_in_set = issubset(k2,k_set)

      if (!k1_in_set $ !k2_in_set)
        !k1_in_set ? k_right = k1 : k_right = k2
        push!(k_new,k_right)
      elseif !k1_in_set & !k2_in_set
        k_right = round(Int64,k)
        push!(k_new,k_right)
      end
    end

    length(k_new) == 0 && break

    for k in k_new
      sol = minimize(RF,k)
      push!(lrho,log10(sol.rho))
      push!(lxi,log10(sol.xi))
      push!(k_values,k)
      union!(k_set,k)
    end

    if !issorted(lrho)
      w = sortperm(lrho)
      lrho = lrho[w]
      lxi = lxi[w]
      k_values = k_values[w]
    end
  end

  i1 = MonotoneCubicSpline(lrho,lxi)
  i2 = PolyharmonicSpline(1,[lrho lxi],float(k_values))

  curvature = circumcircle_curvature(lrho,lxi)
  cmax_ind = indmax(curvature)

  if doplot
    top_k = k_values[cmax_ind]
    fig,ax = plt[:subplots](ncols=2)
    fig[:set_size_inches](12,5)
    lr = linspace(extrema(lrho)...,200)
    lx = interpolate(i1,lr)
    ax[1][:plot](exp10(lr),exp10(lx),"k")
    ax[1][:plot](exp10(lrho),exp10(lxi),"bo",markersize=3)
    ax[1][:plot](exp10(lrho[cmax_ind]),exp10(lxi[cmax_ind]),"ro",label = "K = $(top_k)")
    #ax[1][:set_aspect]("equal")
    ax[1][:set_xscale]("log")
    ax[1][:set_yscale]("log")
    ax[1][:set_xlabel](L"\chi^2(x_{\alpha})",fontsize=18,fontname="cmr10")
    ax[1][:set_ylabel](L"\mathcal{R}(x_{\alpha})",fontsize=18,fontname="cmr10")
    #ax[1][:set_title]("L Curve",fontsize=20,fontname="cmr10")
    ax[1][:legend](loc=3,numpoints=1)
    ax[2][:plot](exp10(lrho),curvature,"k")
    ax[2][:plot](exp10(lrho),curvature,"bo")
    ax[2][:plot](exp10(lrho[cmax_ind]),curvature[cmax_ind],"ro",label=string(top_k))
    ax[2][:set_xlabel](L"\chi^2(x_{\alpha})",fontsize=18,fontname="cmr10")
    ax[2][:set_ylabel](L"\kappa",fontsize=18,fontname="cmr10")
    ax[2][:set_xscale]("log")
    #ax[2][:set_title]("Curvature",fontsize=20,fontname="cmr10")
    ax[2][:legend](numpoints=1,loc=3)
    fig[:tight_layout]()
    if filename != Void
        fig[:savefig](filename,bbox_inches="tight",dpi=1200)
    end
  end
  return k_values[cmax_ind]
end

function lcurve(RF;
                log_alpha_range = (-5,20),
                nseeds=5, maxiter = 50,
                log_tol = 1.0e-5,
                doplot=false, filename=Void,
                verbose=false, debug_plots=false,
                kwargs = Dict{Symbol,Any}())

  # Calculate initial points
  alpha = logspace(log_alpha_range...,nseeds)
  xi = zeros(nseeds)
  rho = zeros(nseeds)
  for i=1:nseeds
    sol = minimize(RF,alpha[i]; kwargs...)
    rho[i] = sol.rho
    xi[i]  = sol.xi
  end

  # Do interpolation in log space
  lal  = log10(alpha)
  lxi  = log10(xi)
  lrho = log10(rho)

  # Make sure everything is sorted
  if !issorted(lrho)
    w = sortperm(lrho)
    lal = lal[w]
    lrho = lrho[w]
    lxi = lxi[w]
  end
  n = 1000
  alpha0 = 0.0
  xi0 = 0.0
  rho0 = 0.0
  num_lrho = length(lrho)
  alpha_vec = Float64[]
  count=0
  for i=1:maxiter
    # Interpolate points
    i1 = MonotoneCubicSpline(lrho,lxi)
    i2 = PolyharmonicSpline(2,hcat(lrho,lxi),lal)

    curvature = circumcircle_curvature(lrho,lxi)

    ind = indmax(curvature)
    rho0 = exp10(lrho[ind])
    xi0 = exp10(lxi[ind])
    alpha0 = exp10(lal[ind])
    if all(curvature .<= 0.0)
      warn("L Curve is concave. Increase nseeds or log_alpha_range.")
      break
    end

    # Debugging plots
    if debug_plots
      fig,ax = plt[:subplots](ncols=2)
      ax[1][:plot](lrho,lxi,"k")
      ax[1][:plot](lrho,lxi,"bo",markersize=3)
      ax[1][:plot](lrho[ind],lxi[ind],"ro")
      ax[1][:set_aspect]("equal")
      ax[1][:set_xlabel](L"log10(\chi^2(x_{\alpha}))")
      ax[1][:set_ylabel](L"log10(\mathcal{R}(x_{\alpha}))")
      ax[2][:plot](lrho,curvature,"k")
      ax[2][:plot](lrho[ind],curvature[ind],"ro")
      ax[2][:set_yscale]("symlog")
      ax[2][:set_xlabel](L"log10(\chi^2(x_{\alpha}))")
      ax[2][:set_ylabel](L"\kappa")
    end

    if ind != 1
      lrho1 = lrho[ind] - 0.5*abs(lrho[ind] - lrho[ind-1])
    end

    if ind != num_lrho
      lrho2 = lrho[ind] + 0.5*abs(lrho[ind+1] - lrho[ind])
    end

    if ind != 1
      lxi1 = interpolate(i1,lrho1)[1]
      alpha1 = exp10(interpolate(i2,lrho1,lxi1))[1]
      sol = minimize(RF,alpha1; kwargs...)
      rho1 = sol.rho
      xi1 = sol.xi
      push!(lrho,log10(rho1))
      push!(lxi,log10(xi1))
      push!(lal,log10(alpha1))
    end

    if ind != num_lrho
      lxi2 = interpolate(i1,lrho2)[1]
      alpha2 = exp10(interpolate(i2,lrho2,lxi2))[1]
      sol = minimize(RF,alpha2; kwargs...)
      rho2 = sol.rho
      xi2 = sol.xi
      push!(lrho,log10(rho2))
      push!(lxi,log10(xi2))
      push!(lal,log10(alpha2))
    end

    # Make sure everything is sorted
    if !issorted(lrho)
      w = sortperm(lrho)
      lal = lal[w]
      lrho = lrho[w]
      lxi = lxi[w]
    end

    if verbose
      println("Iteration: $(i)")
      println(@sprintf(" α = %3.4e",alpha0))
      println(@sprintf(" ρ = %3.4e",rho0))
      println(@sprintf(" ξ = %3.4e",xi0))
      println(@sprintf(" Δlog(ρ) = %3.4e", log10(rho2/rho1)))
    end

    # Increase counters
    num_lrho = length(lrho)
    count = count+1
    push!(alpha_vec,alpha0)

    # Check if converged in rho
    abs(log10(rho2/rho1)) < log_tol && break
  end

  curvature = circumcircle_curvature(lrho,lxi)
  ind = indmax(curvature)
  rho0 = exp10(lrho[ind])
  xi0 = exp10(lxi[ind])
  alpha0 = exp10(lal[ind])
  push!(alpha_vec,alpha0)
  count = count+1

  # Make sure everything is sorted
  if !issorted(lrho)
    w = sortperm(lrho)
    lal = lal[w]
    lrho = lrho[w]
    lxi = lxi[w]
  end

  if doplot
    i1 = MonotoneCubicSpline(lrho,lxi)
    lrho_p = linspace(minimum(lrho),maximum(lrho),n)
    lxi_p = interpolate(i1,lrho_p)
    fig,ax = plt[:subplots](ncols=2)
    fig[:set_size_inches](12,5)
    ax[1][:plot](exp10(lrho_p),exp10(lxi_p),"k")
    ax[1][:plot](exp10(lrho),exp10(lxi),"bo",markersize=3)
    ax[1][:plot](rho0,xi0,"ro",label=latexstring(@sprintf("\\alpha = %3.4e",alpha0)))
    ax[1][:set_xscale]("log")
    ax[1][:set_yscale]("log")
    ax[1][:set_xlabel](L"\chi^2(x_{\alpha})",fontsize=18)
    ax[1][:set_ylabel](L"\mathcal{R}(x_{\alpha})",fontsize=18)
    #ax[1][:set_title]("L Curve",fontname="cmr10",fontsize=20)
    ax[1][:legend](loc=3,numpoints=1)
    ax[2][:plot](collect(1:count),alpha_vec,"k")
    ax[2][:set_ylabel](L"\alpha",fontsize=18)
    ax[2][:set_xlabel]("Iteration",fontname="cmr10",fontsize=18)
    #ax[2][:set_title]("Alpha Convergence",fontname="cmr10",fontsize=18)
    fig[:tight_layout]()
    if filename != Void
      fig[:savefig](filename,bbox_inches="tight",dpi=1200)
    end
    #plt[:close]("all")
  end

  return alpha0
end
