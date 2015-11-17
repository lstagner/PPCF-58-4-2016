type CubicSpline{T} # Cubic Hermite Spline
  n::Int        # Number of knots
  x::Array{T,1} # x position of knots
  y::Array{T,1} # y position of knots
  m::Array{T,1} # Tangent at knots
end

function CubicSpline{T<:Real}(x::AbstractArray{T,1},y::AbstractArray{T,1})
  issorted(x) || throw(ArgumentError("x points must be in ascending order"))
  nx = length(x)
  m = zeros(nx)
  m[1] = (y[2] - y[1])/(x[2]-x[1])
  for i=2:nx-1
    m[i] = 0.5*((y[i+1]-y[i])/(x[i+1]-x[i]) + (y[i]-y[i-1])/(x[i]-x[i-1]))
  end
  m[nx] = (y[nx]-y[nx-1])/(x[nx]-x[nx-1])
  return CubicSpline(nx,x,y,m)
end

function MonotoneCubicSpline{T<:Real}(x::AbstractArray{T,1},y::AbstractArray{T,1})
  issorted(x) || throw(ArgumentError("x points must be in ascending order"))
  nx = length(x)
  m = zeros(nx)
  m[1] = (y[2] - y[1])/(x[2]-x[1])

  for i=2:nx-1
    hi = x[i+1]-x[i]
    hi_1 = x[i]-x[i-1]
    di = (y[i+1]-y[i])/hi
    di_1 = (y[i]-y[i-1])/hi_1
    m[i] = sign(di) != sign(di_1) ? 0.0 : 3*(hi_1+hi)*(((2*hi+hi_1)/di_1)+((hi+2*hi_1)/di))^(-1.0)
  end

  m[nx] = (y[nx]-y[nx-1])/(x[nx]-x[nx-1])

  return CubicSpline(nx,x,y,m)
end

function interpolate{T<:Real}(S::CubicSpline{T},x::Union{T,AbstractArray{T,1}})
  xrange = extrema(S.x)
  if any((x .< xrange[1]) | (x .> xrange[2]))
    throw(ArgumentError("Outside of Range"))
  end
  nx = length(x)
  yout = zeros(nx)
  for i = 1:nx
    xr = searchsorted(S.x,x[i])
    i1 = xr.stop
    i2 = xr.start

    if i1 != i2
      dx = (S.x[i2] - S.x[i1])
      t = (x[i] - S.x[i1])/dx
      h00 = 2*t^3 - 3*t^2 + 1
      h10 = t^3 - 2*t^2 + t
      h01 = -2*t^3 + 3*t^2
      h11 = t^3 - t^2
      yout[i] = h00*S.y[i1] + h10*dx*S.m[i1] + h01*S.y[i2] + h11*dx*S.m[i2]
    else
      yout[i] = S.y[i1]
    end

  end
  return yout
end
