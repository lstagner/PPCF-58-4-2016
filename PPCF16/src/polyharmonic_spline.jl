type PolyharmonicSpline
  dim::Int64
  order::Int64
  coeff::Vector{Float64}
  centers::Array{Float64,2}
  error::Float64
end

function polyharmonicK(r,K)
  if iseven(K)
    return r < 1 ? (r.^(K-1))*log(r.^r) : (r^K)*log(r)
  else
    return r^K
  end
end

function PolyharmonicSpline(K::Int64,
                            centers::AbstractArray{Float64,2},
                            values::AbstractArray{Float64};
                            s = 0.0)
  m,n = size(centers)
  m != length(values) && throw(DimensionMismatch())

  M = zeros(m,m)
  N = zeros(m,n+1)

  for i=1:m
    N[i,1] = 1
    N[i,2:end] = centers[i,:]
    for j=1:m
      M[i,j] = polyharmonicK(norm(centers[i,:] .- centers[j,:]),K)
    end
  end
  M = M .+ s.*eye(m)
  L = [[M N]; [N' zeros(n+1,n+1)]]

  w = pinv(L)*[values;zeros(n+1)]

  ivalues = zeros(m)
  for i=1:m
    tmp = 0.0
    for j=1:m
      tmp = tmp + w[j]*polyharmonicK(norm(centers[i,:] .- centers[j,:]),K)
    end
    tmp = tmp + w[m+1]
    for j=2:n+1
      tmp = tmp + w[m+j]*centers[i,j-1]
    end
    ivalues[i] = tmp
  end
  error = norm(values .- ivalues)^2

  return PolyharmonicSpline(n,K,w,centers,error)
end

function PolyharmonicSpline(K::Int64,
                            centers::Vector{Float64},
                            values::Vector{Float64};
                            s = 0.0)

  PolyharmonicSpline(K,centers'',values,s=s)
end

function interpolate(S::PolyharmonicSpline,x::AbstractArray{Float64,2})
  m,n = size(x)

  n != S.dim && throw(DimensionMismatch("$m != $(S.dim)"))

  interpolates = zeros(m)
  for i=1:m
    tmp = 0.0
    l = length(S.coeff)-(n+1)
    for j=1:l
      tmp = tmp + S.coeff[j]*polyharmonicK(norm(x[i,:] .- S.centers[j,:]),S.order)
    end
    tmp = tmp + S.coeff[l+1]
    for j=2:n+1
      tmp = tmp + S.coeff[l+j]*x[i,j-1]
    end
    interpolates[i] = tmp
  end
  return interpolates
end

function interpolate(S::PolyharmonicSpline,x::Float64)
  return interpolate(S,[x])
end

function interpolate(S::PolyharmonicSpline,x::Float64,y::Float64)
  return interpolate(S,[x y])
end

function interpolate(S::PolyharmonicSpline,x::Float64,y::Float64,z::Float64)
  return interpolate(S,[x y z])
end

function interpolate(S::PolyharmonicSpline,x::Vector{Float64})
  return interpolate(S,x'')
end

function interpolate(S::PolyharmonicSpline,x::Vector{Float64},y::Vector{Float64})
  return interpolate(S,[x y])
end

function interpolate(S::PolyharmonicSpline,x::Vector{Float64},y::Vector{Float64},z::Vector{Float64})
  return interpolate(S,[x y z])
end
