type TSVDFunctional{T}
  A::Array{T,2}
  b::Array{T,1}
  Ascale::T
  U::Array{T,2}
  V::Array{T,2}
  S::Array{T,1}
end

function TSVDFunctional{T<:Real}(A::AbstractArray{T,2},
                                 b::AbstractArray{T,1})

  Amax = maximum(A)
  A = A .* (1.0/Amax)
  Ascale = (1.0/Amax)
  U, S, V = svd(A)
  return TSVDFunctional(A,b,Ascale,U,V,S)
end

function TSVDFunctional{T<:Real}(A0::AbstractArray{T,2},
                                 b0::AbstractArray{T,1},
                                 err::AbstractArray{T,1})

  A, b = scale_by_error(A0,b0,err)
  return TSVDFunctional(A,b)
end

function minimize{T<:Real}(RF::TSVDFunctional{T}, n::Int)

  Sinv = zeros(length(RF.S))

  for i=1:min(length(RF.S),n)
      if RF.S[i] > 0.0
        Sinv[i] = 1.0/RF.S[i]
      end
  end

  Adag = RF.V*diagm(Sinv)*RF.U'
  x = Adag*RF.b

  rho = norm(RF.A*x .- RF.b)^2
  xi = norm(x)

  sigma = sqrt(diag(Adag*Adag'))
  
  return RegularizedSolution(RF.Ascale*x,RF.Ascale*sigma,n,rho,xi)
end
