type TSVDFunctional{T}
  A::Array{T,2}
  b::Array{T,1}
  Ascale::T
  U::Array{T,2}
  V::Array{T,2}
  S::Array{T,1}
end

function TSVDFunctional{T<:Real}(A::Array{T,2}, b::Array{T,1})
  Amax = maximum(A)
  A = A .* (1.0/Amax)
  Ascale = (1.0/Amax)
  U, S, V = svd(A)
  return TSVDFunctional(A,b,Ascale,U,V,S)
end

function TSVDFunctional{T<:Real}(A0::Array{T,2}, b0::Array{T,1}, err::Array{T,1})
  A, b = scale_by_error(A0,b0,err)
  return TSVDFunctional(A,b)
end

function minimize{T<:Real}(RF::TSVDFunctional{T}, n::Int; true_sol = Void)

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
  
  if true_sol != Void
    x_true = vec(true_sol)/RF.Ascale
    reg_err = (x_true .- Adag*(RF.A*x_true))
  else
    reg_err = zeros(length(x))
  end

  return RegularizedSolution(RF.Ascale*x,RF.Ascale*sigma,RF.Ascale*reg_err,n,rho,xi)
end

# function tsvd{T}(A0::Array{T,2},b0::Array{T,1},err::Array{Float64,1}; rChisq = 1.0, Ascale = 1.0, Anorm = true,p = 0.8)
#     nr,nc = size(A0)
#
#     # Set a minimum error level
#     err[err .< 0.02*mean(b0)] = 0.02*mean(b0)
#
#     # Scale matrices by error
#     A = zeros(nr,nc)
#     b = zeros(nr)
#
#     for i=1:nr
#         A[i,:] = (Ascale/err[i]).*A0[i,:]
#         b[i] = b0[i]/err[i]
#     end
#
#     if Anorm
#         Amax = maximum(A)
#         A = A .* (1.0/Amax)
#         Ascale = Ascale*(1.0/Amax)
#     end
#
#     # Calculate Singular Value Decomposition
#     U, S, V = svd(A)
#
#     Sinv = zeros(length(S),length(S))
#     f = zeros(nc,length(S))
#     rcs = zeros(length(S))
#     nullrcs = zeros(length(S))
#
#     # Find optimal number of singular values
#     for i=1:length(S)
#         if S[i] <= 0.0
#             # Trim arrays
#             f = f[:,1:i-1]
#             rcs = rcs[1:i-1]
#             nullrcs = nullrcs[1:i-1]
#             break
#         end
#         Sinv[i,i] = 1.0/S[i]
#
#         fp = (V*Sinv*U')*b
#         fp[fp .< 0.0] = 0.0
#         rcs[i] = (norm(A*fp-b)^2.0)/nr
#         nullrcs[i] = (norm(A*fp)^2.0)/nr
#         f[:,i] = fp
#     end
#
#     # Find smallest reduced chi squared that is smaller than p*nullrcs p around 0.8
#     rchi2 = 1e13
#     ind=0
#
#     for i=1:length(rcs)
#         if (rcs[i] <= p*nullrcs[i]) && (rcs[i] <= rchi2)
#           rchi2 = rcs[i]
#           ind = i
#         end
#     end
#
#     if ind == 0
#         rchi2,ind = findmin(rcs)
#     end
#
#     fp = f[:,ind]
#
#     if rchi2 <= rChisq
#         status = :Optimal
#     else
#         status = :Infeasible
#     end
#
#     return fp.*Ascale, rchi2, status
#
# end
