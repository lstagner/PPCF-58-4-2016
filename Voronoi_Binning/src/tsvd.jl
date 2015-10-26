function tsvd(A0::Array{Float64,2},b0::Array{Float64,1},err::Array{Float64,1},n::Int64; Ascale = 1.0, Anorm = true)
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

    if Anorm
        Amax = maximum(A)
        A = A .* (1.0/Amax)
        Ascale = Ascale*(1.0/Amax)
    end

    U, S, V = svd(A)

    Sinv = zeros(length(S))

    for i=1:min(length(S),n)
        if S[i] > 0.0
          Sinv[i] = 1.0/S[i]
        end
    end

    fp = (V*diagm(Sinv)*U')*b

    fp[fp .< 0] = 0.0

    rchisq = (norm(A*fp-b)^2.0)/nr

    status = :Optimal

    return fp.*Ascale, rchisq, status

end

function tsvd(A0::Array{Float64,2},b0::Array{Float64,1},err::Array{Float64,1}; rChisq = 1.0, Ascale = 1.0, Anorm = true,p = 0.8)
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

    if Anorm
        Amax = maximum(A)
        A = A .* (1.0/Amax)
        Ascale = Ascale*(1.0/Amax)
    end

    # Calculate Singular Value Decomposition
    U, S, V = svd(A)

    Sinv = zeros(length(S),length(S))
    f = zeros(nc,length(S))
    rcs = zeros(length(S))
    nullrcs = zeros(length(S))

    # Find optimal number of singular values
    for i=1:length(S)
        if S[i] <= 0.0
            # Trim arrays
            f = f[:,1:i-1]
            rcs = rcs[1:i-1]
            nullrcs = nullrcs[1:i-1]
            break
        end
        Sinv[i,i] = 1.0/S[i]

        fp = (V*Sinv*U')*b
        fp[fp .< 0.0] = 0.0
        rcs[i] = (norm(A*fp-b)^2.0)/nr
        nullrcs[i] = (norm(A*fp)^2.0)/nr
        f[:,i] = fp
    end

    # Find smallest reduced chi squared that is smaller than p*nullrcs p around 0.8
    rchi2 = 1e13
    ind=0

    for i=1:length(rcs)
        if (rcs[i] <= p*nullrcs[i]) && (rcs[i] <= rchi2)
          rchi2 = rcs[i]
          ind = i
        end
    end

    if ind == 0
        rchi2,ind = findmin(rcs)
    end

    fp = f[:,ind]

    if rchi2 <= rChisq
        status = :Optimal
    else
        status = :Infeasible
    end

    return fp.*Ascale, rchi2, status

end
