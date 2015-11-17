function linReg(A0::Array{Float64,2},b0::Array{Float64,1},err::Array{Float64,1},R::Array{Float64,2},alpha::Real; rChisq = 1.0, Ascale = 1.0, Anorm = true )

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

    ATA = A'*A
    ATB = A'*b

    fp = (ATA + alpha.*R)\ATB
    fp[fp .< 0] = 0.0

    rchi2 = (norm(A*fp-b)^2.0)/nr

    if rchi2 <= rChisq + 0.1
        status = :Optimal
    else
        status = :Infeasible
    end

    return fp.*Ascale, rchi2, status
end

function linReg(A0::Array{Float64,2},b0::Array{Float64,1},err::Array{Float64,1},R::Array{Float64,2}; rChisq = 1.0, Ascale = 1.0, Anorm = true)

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

    ATA = A'*A
    ATB = A'*b

    # Try to find best alpha parameter
    alpha_min = 1e-5
    alpha_max = 1e35
    alpha = 0.0
    fp = (ATA + alpha_max.*R)\ATB
    rchi2 = (norm(A*fp-b)^2.0)/nr

    tau = 0.1
    caim = rChisq
    count = 0
    count_max = 2000
    while abs(rchi2 - caim) > tau && count < count_max
        alpha = 0.5*(alpha_min+alpha_max)
        fp = (ATA + alpha.*R)\ATB

        rchi2 = (norm(A*fp-b)^2.0)/nr
        if rchi2 < caim - tau
            alpha_min = alpha
        elseif rchi2 > caim + tau
            alpha_max = alpha
        end
        count = count+1
    end
    fp[fp .< 0] = 0
    rchi2 = (norm(A*fp-b)^2.0)/nr

    count >= count_max && warn(@sprintf("Did not converge: rchi2 = %3.2e",rchi2))

        # if alpha_max - alpha_min < 2e-5
        #     warn("Did not converge")
        #     caim = ceil(rchi2,1)
        #     info("Increasing reduced chi-squared to "*string(round(caim,2)))
        #     alpha_min = 1e-5
        #     alpha_max = 1e15
        #     rchi2 = 1e20
        #     count+=1
        # end
    if rchi2 <= caim + tau
        status = :Optimal
    else
        status = :Infeasible
    end

    return fp.*Ascale, rchi2, status
end

function linRegEP(A0::Array{Float64,2},b0::Array{Float64,1},err::Array{Float64,1},energy::Array{Float64,1},pitch::Array{Float64,1}; method = :Tikhonov0, kwargs...)

    nr,nc = size(A0)
    if method == :Tikhonov0
        L0 = eye(nc)
        R = L0'*L0
        fp, rchi2, status = linReg(A0, b0, err, R; kwargs...)
    elseif method == :Tikhonov1
        L1E, L1p = gradEP(energy,pitch)
        R = L1p'*L1p .+ L1E'*L1E
        fp, rchi2, status = linReg(A0, b0, err, R; kwargs...)
    #elseif method == :TikhonovHybrid
    #     L1E, L1p = gradEP(energy,pitch,kwargs...)
    #     coverage = sum(A0,1)
    #     for i=1:nc
    #         if coverage[i] < eps()
    #             L1E[i,:] = 0.0
    #             L1E[i,i] = 1.0
    #             L1p[i,:] = 0.0
    #             L1p[i,i] = 1.0
    #         end
    #     end
    #     R = L1p'*L1p .+ L1E'*L1E
    #     fp, rchi2, status = linReg(A0, b0, err, R; kwargs...)
    elseif method == :MinFisher
        L1E, L1p = gradEP(energy,pitch)
        # coverage = sum(A0,1)
        #  for i=1:nc
        #       if coverage[i] < eps()
        #           L1E[i,:] = 0.0
        #           L1E[i,i] = 1.0*L1E[i,i]
        #           L1p[i,:] = 0.0
        #           L1p[i,i] = 1.0*L1p[i,i]
        #       end
        #  end
        R = L1p'*L1p .+ L1E'*L1E
        fp, rchi2, status = linReg(A0, b0, err, R; kwargs...)

        W = zeros(nc,nc)
        for i=1:1
            W[:,:] = 0.0
            fp_min = minimum(fp[fp .> eps()])

            for i=1:nc
                fp[i] > 0 ? W[i,i] = 1.0/fp[i] : W[i,i] = 1.0/fp_min
            end

            R = L1p'*W*L1p .+ L1E'*W*L1E
            fp, rchi2, status = linReg(A0, b0, err, R; kwargs...)
        end
    else
        error("Linear regularization method not implemented")
    end

    return fp, rchi2, status
end
