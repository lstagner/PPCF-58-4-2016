function binAccretion(z::Array{Float64,2},target::Real)

    nr, nc = size(z)
    nind = length(z)
    zz = reshape(z,nind)
    rows = reshape([1:nr]*ones(Float64,nc)',nind)
    cols = reshape(ones(Float64,nr)*[1:nc]',nind)

    bins = zeros(Int,nind)
    candidateBins = zeros(Int,nind)
    binValues = zeros(nind)
    rCenters = zeros(nind)
    cCenters = zeros(nind)
    unassigned = ones(Bool,nind)
    goodIndex  = ones(Bool,nind)
    goodBin = zeros(Bool,nind)

    # Start at largest value
    startValue, startIndex = findmax(zz)

    for i=1:nind
        currentBinValue = startValue
        unassigned[startIndex] = false
        rCen, cCen = rows[startIndex], cols[startIndex]
        bins[startIndex] = i

        # Add pixels to bin while critera is met
        while true
            sum(unassigned) == 0 && break
            # Find pixel closest unassigned and non-zero bin to bin generator
            tmp1 = 1e6
            candidateIndex=1
            @inbounds for j = 1:nind
                (unassigned[j] && goodIndex[j]) || continue
                tmp2 = ((rows[j]-rCen)^2 + (cols[j]-cCen)^2)
                if tmp2 < tmp1
                    candidateIndex = j
                    tmp1 = tmp2
                end
            end
            ri = rows[candidateIndex]
            ci = cols[candidateIndex]
            candidateValue = zz[candidateIndex]

            # Find the minimum distance of candidate pixel to the pixels in the bin
            minDist = 1e6
            @inbounds for j=1:nind
                bins[j] == i || continue
                tmp1 = ((rows[j]-ri)^2 + (cols[j]-ci)^2)
                if tmp1 < minDist
                    minDist = tmp1
                end
            end

            # Calculate bin values if index is added
            candidateBins[:] = bins

            candidateBins[candidateIndex] = i

            # Calculate "roundness" of bin
            candidateBinArea = 0.0
            candidateRCen = 0.0
            candidateCCen = 0.0
            @inbounds for j=1:nind
                candidateBins[j] == i || continue
                candidateRCen += rows[j]
                candidateCCen += cols[j]
                candidateBinArea += 1.0
            end
            candidateRCen = candidateRCen/candidateBinArea
            candidateCCen = candidateCCen/candidateBinArea
            candidateBinRadius = sqrt(candidateBinArea/pi)

            maxDistSq = 0.0
            @inbounds for j=1:nind
                candidateBins[j] == i || continue
                tmp1 = ((rows[j] - candidateRCen)^2 + (cols[j] - candidateCCen)^2)
                if tmp1 > maxDistSq
                    maxDistSq = tmp1
                end
            end

            roundness = sqrt(maxDistSq)/candidateBinRadius - 1.0

            # Calculate bin value if index is added
            candidateBinValue = currentBinValue + candidateValue

            # If adding the index to the bin violates the criteria
            if sqrt(minDist) > 1.2 || roundness > 0.3 ||
               abs(candidateBinValue - target) > abs(currentBinValue-target)
                if currentBinValue > 0.8*target
                    goodBin[i] = true
                end
               break
            end

            # Add index to bin
            currentBinValue = candidateBinValue
            bins[candidateIndex] = i

            # Update bin generator
            rCen, cCen = candidateRCen, candidateCCen

            unassigned[candidateIndex] = false
        end

        # If the bin is a good bin then keep it else
        if !goodBin[i]
            @inbounds for j=1:nind
                bins[j] == i || continue
                goodIndex[j] = false
            end
        end
        binValues[i] = currentBinValue
        rCenters[i] = rCen
        cCenters[i] = cCen

        sum(unassigned) == 0 && break

        # Find closest unassigned index to the centroid
        # of all the assigned indices and start a new bin
        rcen = 0.0
        ccen = 0.0
        areaTot = 0.0
        @inbounds for j=1:nind
            (!unassigned[j] && goodIndex[j]) || continue
            rcen += rows[j]
            ccen += cols[j]
            areaTot += 1.0
        end
        rcen = rcen/areaTot
        ccen = ccen/areaTot

        tmp1 = 1e6
        @inbounds for j=1:nind
            (unassigned[j] && goodIndex[j]) || continue
            tmp2 = (rows[j]-rcen)^2 + (cols[j]-ccen)^2
            if tmp2 < tmp1
                startIndex = j
                tmp1 = tmp2
            end
        end
        rCen = rows[startIndex]
        cCen = cols[startIndex]
        startValue = zz[startIndex]
    end

    # If there are any remaining indices are not good
    if sum(unassigned) != 0
        goodIndex[unassigned] = false
    end

    #bins = bins[:,goodBin]
    binNums = [1:nind][goodBin]
    rCenters = rCenters[goodBin]
    cCenters = cCenters[goodBin]
    binValues = binValues[goodBin]
    nbins = length(rCenters)
    newBins = zeros(Int,nind)

    # Renumber bins
    for i=1:nbins
        newBins[bins .== binNums[i]] = i
    end

    bins[:] = newBins

    # Assign bad indices to closest bin
    @inbounds for i=1:nind
        goodIndex[i] && continue
        binIndex = 1
        tmp1 = 1e6
        for j=1:nbins
            tmp2 = (rCenters[j] - rows[i])^2 + (cCenters[j] - cols[i])^2
            if tmp2 < tmp1
                tmp1 = tmp2
                binIndex = j
            end
        end
        zz[i] == 0 && continue
        bins[i] = binIndex
        binValues[binIndex] += zz[i]
        unassigned[i] = false
    end

    # Recalculate bin centroids
    for i=1:nbins
        rCenters[i] = 0.0
        cCenters[i] = 0.0
        tot = 0.0
        @inbounds for j = 1:nind
            bins[j] == i || continue
            rCenters[i] += rows[j]
            cCenters[i] += cols[j]
            tot += 1.0
        end
        rCenters[i] = rCenters[i]/tot
        cCenters[i] = cCenters[i]/tot
    end

    return bins,binValues,rCenters,cCenters,rows,cols
end

function reassign!(bins::Array{Int},rCenters::Array{Float64},
                   cCenters::Array{Float64},scales::Array{Float64},
                   rows::Array{Float64},cols::Array{Float64})

    nind = length(bins)
    nbins = length(rCenters)

    @inbounds for i=1:nind
        minValue = 1e6
        k = 1
        ri = rows[i]
        ci = cols[i]
        for j=1:nbins
            v = ((rCenters[j] - ri)^2 + (cCenters[j] - ci)^2)/(scales[j])^2
            if v < minValue
                minValue = v
                k = j
            end
        end
        bins[i] = k
    end
end

function optimize!(bins::Array{Int},rCenters::Array{Float64},
                   cCenters::Array{Float64},rows::Array{Float64},
                   cols::Array{Float64},z::Array{Float64}; wvt=false)

    nind = length(z)
    nr, nc = size(z)
    zz = reshape(z,nind)
    nbins = length(rCenters)

    rCentersOld = copy(rCenters)
    cCentersOld = copy(cCenters)

    scales = ones(nbins)
    if wvt
        dens = ones(nind)
    else
        dens = zz
    end

    difference = 0.0

    # Try to optimize for at most ninds times
    for i=1:nind
        reassign!(bins,rCenters,cCenters,scales,rows,cols)
        difference = 0.0
        for j = 1:nbins
            # Compute new bin centroids
            rCenters[j] = 0.0
            cCenters[j] = 0.0
            binMass = 0.0
            binArea = 0.0
            binValue = 0.0
            for k = 1:nind
               bins[k] == j || continue
               rCenters[j] += rows[k]*dens[k]^2
               cCenters[j] += cols[k]*dens[k]^2
               binMass += dens[k]^2
               binValue += zz[k]
               binArea += 1.0
            end
            rCenters[j] = rCenters[j]/binMass
            cCenters[j] = cCenters[j]/binMass
            if wvt
                scales[j] = sqrt(binArea/binValue)
            end

            difference += (rCenters[j] - rCentersOld[j])^2 +
                          (cCenters[j] - cCentersOld[j])^2
            rCentersOld[j] = rCenters[j]
            cCentersOld[j] = cCenters[j]
        end

        # Check for convergence
        if difference < eps()
            break
        end
    end

    if difference > eps()
        reassign!(bins,rCenters,cCenters,scales,nr,nc)
        warn("Optimization stage did not converge: diff = $diff")
    end
end

function voronoiBinning(z::Array{Float64,2},target::Real,;wvt = false,threshold=1.0)

    if sum(z) < target || minimum(z) > target
        throw(ArgumentError("Target value infeasible"))
    end

    bins,values,rowCenters,colCenters,rows,cols = binAccretion(z,target)
    nbins = length(values)

    optimize!(bins,rowCenters,colCenters,rows,cols,z,wvt=wvt)

    zz = reshape(z,length(z))
    bins[zz .<= threshold] = 0
    values = zeros(Float64,nbins)
    areas = zeros(Int64,nbins)
    for i=1:nbins
        values[i] = sum(zz[bins .== i])
        areas[i] = sum(bins .== i)
    end

    return BinnedArray(BinMapping(size(z),bins,rowCenters,colCenters,nbins,areas),values./areas)
end

function plotBins(A::BinnedArray;n=11,bin=nothing)
    plotBins(A.map,n=n,bin=bin)
end

function plotBins(B::BinMapping; n=11,bin=nothing)

    nx,ny = B.dims
    A = zeros(nx*ny)
    bins = B.bins
    nbins = B.nbins
    for i=1:nbins
        if bin != nothing
            i != bin && continue
        end
        A[bins .== i] = i % n
    end

    pcolor(reshape(A,nx,ny))
    colorbar()
end

function plotBins(bins::Array{Int,1},nx,ny; n=11,bin=nothing)

    A = zeros(nx*ny)

    for i=1:length(A)
        if bin != nothing
            bins[i] != bin && continue
        end
        A[i] = bins[i] % n
    end

    pcolor(reshape(A,nx,ny))
    colorbar()
end

function plotBinValues(B::BinnedArray)

    nx,ny = size(B)

    A = zeros(nx*ny)
    for i = 1:length(A)
        A[i] = B[i]
    end

    pcolor(reshape(A,nx,ny))
    colorbar()
end

function plotBinIntegral(B::BinnedArray)

    nx,ny = size(B)

    A = zeros(nx*ny)
    for i = 1:length(A)
        val = B[i]
        val == 0 && continue
        area = B.map.nelements[val .== B.values][1]
        A[i] = val*area
    end

    pcolor(reshape(A,nx,ny))
    colorbar()
end
