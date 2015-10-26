# Type Definitions
type BinMapping
    dims::Dims
    bins::Array{Int,1}
    rows::Array{FloatingPoint,1}
    cols::Array{FloatingPoint,1}
    nbins::Int
    nelements::Array{Int,1}
end

type BinnedArray{T,N} <: AbstractArray{T,N}
    map::BinMapping
    values::Array{T,1}
end

function BinnedArray(map,values)
    BinnedArray{eltype(values),length(map.dims)}(map,values)
end

typealias BinnedVector{T} BinnedArray{T,1}
typealias BinnedMatrix{T} BinnedArray{T,2}

# Utility functions
Base.length(B::BinnedArray) = length(B.values)
Base.length(B::BinMapping) = B.nbins

Base.size(B::BinnedArray) = B.map.dims
Base.size(B::BinnedArray,d) = (d > length(B.map.dims) ? 1 : B.map.dims[d])
Base.size(B::BinMapping) = B.dims
Base.size(B::BinMapping,d) = (d > length(B.dims) ? 1 : B.dims[d])

Base.similar{T}(B::BinnedArray{T}) = BinnedArray(B.map, Array(T,length(B)))
Base.similar(B::BinnedArray,t) = BinnedArray(B.map, Array(t,length(B)))

# Get Index functions
function Base.getindex(B::BinnedArray,i::Int)

    binNumber = B.map.bins[i]
    if binNumber != 0
        value = B.values[binNumber]
        return value
    end

    return zero(eltype(B.values))
end

function Base.getindex(B::BinnedArray,inds::Int...)
     i = sub2ind(size(B),inds...)
     return getindex(B,i)
end

function Base.getindex(B::BinMapping,i::Int,j::Int)
    return B.bins[j] == i
end

function Base.getindex(B::BinMapping,i::Int,inds::Int...)

    if sum(map(.<=,inds,size(B))) != length(inds)
        throw(ArgumentError("Dimension mismatch [bin,dims...]"))
    end
    j = sub2ind(B.dims,inds...)
    return getindex(B,i,j)
end

Base.sum(B::BinnedArray) = sum(B.values .* B.binMapping.nelements)

function mapBins{T}(b::BinMapping,v::Vector{T})
    length(v) == b.nbins || throw(ArgumentError("Dimension mismatch"))
    return BinnedArray(b,v)
end

function mapBins(b::BinMapping,A)
    b.dims == size(A) || throw(ArgumentError("Dimension mismatch"))
    v = zeros(eltype(A),b.nbins)

    for j = 1:length(A)
        b.bins[j] == 0 && continue
        v[b.bins[j]] += A[j]
    end

    return BinnedArray(b,v./b.nelements)
end
