"""
    KmerVectors

A container for K-mer vectors, where the values are stored together as columns or rows in a matrix.

`D` is the dimension of the vectors (i.e. `1` for column-vectors, `2` for row-vectors).
`S` is the alphabet size,
`k` is the K-mer length,
and `T` is the element type of the underlying `values` field,
which in turn has type `M`.
"""
struct KmerVectors{D, S, K, T, M} <: AbstractKmerMatrix{S, K, T, M}
    values::M

    function KmerVectors{D, S, K}(values::M) where {D, S, K, T, M <: AbstractMatrix{T}}
        @assert D isa Integer && 1 <= D <= 2 "Type parameter `D` must be an integer between 1 and 2."
        @assert size(values, D) == S^K
        return new{D, S, K, T, M}(values)
    end

    function KmerVectors{D, S, K}(n::Integer; T::Type{<:Real}=Int, zeros::Function=zeros) where {D, S, K}
        dims = (D == 1) ? (S^K, n) : (n, S^K)
        return KmerVectors{D, S, K}(zeros(T, dims))
    end
end

const KmerColumns = KmerVectors{1}
const KmerRows = KmerVectors{2}

# used to distinguish between KmerColumns{4, 2}(4) and KmerRows{4, 1}(16)
Base.hash(kvs::KmerVectors{D, S}, h::UInt) where {D, S} = hash(kvs.values, h ⊻ S ⊻ D)

@inline Base.length(kvs::KmerVectors{D}) where D = size(kvs, 3-D)

@inline function Base.getindex(krs::KmerRows{S, K}, i) where {S, K}
    return KmerVector{S, K}(view(krs.values, i, :))
end

@inline function Base.getindex(kcs::KmerColumns{S, K}, i) where {S, K}
    return KmerVector{S, K}(view(kcs.values, :, i))
end

# could make a proper type for this that isn't just a generator
@inline _eachvec_iter(kvs::KmerVectors{D}) where D = (D == 1) ? eachcol(kvs.values) : eachrow(kvs.values)
@inline eachvec(kvs::KmerVectors{D, S, K}) where {D, S, K} = (KmerVector{S, K}(v) for v in _eachvec_iter(kvs))