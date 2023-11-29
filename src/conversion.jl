@inline Base.transpose(kvs::KmerVectors{D, S, K}) where {D, S, K} = KmerVectors{3-D, S, K}(transpose(kvs.values))
@inline Base.adjoint(kvs::KmerVectors) = transpose(kvs)

function KmerVector(f::Function, kvs::KmerVectors{D, S, K}) where {D, S, K}
    return KmerVector{S, K}(reshape(mapslices(f, kvs.values, dims=3-D), :))
end