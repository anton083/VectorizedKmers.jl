"""
    KmerVector{S, K, T, V} <: AbstractKmerVector{S, K, T, V}

`A` is the alphabet size,
`K` is the K-mer size,
and `T` is the element type of the underlying `values` field,
which in turn has type `V`.
"""
struct KmerVector{S, K, T, V} <: AbstractKmerVector{S, K, T, V}
    values::V

    function KmerVector{S, K}(values::V) where {S, K, T, V <: AbstractVector{T}}
        @assert length(values) == S^K
        return new{S, K, T, V}(values)
    end

    function KmerVector{S, K}(; T::Type{<:Real}=Int, zeros::Function=zeros) where {S, K}
        return KmerVector{S, K}(zeros(T, S^K))
    end
end

Base.summary(kv::KmerVector) = string(typeof(kv))

function Base.show(io::IO, kv::KmerVector)
    max_elements_to_show = 10
    len_values = length(kv.values)
    
    print(io, summary(kv) * "([")
    elements_to_show = min(len_values, max_elements_to_show)
    print(io, join(kv.values[1:elements_to_show], ", "))
    elements_to_show < len_values && print(io, " … (", len_values - max_elements_to_show, " more)")
    print(io, "])")

    return nothing
end
