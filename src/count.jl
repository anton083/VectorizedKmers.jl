"""
    count_kmers!(kmer_array, sequence; reset=true)

Requires method `axis_index(::KmerArray{N}, ::eltype(sequence)) where N` to be defined
"""
function count_kmers!(kmer_array::KmerArray{N, K}, sequence; reset::Bool = true) where {N, K}
    reset && zeros!(kmer_array)
    kmer_array_values = kmer_array.values
    mask = UInt(N)^K
    kmer = zero(UInt)
    for (i, m) in enumerate(sequence)
        kmer = (kmer * N + axis_index(kmer_array, m)) % mask
        kmer_array_values[kmer + 1] += K <= i
    end
    return kmer_array
end

"""
    count_kmers(sequence, [N,] K, T=Int, zeros=zeros)
"""
count_kmers(sequence, ::Val{N}, ::Val{K}, T::Type{<:Real}=Int, zeros=zeros) where {N, K} = count_kmers!(KmerArray(N, K, T, zeros), sequence)
count_kmers(sequence, N::Integer, K::Integer, args...) = count_kmers(sequence, Val(N), Val(K), args...)
count_kmers(sequence, K::Integer, args...) = count_kmers(sequence, default_alphabet_size(eltype(sequence)), K, args...)