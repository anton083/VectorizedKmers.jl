module SparseArraysExt

using VectorizedKmers, SparseArrays
import VectorizedKmers: AbstractKmerCount

# this probably doesn't even deserve to be an extension. was really easy to implement and build on default types
# SparseKmerCount type not necessary since KmerCount now wraps AbstractVector
# this should be scrapped unless it can be generalized

struct SparseKmerCount{A, K, T, V} <: AbstractKmerCount{A, K, T, V}
    count::SparseVector{T}

    function SparseKmerCount{A, K, T}() where {A, K, T}
        new{A, K, T, SparseVector}(spzeros(T, A^K))
    end
end

#= Potential use-case:

"""
    sparse_aa_kmer_count(seq, K)

Count amino acid kmers of size `K` in seq, returns a sparse array, where each kmer has been mapped
to a cell and the value in that cell is the number of occurences of that kmer within `seq`.
"""
function sparse_aa_kmer_count(seq::LongDNA, K::Integer)
    sp_kmer_count = KmerCount{28, K, Int}(spzeros(T, A^K))
    mask = UInt(28)^K
    kmer = UInt(0)
    for aa_seq in translate_all_frames(seq)
        for aa in view(aa_seq, 1:K-1)
            kmer = kmer * 28 + reinterpret(UInt8, aa)
        end
        for aa in view(aa_seq, K:length(aa_seq))
            kmer = kmer * 28 % mask + reinterpret(UInt8, aa)
            sp_kmer_count[kmer + 1] += 1
        end
    end
    sp_kmer_count
end
=#

end