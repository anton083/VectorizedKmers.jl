module BioSequencesExt
 
using VectorizedKmers, BioSequences

const LongNucOrView{N} = Union{LongSequence{<:NucleicAcidAlphabet{N}}, LongSubSeq{<:NucleicAcidAlphabet{N}}}

VectorizedKmers.alphabet_size(::Type{<:LongNucOrView}) = 4

@inline function VectorizedKmers.count_kmers!(
    kmer_vector::KmerVector{4, k, T, A},
    sequence::LongNucOrView{2};
    reset::Bool = true,
) where {k, T, A <: AbstractVector{T}}
    reset && VectorizedKmers.zeros!(kmer_vector)
    values = kmer_vector.values
    len = length(sequence)
    len < k && return kmer_vector
    mask = UInt(4^k - 1)
    kmer = UInt(0)
    start, stop = sequence isa LongSubSeq ? (sequence.part.start, sequence.part.stop) : (1, len)
    data_start, data_stop = (start - 1) ÷ 32 + 1, (stop - 1) ÷ 32 + 1
    first_count_index = k + start - 1
    i = 32 * (data_start - 1)
    @inbounds for data_int in @view sequence.data[data_start:data_stop]
        for j in 0:2:63
            i += 1
            i > stop && break
            kmer = ((kmer << 2) & mask) | ((data_int >> j) & 0b11)
            values[kmer + 1] += first_count_index <= i
        end
    end
    return kmer_vector
end

@inline function VectorizedKmers.count_kmers!(
    kmer_vector::KmerVector{4, k, T, A},
    sequence::LongNucOrView{4};
    reset::Bool = true,
) where {k, T, A <: AbstractVector{T}}
    reset && VectorizedKmers.zeros!(kmer_vector)
    values = kmer_vector.values
    len = length(sequence)
    len < k && return kmer_vector
    mask = UInt(4^k - 1)
    kmer = UInt(0)
    start, stop = sequence isa LongSubSeq ? (sequence.part.start, sequence.part.stop) : (1, len)
    data_start, data_stop = (start - 1) ÷ 16 + 1, (stop - 1) ÷ 16 + 1
    first_count_index = k + start - 1
    i = 16 * (data_start - 1)
    @inbounds for data_int in @view sequence.data[data_start:data_stop]
        for j in 0:4:63
            i += 1
            i > stop && break
            kmer = ((kmer << 2) & mask) | (trailing_zeros(data_int >> j) & 0b11)
            values[kmer + 1] += first_count_index <= i
        end
    end
    return kmer_vector
end

@inline function VectorizedKmers.count_kmers!(
    kmer_vector::KmerVector{4, k, Bool, BitVector},
    sequence::LongNucOrView{4};
    reset::Bool = true,
) where k
    chunks = kmer_vector.values.chunks
    reset && fill!(chunks, zero(UInt))
    len = length(sequence)
    len < k && return kmer_vector
    mask = UInt(4^k - 1)
    kmer = UInt(0)
    start, stop = sequence isa LongSubSeq ? (sequence.part.start, sequence.part.stop) : (1, len)
    data_start, data_stop = (start - 1) ÷ 16 + 1, (stop - 1) ÷ 16 + 1
    first_count_index = k + start - 1
    i = 16 * (data_start - 1)
    @inbounds for data_int in @view sequence.data[data_start:data_stop]
        for j in 0:4:63
            i += 1
            i > stop && break
            kmer = ((kmer << 2) & mask) | (trailing_zeros(data_int >> j) & 0b11)
            chunk_index = kmer >> 6 + one(UInt)
            index_in_chunk = kmer & 0x3f
            first_count_index <= i && (chunks[chunk_index] ⊻= one(UInt) << index_in_chunk)
        end
    end
    return kmer_vector
end

end