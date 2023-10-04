module BioSeqCUDAExt

using VectorizedKmers, BioSequences, CUDA

@inline data_length(seq::LongDNA{4}) = length(seq.data) 

function VectorizedKmers.count_kmers!(
    kmer_columns::KmerColumns{4, k, T, M},
    sequences::Vector{LongDNA{4}};
    offset::Integer = 0,
    reset::Bool = true,
) where {k, T, M <: CuMatrix{T}}
    values = kmer_columns.values
    reset && CUDA.fill!(values, 0)
    isempty(sequences) && return kmer_columns

    num_seqs = length(sequences)

    seq_lengths_h = length.(sequences)
    max_seq_len = maximum(seq_lengths_h)
    seq_lengths = CuVector(seq_lengths_h)

    # data vector length could be calculated with (seq_len - 1) ÷ 16 + 1
    data_lengths_h = data_length.(sequences)
    max_data_length = maximum(data_lengths_h)
    data_lengths = CuVector(data_lengths_h)

    data_matrix_h = Matrix{UInt64}(undef, (max_data_length, num_seqs))
    for (i, seq) in enumerate(sequences)
        data_matrix_h[1:data_length(seq), i] = seq.data
    end

    data_matrix = CuMatrix(data_matrix_h)

    @assert length(kmer_columns) >= num_seqs "The k-mer counts of $num_seqs sequences would not fit in $(length(kmer_columns)) columns."
    @assert length(kmer_columns) - offset >= num_seqs "The column offset is too high. The sequences don't fit."
    max_seq_len > k-1 || @warn "All sequences are shorter than k. No k-mers will be counted."

    function count_kmers_column!(
        values,
        data_matrix,
        seq_lengths,
        num_seqs,
        k,
        mask,
        data_lengths,
        offset,
    )
        seq_idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x

        @inbounds if seq_idx <= num_seqs
            seq_length = seq_lengths[seq_idx]
            data_length = data_lengths[seq_idx]

            count_vector = view(values, :, seq_idx + offset)
            data_vector = view(data_matrix, 1:data_length-1, seq_idx)

            kmer = UInt(0)
            i = 0
            for data_int in data_vector
                for j in 0:4:63
                    i += 1
                    kmer = ((kmer << 2) & mask) | (trailing_zeros(data_int >> j) & 0b11)
                    count_vector[kmer + 1] += k <= i
                end
            end

            # TODO: make sure this works correctly. it's not tested.
            data_int = data_matrix[data_length, seq_idx]
            for j in 0:4:(4 * ((seq_length - 1) % 16 + 1) - 1)
                i += 1
                kmer = ((kmer << 2) & mask) | (trailing_zeros(data_int >> j) & 0b11)
                count_vector[kmer + 1] += k <= i
            end
        end

        return nothing
    end

    mask = unsigned(4^k - 1)
    threads = 256
    blocks = ceil(Int, num_seqs / threads)

    @cuda threads=threads blocks=blocks count_kmers_column!(
        values, data_matrix, seq_lengths, num_seqs, k, mask, data_lengths, offset)

    return kmer_columns
end

end