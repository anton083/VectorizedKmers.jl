module VectorizedKmers

export
    # AbstractKmerArray
    AbstractKmerArray,

    # Vector
    AbstractKmerVector,
    KmerVector,

    # Matrix
    AbstractKmerMatrix,
    KmerVectors,
    KmerColumns,
    KmerRows,

    # counting
    count_kmers!,
    count_kmers

include("array.jl")
include("vector.jl")
include("vectors.jl")
include("conversion.jl")
include("counting.jl")

end