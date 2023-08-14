```@meta
CurrentModule = VectorizedKmers
DocTestSetup = quote
    using VectorizedKmers
end
```

# VectorizedKmers

[![Latest Release](https://img.shields.io/github/release/anton083/VectorizedKmers.jl.svg)](https://github.com/anton083/VectorizedKmers.jl/releases/latest)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://opensource.org/license/MIT)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://anton083.github.io/VectorizedKmers.jl/stable/)
[![Status](https://github.com/anton083/VectorizedKmers.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/anton083/VectorizedKmers.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/anton083/VectorizedKmers.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/anton083/VectorizedKmers.jl)

[VectorizedKmers.jl](https://github.com/anton083/VectorizedKmers.jl) is a Julia package for fast $k$-mer counting of biological sequences. The core idea is that $k$-mers with an alphabet size $A$ are essentially integers in base $A$, and can be used as indices in a vector of size $A^k$ to count the corresponding $k$-mers.

The `KmerCount` type is a wrapper for `AbstractVector`, which means that these vector k-mer counts are not limited to Julia's `Base.Vector` type; other kinds of vectors can be used as well, such as `CUDA.CuVector`, `SparseArrays.SparseVector` or even matrix views. To efficiently group k-mer counts together, the `KmerCountVector` stores them in a matrix as rows or columns. It can wrap any `AbstractMatrix`, such as `Matrix` or `CuMatrix`, and accessing its elements by index returns a `KmerCount` wrapped around a view of a row or column of the original matrix.

This data structure can be used to quickly approximate distances between sequences. Most notably, the squared Euclidean distance was used to estimate edit distance in [this paper](https://doi.org/10.1093/nar/gkz657). The dot product has also proven to be a useful metric for comparing correlation between sequences.

## Installation

You can install VectorizedKmers in your environment from the Julia REPL by entering pkg mode with `]` and then running:

```
add VectorizedKmers
```


```@index
```

```@autodocs
Modules = [VectorizedKmers]
```