using VectorizedKmers
using Documenter

DocMeta.setdocmeta!(VectorizedKmers, :DocTestSetup, :(using VectorizedKmers); recursive=true)

makedocs(;
    modules = [VectorizedKmers],
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
    ),
    sitename = "VectorizedKmers.jl",
    doctest = false,
    pages = [
        "Home" => "index.md",
        "Integer representation of K-mers" => "kmer_int_repr.md",
        "K-mer counting" => Any[
            "Vectors of K-mer counts" => "kmer_count_vector.md",
            "Matrices of K-mers counts" => "kmer_count_matrix.md",
        ],
        "API Reference" => "references.md",
    ],
    authors = "Anton O. Sollman",
    checkdocs = :all,
)

deploydocs(;
    repo = "github.com/anton083/VectorizedKmers.jl.git",
    push_preview = true,
    devbranch = "dev",
    deps = nothing,
    make = nothing,
)
