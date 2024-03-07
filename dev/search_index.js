var documenterSearchIndex = {"docs":
[{"location":"references/#API-Reference","page":"API Reference","title":"API Reference","text":"","category":"section"},{"location":"references/","page":"API Reference","title":"API Reference","text":"Modules = [VectorizedKmers]","category":"page"},{"location":"references/#VectorizedKmers.KmerArray","page":"API Reference","title":"VectorizedKmers.KmerArray","text":"KmerArray{N, K, T <: Real, A <: AbstractArray{T, K}} <: StaticArray{NTuple{K, N}, T, K}\n\nN is the alphabet size\nK is the K-mer size\nT is the element type\nA is the array type\n\n\n\n\n\n","category":"type"},{"location":"references/#VectorizedKmers.count_kmers","page":"API Reference","title":"VectorizedKmers.count_kmers","text":"count_kmers(seq, N, K, T=Int, zeros=zeros)\n\n\n\n\n\n","category":"function"},{"location":"references/#VectorizedKmers.count_kmers!","page":"API Reference","title":"VectorizedKmers.count_kmers!","text":"count_kmers!(kmer_array, seq; reset=true)\n\n\n\n\n\n","category":"function"},{"location":"kmer_int_repr/#Integer-representation-of-k-mers","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"","category":"section"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"This package relies on representing k-mers as integers for indexing. This page goes over how that works exactly.","category":"page"},{"location":"kmer_int_repr/#DNA-sequences","page":"Integer representation of k-mers","title":"DNA sequences","text":"","category":"section"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"For DNA, each non-ambiguous nucleotide is assigned a number between 0 and 3:","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"Nucleotide Base-4 Base-2\nA 0 00\nC 1 01\nG 2 10\nT 3 11","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"note: Note\nAny ordering works, but this is the one used by BioSequences.jl. It also has some nice properties, like being in alphabetical order, and that XOR-ing a base with 3 gives you its complement.","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"We could theoretically convert any DNA sequence to an integer, but 64-bit unsigned integers limit us to 32-mers.","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"Consider the DNA sequence GATTACA. If we convert it to an integer using the table above, we get 2033010_4 = 10001111000100_2 = 9156_10, so the integer value of GATTACA is 9156. Since Julia uses 1-based indexing, we would add 1 to this value to get the index for the value in a vector associated with GATTACA.","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"Writing a function for this might look like the following:","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"julia> const DNA_ENCODING_VECTOR = zeros(UInt, 127);\n\njulia> for (i, char) in enumerate(\"ACGT\")\n           DNA_ENCODING_VECTOR[char % Int8] = i - 1\n       end;\n\njulia> function kmer_to_int(kmer::String)\n           kmer_int = zero(UInt)\n           for char in kmer\n               kmer_int = (kmer_int << 2) | DNA_ENCODING_VECTOR[char % Int8]\n           end\n           return kmer_int\n       end;\n\njulia> f(k;i=0)=[i=4i|(c%Int-1-(c=='C'))&3 for c=k][end]; # or if you're into code golf:","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"note: Note\nConsider using types like BioSequences.LongDNA{4} instead of String, as strings are not optimized for performance.","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"This function would not be very efficient in a practical setting (even though we're using super cool bit manipulation), since we're convert each k-mer individually, instead of having some kind of sliding window. Moreover, the function takes the k-mer in the form of a String, which is not ideal. The function should work as intended, though. Let's test it:","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"julia> kmer_to_int(\"GATTACA\")\n9156\n\njulia> f(\"GATTACA\")\n9156","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"This package does not use vectors for this, but rather K-dimensional arrays, with each dimension having a size of 4, or whatever else the alphabet size might be. Linear indexing of these arrays is still possible, and they can be lazily reshaped to vectors.","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = VectorizedKmers\nDocTestSetup = quote\n    using VectorizedKmers\nend","category":"page"},{"location":"#VectorizedKmers","page":"Home","title":"VectorizedKmers","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"(Image: Latest Release) (Image: MIT license) (Image: Documentation) (Image: Documentation) (Image: Status) (Image: Coverage)","category":"page"},{"location":"","page":"Home","title":"Home","text":"VectorizedKmers.jl is a Julia package primarily designed for fast K-mer counting of biological sequences. The core idea is that K-mers with an alphabet size of N are essentially integers in base N, and can be used as indices in a vector of size N^K to count the corresponding K-mers.","category":"page"},{"location":"","page":"Home","title":"Home","text":"This data structure can be used to quickly approximate distances between sequences. Notably, the squared Euclidean distance was used to approximate edit distance in this paper. The dot product has also proven to be a useful metric for comparing correlation between sequences.","category":"page"},{"location":"#Examples","page":"Home","title":"Examples","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"julia> using VectorizedKmers, BioSequences\n\njulia> kmer_array = count_kmers(dna\"ACCGGGTTTT\", 1)\nKmerArray{4, 1, Int64, Vector{Int64}} with size (4,)\n\njulia> kmer_array |> values\n4-element Vector{Int64}:\n 1\n 2\n 3\n 4\n\njulia> count_kmers(dna\"AATT\", 2) |> values # 2-mers of AATT\n4×4 Matrix{Int64}:\n 1  0  0  0\n 0  0  0  0\n 0  0  0  0\n 1  0  0  1\n\njulia> count_kmers(aa\"AY\", 1) |> values\n20-element Vector{Int64}:\n 1\n 0\n 0\n ⋮\n 0\n 1\n 0","category":"page"},{"location":"","page":"Home","title":"Home","text":"For more examples, including some technical stuff, see the documentation.","category":"page"},{"location":"#Limitations","page":"Home","title":"Limitations","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The main downside of counting K-mers this way is that the arrays grow exponentially with respect to K. The 31-mer array of a DNA sequence would have a length of 4^31 = 4611686018427387904, which is equivalent to four exbibytes of memory, if the values are stored with 8-bit integers — which is just not feasible, really. Not only does allocating a lot of memory take up a lot of memory, but it can also take a substantial amount of time! This method of counting K-mers therefore works best for lower K-values.","category":"page"}]
}
