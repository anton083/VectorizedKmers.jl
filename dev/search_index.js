var documenterSearchIndex = {"docs":
[{"location":"kmer_int_repr/#Integer-representation-of-k-mers","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"","category":"section"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"This package relies on representing K-mers as integers for indexing, and understanding how it works is recommended (unless you're only using higher-level API stuff).","category":"page"},{"location":"kmer_int_repr/#DNA-sequences","page":"Integer representation of k-mers","title":"DNA sequences","text":"","category":"section"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"For DNA, each non-ambiguous nucleotide is assigned a number between 0 and 3:","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"Nucleotide Base-4 Base-2\nA 0 00\nC 1 01\nG 2 10\nT 3 11","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"note: Note\nAny ordering works, but this one is the one used by BioSequences.jl, and it also has some nice properties, like being in alphabetical order, and that XOR-ing a base with 3 gives you its complement.","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"We could theoretically convert any DNA sequence to an integer, but 64-bit integers limit us to 32-mers.","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"Consider the DNA sequence GATTACA. If we convert it to an integer using the table above, we get 2033010_4 = 10001111000100_2 = 9156_10, so the integer value of GATTACA is 9156. Since Julia uses 1-based indexing, we would add 1 to this value to get the index of GATTACA in the vector.","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"If we want to write a function for this, we may do something like the following:","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"const DNA_ENCODING_VECTOR = zeros(Int, 127)\n\nfor (i, char) in enumerate(\"ACGT\")\n    # `c % Int8` is the ASCII value of the character `c`.\n    DNA_ENCODING_VECTOR[char % Int8] = i - 1\nend\n\nfunction kmer_to_int(kmer::String)\n    kmer_int = 0\n    for char in kmer\n        kmer_int = (kmer_int << 2) + DNA_ENCODING_VECTOR[charh % Int8]\n    end\n    kmer_int\nend\n\n# or if you like playing code golf:\nf(k;i=0)=[i=(c%Int-1-(c=='C'))&3|i<<2 for c=k][end]","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"note: Note\nStrings are bad! Don't use strings! Please! They're bad! Very bad! Chars have variable length when part of Strings in Julia, so indexing and taking lengths and stuff is kinda slow! Don't use strings! Again: very bad! Use something like LongDNA{4}, please!","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"This function would not be very efficient in a practical setting (even though we're using super cool bit manipulation stuff), and converting each k-mer individually, instead of having some kind of sliding window. Moreover, the function takes the k-mer in the form of a String, which is not ideal. The function should work as intended, though. Let's test it:","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"julia> kmer_to_int(\"GATTACA\")\n9156","category":"page"},{"location":"kmer_int_repr/#Amino-acid-sequences","page":"Integer representation of k-mers","title":"Amino acid sequences","text":"","category":"section"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"Amino acid sequences are a little more difficult to deal with since there are a lot more of them, and the vectors would grow in size even quicker. However, we can still represent them as integers, but we can't use bit manipulation anymore.","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"BioSequences.jl has 28 amino acids in its AminoAcidAlphabet, so we can represent each amino acid as an integer between 0 and 27.","category":"page"},{"location":"kmer_count/","page":"The KmerCount type","title":"The KmerCount type","text":"CurrentModule = VectorizedKmers\nDocTestSetup = quote\n    using VectorizedKmers\nend","category":"page"},{"location":"kmer_count/#The-KmerCount-type","page":"The KmerCount type","title":"The KmerCount type","text":"","category":"section"},{"location":"kmer_count/","page":"The KmerCount type","title":"The KmerCount type","text":"The KmerCount type has four type parameters, but you only really need to care about the first two: A, the alphabet size, and K, the K-mer length. So, to count the 6-mers of a DNA sequence, you would use KmerCount{4, 6}. For each of these K-mer counts, memory for a vector of size A^K is allocated, unless a vector type like SparseVector is used. This brings us to the two other type parameters: T, which is the vector element type, and V, which is the type of the actual vector.","category":"page"},{"location":"kmer_count/","page":"The KmerCount type","title":"The KmerCount type","text":"Let's see it in action:","category":"page"},{"location":"kmer_count/","page":"The KmerCount type","title":"The KmerCount type","text":"julia> kc = KmerCount{4, 2}(); # creates a Vector{Int} of zeros with length 4^2\n\njulia> using BioSequences # a weak dependency that lets us count kmers of LongDNA{4} sequences\n\njulia> count_kmers!(kc, dna\"ACGT\"); # AC is 0001, CG is 0110, GT is 1011\n\njulia> @show kc; # indices offset by 1. thanks julia...\nkc = [0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0]\n\njulia> count_kmers!(kc, dna\"ACGT\"); # count again\n\njulia> @show kc; # the vector gets reset to zeros before counting again\nkc = [0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0]\n\njulia> count_kmers!(kc, dna\"ACGT\", reset=false); # avoid reset with reset=false\n\njulia> @show kc; # look! we counted 2-mers of ACGT twice\nkc = [0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0]","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = VectorizedKmers\nDocTestSetup = quote\n    using VectorizedKmers\nend","category":"page"},{"location":"#VectorizedKmers","page":"Home","title":"VectorizedKmers","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"(Image: Latest Release) (Image: MIT license) (Image: Documentation) (Image: Status) (Image: Coverage)","category":"page"},{"location":"","page":"Home","title":"Home","text":"VectorizedKmers.jl is a Julia package for fast k-mer counting of biological sequences. The core idea is that k-mers with an alphabet size A are essentially integers in base A, and can be used as indices in a vector of size A^k to count the corresponding k-mers.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The KmerCount type is a wrapper for AbstractVector, which means that these vector k-mer counts are not limited to Julia's Base.Vector type; other kinds of vectors can be used as well, such as CUDA.CuVector, SparseArrays.SparseVector or even matrix views. To efficiently group k-mer counts together, the KmerCountVector stores them in a matrix as rows or columns. It can wrap any AbstractMatrix, such as Matrix or CuMatrix, and accessing its elements by index returns a KmerCount wrapped around a view of a row or column of the original matrix.","category":"page"},{"location":"","page":"Home","title":"Home","text":"This data structure can be used to quickly approximate distances between sequences. Most notably, the squared Euclidean distance was used to estimate edit distance in this paper. The dot product has also proven to be a useful metric for comparing correlation between sequences.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"You can install VectorizedKmers in your environment from the Julia REPL by entering pkg mode with ] and then running:","category":"page"},{"location":"","page":"Home","title":"Home","text":"add VectorizedKmers","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [VectorizedKmers]","category":"page"},{"location":"#VectorizedKmers.AbstractKmerCount","page":"Home","title":"VectorizedKmers.AbstractKmerCount","text":"AbstractKmerCount{A, K, T <: Real, V <: AbstractVector{T}}\n\nAbstract type for K-mer counts. A is the alphabet size, K is the K-mer size, and T is the element type of the underlying counts field, which in turn has type V{T}.\n\n\n\n\n\n","category":"type"},{"location":"#VectorizedKmers.AbstractKmerCountVector","page":"Home","title":"VectorizedKmers.AbstractKmerCountVector","text":"AbstractKmerCountVector{A, K, T <: Real, M <: AbstractMatrix{T}} <: AbstractVector{KmerCount{A, K, T, V} where {V <: AbstractVector{T}}}\n\nYup... that's indeed an abomination of a type. A container for K-mer counts, where K-mer counts are stored together as rows or columns in a matrix. A is the alphabet size, K is the K-mer size, T is the element type of the counts, and M is the type of the matrix in which the K-mer counts are stored.\n\n\n\n\n\n","category":"type"},{"location":"#VectorizedKmers.KmerCount","page":"Home","title":"VectorizedKmers.KmerCount","text":"KmerCount{A, K, T, V} <: AbstractKmerCount{A, K, T, V}\n\nA concrete type for K-mer counts with vector type V and element type T.\n\n\n\n\n\n","category":"type"},{"location":"#VectorizedKmers.KmerCountColumns","page":"Home","title":"VectorizedKmers.KmerCountColumns","text":"KmerCountColumns{A, K, T, M} <: AbstractKmerCountVector{A, K, T, M}\n\nA container for K-mer counts, where K-mer counts are stored together as columns in a matrix. This is more efficient than storing K-mer counts as rows in a matrix, since the elements in a column are contiguous in memory.\n\n\n\n\n\n","category":"type"},{"location":"#VectorizedKmers.KmerCountRows","page":"Home","title":"VectorizedKmers.KmerCountRows","text":"KmerCountRows{A, K, T, M} <: AbstractKmerCountVector{A, K, T, M}\n\nA container for K-mer counts, where K-mer counts are stored together as rows in a matrix. This is not as efficient as storing K-mer counts as columns in a matrix, since the elements in a row are not contiguous in memory.\n\n\n\n\n\n","category":"type"},{"location":"#VectorizedKmers.count_kmers!-Union{Tuple{T}, Tuple{K}, Tuple{A}, Tuple{KmerCount{A, K, T}, Vector{<:Integer}}} where {A, K, T}","page":"Home","title":"VectorizedKmers.count_kmers!","text":"count_kmers!(kmer_count, kmers; reset=true)\n\nMutate the counts vector in kmer_count by counting K-mers in kmers. The K-mers in kmers must be represented as integers between 0 and length(kmer_count) - 1.\n\nIf reset is true, the counts vector will be zero-ed before counting.\n\n\n\n\n\n","category":"method"},{"location":"#VectorizedKmers.count_kmers-Union{Tuple{T}, Tuple{K}, Tuple{A}, Tuple{KmerCount{A, K, T}, Vector{<:Integer}}} where {A, K, T}","page":"Home","title":"VectorizedKmers.count_kmers","text":"count_kmers(KmerCount{A, K, T}, kmers; zeros_func=zeros, reset=true)\n\nCreate a new A^K sized vector using zeros_func and count the K-mers in kmers. The K-mers in kmers must be represented as integers between 0 and length(kmer_count) - 1.\n\nIf reset is true, the counts vector will be zero-ed before counting.\n\n\n\n\n\n","category":"method"}]
}
