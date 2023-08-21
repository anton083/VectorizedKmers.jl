var documenterSearchIndex = {"docs":
[{"location":"kmer_int_repr/#Integer-representation-of-k-mers","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"","category":"section"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"This package relies on representing k-mers as integers for indexing. This page goes over how that works exactly.","category":"page"},{"location":"kmer_int_repr/#DNA-sequences","page":"Integer representation of k-mers","title":"DNA sequences","text":"","category":"section"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"For DNA, each non-ambiguous nucleotide is assigned a number between 0 and 3:","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"Nucleotide Base-4 Base-2\nA 0 00\nC 1 01\nG 2 10\nT 3 11","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"note: Note\nAny ordering works, but this is the one used by BioSequences.jl. It also has some nice properties, like being in alphabetical order, and that XOR-ing a base with 3 gives you its complement.","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"We could theoretically convert any DNA sequence to an integer, but 64-bit unsigned integers limit us to 32-mers.","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"Consider the DNA sequence GATTACA. If we convert it to an integer using the table above, we get 2033010_4 = 10001111000100_2 = 9156_10, so the integer value of GATTACA is 9156. Since Julia uses 1-based indexing, we would add 1 to this value to get the index of GATTACA in the vector.","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"Writing a function for this might look like the following:","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"const DNA_ENCODING_VECTOR = zeros(UInt, 127)\n\nfor (i, char) in enumerate(\"ACGT\")\n    DNA_ENCODING_VECTOR[char % Int8] = i - 1\nend\n\nfunction kmer_to_int(kmer::String)\n    kmer_int = zero(UInt)\n    for char in kmer\n        kmer_int = (kmer_int << 2) | DNA_ENCODING_VECTOR[char % Int8]\n    end\n    kmer_int\nend\n\n# or if you're into code golf:\nf(k;i=0)=[i=4i|(c%Int-1-(c=='C'))&3 for c=k][end]","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"note: Note\nIf you care about performance, you should use other types like LongDNA{4} instead of String. A lot of operations on strings take linear time, whereas they're constant time for LongDNA{4}. This includes accessing substrings like k-mers.","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"This function would not be very efficient in a practical setting (even though we're using super cool bit manipulation stuff), since we're convert each k-mer individually, instead of having some kind of sliding window. Moreover, the function takes the k-mer in the form of a String, which is not ideal. The function should work as intended, though. Let's test it:","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"julia> kmer_to_int(\"GATTACA\")\n9156\n\njulia> f(\"GATTACA\")\n9156","category":"page"},{"location":"kmer_int_repr/#Amino-acid-sequences","page":"Integer representation of k-mers","title":"Amino acid sequences","text":"","category":"section"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"Amino acid sequences are a little more difficult to deal with since there are a lot more of them, and the vectors would grow in size even quicker. However, we can still represent them as integers, although we can't use bit manipulation anymore.","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"BioSequences.jl has 28 amino acids in its AminoAcidAlphabet, so we can represent each amino acid as an integer between 0 and 27.","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"julia> using BioSequences\n\njulia> length(AminoAcidAlphabet())\n28","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"The AminoAcidAlphabet consists of amino acids with 8 bits each, so we can reinterpret them as 8-bit integers.","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"julia> reinterpret.(Int8, [AA_A, AA_M, AA_I, AA_N, AA_O])\n5-element Vector{Int8}:\n  0\n 12\n  9\n  2\n 20","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"We can use this to convert amino acid sequences to integers, like we did with DNA in the form of Strings, but without the fancy bit manipulation stuff.","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"Let's say we want to convert the amino acid sequence AMINO to an integer. As seen above, the amino acids in the sequence have values of 0, 12, 9, 2, and 20 respectively. Thus, the integer value of the k-mer should be:","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"0 cdot 28^4 + 12 cdot 28^3 + 9 cdot 28^2 + 2 cdot 28^1 + 20 cdot 28^0 = 270556","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"We can write a function for this:","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"function kmer_to_int(kmer::LongAA)\n    kmer_int = zero(UInt)\n    for aa in kmer\n        kmer_int = kmer_int * 28 + reinterpret(UInt8, aa)\n    end\n    kmer_int\nend","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"To test it, we can use the aa\"...\" string macro to create a LongAA instance:","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"julia> kmer_to_int(aa\"AMINO\")\n270556","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"Marvelous!","category":"page"},{"location":"kmer_count_matrix/","page":"Matrices of k-mers counts","title":"Matrices of k-mers counts","text":"CurrentModule = VectorizedKmers\nDocTestSetup = quote\n    using VectorizedKmers\nend","category":"page"},{"location":"kmer_count_matrix/#Matrices-of-k-mer-counts","page":"Matrices of k-mers counts","title":"Matrices of k-mer counts","text":"","category":"section"},{"location":"kmer_count_matrix/","page":"Matrices of k-mers counts","title":"Matrices of k-mers counts","text":"How can we efficiently store multiple k-mer counts of sequences? We could use a normal vector: Base.Vector{<:KmerCountVector}, but remember that KmerCountVector can wrap any AbstractVector, including rows/columns of matrices, which means that we could store the k-mer counts of multiple sequences next to each other in a matrix (all k-mer counts will have a size of A^k). This is exactly  what the AbstractKmerCountMatrix type is for. It has two subtypes: KmerCountColumns and KmerCountRows, which wrap AbstractMatrix types, and store the k-mer counts as columns or rows of the matrix, respectively.","category":"page"},{"location":"kmer_count_matrix/","page":"Matrices of k-mers counts","title":"Matrices of k-mers counts","text":"Let's create an instance of KmerCountColumns, and configure it for storing the 1-mer counts of three DNA sequences. The alphabet size for DNA is 4, so each KmerCountVector will have a size of A^k=4^1=4. We'll initialize it with a matrix of zeros:","category":"page"},{"location":"kmer_count_matrix/","page":"Matrices of k-mers counts","title":"Matrices of k-mers counts","text":"julia> k = 1;\n\njulia> n = 3;\n\njulia> kcc = KmerCountColumns{4, k}(zeros(Int, 4^k, n))\n3-element KmerCountColumns{4, 1, Int64, Matrix{Int64}}:\n KmerCountVector{4, 1, Int64, SubArray{Int64, 1, Matrix{Int64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}}([0, 0, 0, 0])\n KmerCountVector{4, 1, Int64, SubArray{Int64, 1, Matrix{Int64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}}([0, 0, 0, 0])\n KmerCountVector{4, 1, Int64, SubArray{Int64, 1, Matrix{Int64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}}([0, 0, 0, 0])","category":"page"},{"location":"kmer_count_matrix/","page":"Matrices of k-mers counts","title":"Matrices of k-mers counts","text":"Oof! That does not look pretty... Let's break down what's happening here:","category":"page"},{"location":"kmer_count_matrix/","page":"Matrices of k-mers counts","title":"Matrices of k-mers counts","text":"The matrix consists of three different 1-mer counts, stored in columns.\nKmerCountColumns is an AbstractVector, so it gets displayed in the same way that a vector normally would, with one element per row.\nEach element of the KmerCountColumns is a KmerCountVector wrapped around a view of a column of the underlying matrix, hence the SubArray type.","category":"page"},{"location":"kmer_count_matrix/","page":"Matrices of k-mers counts","title":"Matrices of k-mers counts","text":"We can easily access each KmerCountVector by index:","category":"page"},{"location":"kmer_count_matrix/","page":"Matrices of k-mers counts","title":"Matrices of k-mers counts","text":"julia> kcc[1]\n4-element KmerCountVector{4, 1, Int64, SubArray{Int64, 1, Matrix{Int64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}}:\n 0\n 0\n 0\n 0","category":"page"},{"location":"kmer_count_matrix/","page":"Matrices of k-mers counts","title":"Matrices of k-mers counts","text":"We can see the actual underlying matrix by looking at the counts field:","category":"page"},{"location":"kmer_count_matrix/","page":"Matrices of k-mers counts","title":"Matrices of k-mers counts","text":"julia> kcc.counts\n4×3 Matrix{Int64}:\n 0  0  0\n 0  0  0\n 0  0  0\n 0  0  0","category":"page"},{"location":"kmer_count_matrix/","page":"Matrices of k-mers counts","title":"Matrices of k-mers counts","text":"That looks much better!","category":"page"},{"location":"kmer_count_matrix/#Let's-count-some-k-mers","page":"Matrices of k-mers counts","title":"Let's count some k-mers","text":"","category":"section"},{"location":"kmer_count_matrix/","page":"Matrices of k-mers counts","title":"Matrices of k-mers counts","text":"Now that we understand how AbstractKmerCountMatrix works, we can count k-mers of sequences using the count_kmers! function, which takes a KmerCountVector and a sequence as input.","category":"page"},{"location":"kmer_count_matrix/","page":"Matrices of k-mers counts","title":"Matrices of k-mers counts","text":"Here, our KmerCountColumns instance, which we call kcc, is configured to count the 1-mers of DNA. Let's count the 1-mers of the sequence GATTACA in the first column of kcc:","category":"page"},{"location":"kmer_count_matrix/","page":"Matrices of k-mers counts","title":"Matrices of k-mers counts","text":"julia> using BioSequences\n\njulia> count_kmers!(kcc[1], dna\"GATTACA\");\n\njulia> kcc.counts\n4×3 Matrix{Int64}:\n 3  0  0\n 1  0  0\n 1  0  0\n 2  0  0","category":"page"},{"location":"kmer_count_matrix/","page":"Matrices of k-mers counts","title":"Matrices of k-mers counts","text":"Splendid! We can see that the 1-mers A, C, G, and T occur 3, 1, 1, and 2 times, respectively, in the sequence GATTACA.","category":"page"},{"location":"kmer_count_matrix/","page":"Matrices of k-mers counts","title":"Matrices of k-mers counts","text":"Although easy to visualize, 1-mer-counts of short DNA sequences are not very interesting, and are also not very useful...","category":"page"},{"location":"kmer_count_vector/","page":"Vectors of k-mer counts","title":"Vectors of k-mer counts","text":"CurrentModule = VectorizedKmers\nDocTestSetup = quote\n    using VectorizedKmers\nend","category":"page"},{"location":"kmer_count_vector/#Vectors-of-k-mer-counts","page":"Vectors of k-mer counts","title":"Vectors of k-mer counts","text":"","category":"section"},{"location":"kmer_count_vector/","page":"Vectors of k-mer counts","title":"Vectors of k-mer counts","text":"The KmerCountVector type has four type parameters, but you only really need to care about the first two: A, the alphabet size, and k, the k-mer length. So, to count the 6-mers of a DNA sequence, you would use KmerCountVector{4, 6}. For each of these k-mer counts, memory for a vector of size A^k is allocated, unless a vector type like SparseVector is used. This brings us to the two other type parameters: T, which is the vector element type, and V, which is the type of the actual vector.","category":"page"},{"location":"kmer_count_vector/","page":"Vectors of k-mer counts","title":"Vectors of k-mer counts","text":"Let's see it in action! Here we import BioSequences to unlock a method of count_kmers that works on the LongDNA type. In this example, we count the 1-mers of the sequence GATTACA. The result is a KmerCountVector{4, 1, Int64, Vector{Int64}}, which is a vector of 4 Int64 elements.","category":"page"},{"location":"kmer_count_vector/","page":"Vectors of k-mer counts","title":"Vectors of k-mer counts","text":"julia> using BioSequences\n\njulia> count_kmers(dna\"GATTACA\", 1)\n4-element KmerCountVector{4, 1, Int64, Vector{Int64}}:\n 3\n 1\n 1\n 2","category":"page"},{"location":"kmer_count_vector/","page":"Vectors of k-mer counts","title":"Vectors of k-mer counts","text":"We can also set the element type to UInt16 to save memory, at the cost of a smaller maximum count value.","category":"page"},{"location":"kmer_count_vector/","page":"Vectors of k-mer counts","title":"Vectors of k-mer counts","text":"julia> count_kmers(dna\"GATTACA\", 1, UInt16)\n4-element KmerCountVector{4, 1, UInt16, Vector{UInt16}}:\n 0x0003\n 0x0001\n 0x0001\n 0x0002","category":"page"},{"location":"kmer_count_vector/","page":"Vectors of k-mer counts","title":"Vectors of k-mer counts","text":"note: Note\nBe careful when using element types with fewer bits, such as UInt16. Even though you might not expect any one k-mer to occur more than 65,535 times, some vector operations such as LinearAlgebra.dot and Distances.sqeuclidean will still use UInt16 when summing up terms, which might lead to integer overflow.","category":"page"},{"location":"kmer_count_vector/","page":"Vectors of k-mer counts","title":"Vectors of k-mer counts","text":"The default count_kmers method is a bit different. Under the hood, it uses count_kmers! to modify a KmerCountVector instance in-place. This is useful when you want to count k-mers in a sequence without allocating a new vector. count_kmers doesn't modify a vector though, so it needs to create a new instance of a given type:","category":"page"},{"location":"kmer_count_vector/","page":"Vectors of k-mer counts","title":"Vectors of k-mer counts","text":"julia> kc1 = count_kmers(KmerCountVector{4, 1, Int64}, [2, 0, 3, 3, 0, 1, 0])\n4-element KmerCountVector{4, 1, Int64, Vector{Int64}}:\n 3\n 1\n 1\n 2","category":"page"},{"location":"kmer_count_vector/","page":"Vectors of k-mer counts","title":"Vectors of k-mer counts","text":"This default method is supposed to be as generic as possible, which is why it takes the k-mers in the form of integers already, but that's not very efficient, since the entire vector of k-mers would have to be allocated. Ideally, k-mers would be procedurally generated in constant memory, as is the case for the KmerCountVector(::LongDNA, ::Integer) method.","category":"page"},{"location":"kmer_count_vector/","page":"Vectors of k-mer counts","title":"Vectors of k-mer counts","text":"To reiterate, the default count_kmers method takes a type as its first argument so that it can create a new instance of that type, which includes creating a new vector of zeros. The count_kmers! method on the other hand, takes an instance of KmerCountVector as its first argument, and modifies it in-place, which is more flexible.","category":"page"},{"location":"kmer_count_vector/","page":"Vectors of k-mer counts","title":"Vectors of k-mer counts","text":"Let's create a KmerCountVector instance with SparseVector as its vector type. We can do this by passing spzeros to the KmerCountVector{A, k, T} constructor (it's just zeros by default):","category":"page"},{"location":"kmer_count_vector/","page":"Vectors of k-mer counts","title":"Vectors of k-mer counts","text":"julia> using SparseArrays\n\njulia> kc2 = KmerCountVector{4, 1, Int64}(spzeros)\n4-element KmerCountVector{4, 1, Int64, SparseVector{Int64, Int64}}:\n 0\n 0\n 0\n 0\n\njulia> count_kmers!(kc2, [2, 0, 3, 3, 0, 1, 0])\n4-element KmerCountVector{4, 1, Int64, SparseVector{Int64, Int64}}:\n 3\n 1\n 1\n 2\n\njulia> kc1 == kc2\ntrue","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = VectorizedKmers\nDocTestSetup = quote\n    using VectorizedKmers\nend","category":"page"},{"location":"#VectorizedKmers","page":"Home","title":"VectorizedKmers","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"(Image: Latest Release) (Image: MIT license) (Image: Documentation) (Image: Documentation) (Image: Status) (Image: Coverage)","category":"page"},{"location":"","page":"Home","title":"Home","text":"VectorizedKmers.jl is a Julia package for fast k-mer counting of biological sequences. The core idea is that k-mers with an alphabet size A are essentially integers in base A, and can be used as indices in a vector of size A^k to count the corresponding k-mers.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The KmerCountVector type is a wrapper for AbstractVector, which means that these vector k-mer counts are not limited to Julia's Base.Vector type; other kinds of vectors can be used as well, such as CUDA.CuVector, SparseArrays.SparseVector or even matrix views.","category":"page"},{"location":"","page":"Home","title":"Home","text":"To efficiently group k-mer counts together, the AbstractKmerCountMatrix stores them in a matrix as rows or columns. It can wrap any AbstractMatrix, such as Matrix or CuMatrix, and accessing its elements by index returns a KmerCountVector wrapped around a view of a row or column of the original matrix.","category":"page"},{"location":"","page":"Home","title":"Home","text":"This data structure can be used to quickly approximate distances between sequences. Most notably, the squared Euclidean distance was used to approximate edit distance in this paper. The dot product has also proven to be a useful metric for comparing correlation between sequences.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"You can install VectorizedKmers in your environment from the Julia REPL by entering pkg mode with ] and then running:","category":"page"},{"location":"","page":"Home","title":"Home","text":"add VectorizedKmers","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [VectorizedKmers]","category":"page"},{"location":"#VectorizedKmers.AbstractKmerCountMatrix","page":"Home","title":"VectorizedKmers.AbstractKmerCountMatrix","text":"AbstractKmerCountMatrix{A, k, T <: Real, M <: AbstractMatrix{T}} <: AbstractVector{KmerCountVector{A, k, T, V} where {V <: AbstractVector{T}}}\n\nYup... that's indeed an abomination of a type. A container for k-mer counts, where k-mer counts are stored together as rows or columns in a matrix. A is the alphabet size, k is the k-mer size, T is the element type of the counts, and M is the type of the matrix in which the k-mer counts are stored.\n\n\n\n\n\n","category":"type"},{"location":"#VectorizedKmers.KmerCountColumns","page":"Home","title":"VectorizedKmers.KmerCountColumns","text":"KmerCountColumns{A, k, T, M} <: AbstractKmerCountMatrix{A, k, T, M}\n\nA container for k-mer counts, where k-mer counts are stored together as columns in a matrix. This is more efficient than storing k-mer counts as rows in a matrix, since the elements in a column are contiguous in memory.\n\n\n\n\n\n","category":"type"},{"location":"#VectorizedKmers.KmerCountRows","page":"Home","title":"VectorizedKmers.KmerCountRows","text":"KmerCountRows{A, k, T, M} <: AbstractKmerCountMatrix{A, k, T, M}\n\nA container for k-mer counts, where k-mer counts are stored together as rows in a matrix. This is not as efficient as storing k-mer counts as columns in a matrix, since the elements in a row are not contiguous in memory.\n\n\n\n\n\n","category":"type"},{"location":"#VectorizedKmers.KmerCountVector","page":"Home","title":"VectorizedKmers.KmerCountVector","text":"KmerCountVector{A, k, T, V}\n\nA is the alphabet size, k is the k-mer size, and T is the element type of the underlying counts field, which in turn has type V.\n\n\n\n\n\n","category":"type"},{"location":"#VectorizedKmers.count_kmers!-Union{Tuple{k}, Tuple{A}, Tuple{KmerCountVector{A, k, T} where T<:Real, Vector{<:Integer}}} where {A, k}","page":"Home","title":"VectorizedKmers.count_kmers!","text":"count_kmers!(kmer_count, kmers; reset=true)\n\nMutate the counts vector in kmer_count by counting k-mers in kmers. The k-mers in kmers must be represented as integers between 0 and length(kmer_count) - 1.\n\nIf reset is true, the counts vector will be zero-ed before counting.\n\n\n\n\n\n","category":"method"},{"location":"#VectorizedKmers.count_kmers-Union{Tuple{T}, Tuple{k}, Tuple{A}, Tuple{Type{KmerCountVector{A, k, T, V} where V<:AbstractVector{T}}, Vector{<:Integer}}} where {A, k, T}","page":"Home","title":"VectorizedKmers.count_kmers","text":"count_kmers(KmerCountVector{A, k, T}, kmers; zeros_func=zeros)\n\nCreate a new A^k sized vector using zeros_func and count the k-mers in kmers. The k-mers in kmers must be represented as integers between 0 and length(kmer_count) - 1.\n\nIf reset is true, the counts vector will be zero-ed before counting.\n\n\n\n\n\n","category":"method"}]
}
