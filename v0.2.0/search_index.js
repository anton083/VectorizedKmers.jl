var documenterSearchIndex = {"docs":
[{"location":"kmer_int_repr/#Integer-representation-of-k-mers","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"","category":"section"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"This package relies on representing k-mers as integers for indexing, and understanding how it works is recommended (unless you're only using higher-level API stuff).","category":"page"},{"location":"kmer_int_repr/#DNA-sequences","page":"Integer representation of k-mers","title":"DNA sequences","text":"","category":"section"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"For DNA, each non-ambiguous nucleotide is assigned a number between 0 and 3:","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"Nucleotide Base-4 Base-2\nA 0 00\nC 1 01\nG 2 10\nT 3 11","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"note: Note\nAny ordering works, but this one is the one used by BioSequences.jl, and it also has some nice properties, like being in alphabetical order, and that XOR-ing a base with 3 gives you its complement.","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"We could theoretically convert any DNA sequence to an integer, but 64-bit integers limit us to 32-mers.","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"Consider the DNA sequence GATTACA. If we convert it to an integer using the table above, we get 2033010_4 = 10001111000100_2 = 9156_10, so the integer value of GATTACA is 9156. Since Julia uses 1-based indexing, we would add 1 to this value to get the index of GATTACA in the vector.","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"Writing a function for this might look like the following:","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"const DNA_ENCODING_VECTOR = zeros(Int, 127)\n\nfor (i, char) in enumerate(\"ACGT\")\n    DNA_ENCODING_VECTOR[char % Int8] = i - 1\nend\n\nfunction kmer_to_int(kmer::String)\n    kmer_int = 0\n    for char in kmer\n        kmer_int = (kmer_int << 2) + DNA_ENCODING_VECTOR[char % Int8]\n    end\n    kmer_int\nend\n\n# or if you're into code golf:\nf(k;i=0)=[i=4i|(c%Int-1-(c=='C'))&3 for c=k][end]","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"note: Note\nStrings are bad! Don't use strings! Please! They're bad! Very bad! Chars have variable length when part of Strings in Julia, so indexing and taking lengths and stuff is kinda slow! Don't use strings! Again: very bad! Use something like LongDNA{4}, please!","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"This function would not be very efficient in a practical setting (even though we're using super cool bit manipulation stuff), since we're convert each k-mer individually, instead of having some kind of sliding window. Moreover, the function takes the k-mer in the form of a String, which is not ideal. The function should work as intended, though. Let's test it:","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"julia> kmer_to_int(\"GATTACA\")\n9156\n\njulia> f(\"GATTACA\")\n9156","category":"page"},{"location":"kmer_int_repr/#Amino-acid-sequences","page":"Integer representation of k-mers","title":"Amino acid sequences","text":"","category":"section"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"Amino acid sequences are a little more difficult to deal with since there are a lot more of them, and the vectors would grow in size even quicker. However, we can still represent them as integers, but we can't use bit manipulation anymore.","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"BioSequences.jl has 28 amino acids in its AminoAcidAlphabet, so we can represent each amino acid as an integer between 0 and 27.","category":"page"},{"location":"kmer_int_repr/","page":"Integer representation of k-mers","title":"Integer representation of k-mers","text":"","category":"page"},{"location":"kmer_count/","page":"k-mer counting","title":"k-mer counting","text":"CurrentModule = VectorizedKmers\nDocTestSetup = quote\n    using VectorizedKmers\nend","category":"page"},{"location":"kmer_count/#Counting-k-mers-with-the-KmerCount-type","page":"k-mer counting","title":"Counting k-mers with the KmerCount type","text":"","category":"section"},{"location":"kmer_count/","page":"k-mer counting","title":"k-mer counting","text":"The KmerCount type has four type parameters, but you only really need to care about the first two: A, the alphabet size, and k, the k-mer length. So, to count the 6-mers of a DNA sequence, you would use KmerCount{4, 6}. For each of these k-mer counts, memory for a vector of size A^k is allocated, unless a vector type like SparseVector is used. This brings us to the two other type parameters: T, which is the vector element type, and V, which is the type of the actual vector.","category":"page"},{"location":"kmer_count/","page":"k-mer counting","title":"k-mer counting","text":"Let's see it in action! Here we import BioSequences to unlock a method of count_kmers that works on the LongDNA type. In this example, we count the 1-mers of the sequence GATTACA. The result is a KmerCount{4, 1, Int64, Vector{Int64}}, which is a vector of 4 Int64 elements.","category":"page"},{"location":"kmer_count/","page":"k-mer counting","title":"k-mer counting","text":"julia> using BioSequences\n\njulia> count_kmers(dna\"GATTACA\", 1)\n4-element KmerCount{4, 1, Int64, Vector{Int64}}:\n 3\n 1\n 1\n 2","category":"page"},{"location":"kmer_count/","page":"k-mer counting","title":"k-mer counting","text":"We can also set the element type to UInt16 to save memory, at the cost of a smaller maximum count value.","category":"page"},{"location":"kmer_count/","page":"k-mer counting","title":"k-mer counting","text":"julia> count_kmers(dna\"GATTACA\", 1, UInt16)\n4-element KmerCount{4, 1, UInt16, Vector{UInt16}}:\n 0x0003\n 0x0001\n 0x0001\n 0x0002","category":"page"},{"location":"kmer_count/","page":"k-mer counting","title":"k-mer counting","text":"note: Note\nBe careful when using element types with fewer bits, such as UInt16. Even though you might not expect any one k-mer to occur more than 65,535 times, some vector operations such as LinearAlgebra.dot and Distances.sqeuclidean will still use UInt16 when summing up terms, which might lead to integer overflow.","category":"page"},{"location":"kmer_count/","page":"k-mer counting","title":"k-mer counting","text":"The default count_kmers method is a bit different. Under the hood, it uses count_kmers! to modify a KmerCount instance in-place. This is useful when you want to count k-mers in a sequence without allocating a new vector. count_kmers doesn't modify a vector though, so it needs to create a new instance of a given type:","category":"page"},{"location":"kmer_count/","page":"k-mer counting","title":"k-mer counting","text":"julia> kc1 = count_kmers(KmerCount{4, 1, Int64}, [2, 0, 3, 3, 0, 1, 0])\n4-element KmerCount{4, 1, Int64, Vector{Int64}}:\n 3\n 1\n 1\n 2","category":"page"},{"location":"kmer_count/","page":"k-mer counting","title":"k-mer counting","text":"This default method is supposed to be as generic as possible, which is why it takes the k-mers in the form of integers already, but that's not very efficient, since that whole array would have to be allocated. Ideally, k-mers would be procedurally generated in constant memory, as is the case for the KmerCount(::LongDNA, ::Integer) method.","category":"page"},{"location":"kmer_count/","page":"k-mer counting","title":"k-mer counting","text":"To reiterate, the default count_kmers method takes a type as its first argument so that it can create a new instance of that type, which includes creating a new vector of zeros. The count_kmers! method on the other hand, takes an instance of KmerCount as its first argument, and modifies it in-place, which is more flexible.","category":"page"},{"location":"kmer_count/","page":"k-mer counting","title":"k-mer counting","text":"Let's create a KmerCount instance with SparseVector as its vector type. We can do this by passing spzeros to the KmerCount{A, k, T} constructor (it's just zeros by default):","category":"page"},{"location":"kmer_count/","page":"k-mer counting","title":"k-mer counting","text":"julia> using SparseArrays\n\njulia> kc2 = KmerCount{4, 1, Int64}(spzeros)\n4-element KmerCount{4, 1, Int64, SparseVector{Int64, Int64}}:\n 0\n 0\n 0\n 0\n\njulia> count_kmers!(kc2, [2, 0, 3, 3, 0, 1, 0])\n4-element KmerCount{4, 1, Int64, SparseVector{Int64, Int64}}:\n 3\n 1\n 1\n 2\n\njulia> kc1 == kc2\ntrue","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = VectorizedKmers\nDocTestSetup = quote\n    using VectorizedKmers\nend","category":"page"},{"location":"#VectorizedKmers","page":"Home","title":"VectorizedKmers","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"(Image: Latest Release) (Image: MIT license) (Image: Documentation) (Image: Status) (Image: Coverage)","category":"page"},{"location":"","page":"Home","title":"Home","text":"VectorizedKmers.jl is a Julia package for fast k-mer counting of biological sequences. The core idea is that k-mers with an alphabet size A are essentially integers in base A, and can be used as indices in a vector of size A^k to count the corresponding k-mers.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The KmerCount type is a wrapper for AbstractVector, which means that these vector k-mer counts are not limited to Julia's Base.Vector type; other kinds of vectors can be used as well, such as CUDA.CuVector, SparseArrays.SparseVector or even matrix views.","category":"page"},{"location":"","page":"Home","title":"Home","text":"To efficiently group k-mer counts together, the KmerCountVector stores them in a matrix as rows or columns. It can wrap any AbstractMatrix, such as Matrix or CuMatrix, and accessing its elements by index returns a KmerCount wrapped around a view of a row or column of the original matrix.","category":"page"},{"location":"","page":"Home","title":"Home","text":"This data structure can be used to quickly approximate distances between sequences. Most notably, the squared Euclidean distance was used to approximate edit distance in this paper. The dot product has also proven to be a useful metric for comparing correlation between sequences.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"You can install VectorizedKmers in your environment from the Julia REPL by entering pkg mode with ] and then running:","category":"page"},{"location":"","page":"Home","title":"Home","text":"add VectorizedKmers","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [VectorizedKmers]","category":"page"},{"location":"#VectorizedKmers.AbstractKmerCountVector","page":"Home","title":"VectorizedKmers.AbstractKmerCountVector","text":"AbstractKmerCountVector{A, k, T <: Real, M <: AbstractMatrix{T}} <: AbstractVector{KmerCount{A, k, T, V} where {V <: AbstractVector{T}}}\n\nYup... that's indeed an abomination of a type. A container for k-mer counts, where k-mer counts are stored together as rows or columns in a matrix. A is the alphabet size, k is the k-mer size, T is the element type of the counts, and M is the type of the matrix in which the k-mer counts are stored.\n\n\n\n\n\n","category":"type"},{"location":"#VectorizedKmers.KmerCount","page":"Home","title":"VectorizedKmers.KmerCount","text":"KmerCount{A, k, T, V}\n\nA is the alphabet size, k is the k-mer size, and T is the element type of the underlying counts field, which in turn has type V.\n\n\n\n\n\n","category":"type"},{"location":"#VectorizedKmers.KmerCountColumns","page":"Home","title":"VectorizedKmers.KmerCountColumns","text":"KmerCountColumns{A, k, T, M} <: AbstractKmerCountVector{A, k, T, M}\n\nA container for k-mer counts, where k-mer counts are stored together as columns in a matrix. This is more efficient than storing k-mer counts as rows in a matrix, since the elements in a column are contiguous in memory.\n\n\n\n\n\n","category":"type"},{"location":"#VectorizedKmers.KmerCountRows","page":"Home","title":"VectorizedKmers.KmerCountRows","text":"KmerCountRows{A, k, T, M} <: AbstractKmerCountVector{A, k, T, M}\n\nA container for k-mer counts, where k-mer counts are stored together as rows in a matrix. This is not as efficient as storing k-mer counts as columns in a matrix, since the elements in a row are not contiguous in memory.\n\n\n\n\n\n","category":"type"},{"location":"#VectorizedKmers.count_kmers!-Tuple{KmerCount, Vector{<:Integer}}","page":"Home","title":"VectorizedKmers.count_kmers!","text":"count_kmers!(kmer_count, kmers; reset=true)\n\nMutate the counts vector in kmer_count by counting k-mers in kmers. The k-mers in kmers must be represented as integers between 0 and length(kmer_count) - 1.\n\nIf reset is true, the counts vector will be zero-ed before counting.\n\n\n\n\n\n","category":"method"},{"location":"#VectorizedKmers.count_kmers-Union{Tuple{T}, Tuple{k}, Tuple{A}, Tuple{Type{KmerCount{A, k, T, V} where V<:AbstractVector{T}}, Vector{<:Integer}}} where {A, k, T}","page":"Home","title":"VectorizedKmers.count_kmers","text":"count_kmers(KmerCount{A, k, T}, kmers; zeros_func=zeros)\n\nCreate a new A^k sized vector using zeros_func and count the k-mers in kmers. The k-mers in kmers must be represented as integers between 0 and length(kmer_count) - 1.\n\nIf reset is true, the counts vector will be zero-ed before counting.\n\n\n\n\n\n","category":"method"}]
}