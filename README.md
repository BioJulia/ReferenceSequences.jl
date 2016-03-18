# ReferenceSequences

[![Build Status](https://travis-ci.org/bicycle1885/ReferenceSequences.jl.svg?branch=master)](https://travis-ci.org/bicycle1885/ReferenceSequences.jl)

ReferenceSequences.jl provides a data structure for reference sequences.  It is
common that a reference sequence contains only five kinds of nucleotides
("ACGTN") and the occurrence of 'N' is sparse and clustered. In such a case,
`ReferenceSequence` can compress positions of 'N' and aggressively save memory
space.


```julia
julia> using Bio.Seq

julia> using ReferenceSequences

# create ReferenceSequence from DNASequence of Bio.Seq
julia> seq = ReferenceSequence(dna"ACGT"^5 * dna"N"^10)
30nt Reference Sequence:
 ACGTACGTACGTACGTACGTNNNNNNNNNN

julia> DNASequence(seq)  # round trip
30nt DNA Sequence:
ACGTACGTACGTACGTACGTNNNNNNNNNN

julia> seq[4]  # access an element
T

julia> seq[15:25]  # make a subsequence (copy-free)
11nt Reference Sequence:
 GTACGTNNNNN

julia> seq[1:4] == dna"ACGT"  # comparison
true

```


## TODO

* FASTA parser
