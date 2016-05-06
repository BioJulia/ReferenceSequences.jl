"""
Reference Sequence

Reference sequence is a sequence of A/C/G/T/N. In the internals, it compresses
`N` positions and consumes less than three bits per base. Unlike `BioSequence`
in Bio.jl, reference sequences are immutable and hence no modifyting
operators are provided.
"""
immutable ReferenceSequence <: Sequence
    data::Vector{UInt64}
    nmask::NMask
    part::UnitRange{Int}
end

function ReferenceSequence{T<:Integer}(seq::ReferenceSequence, part::UnitRange{T})
    ReferenceSequence(seq.data, seq.nmask, part)
end

function Base.convert{A<:DNAAlphabet}(::Type{ReferenceSequence}, seq::BioSequence{A})
    data = Vector{UInt64}(cld(length(seq), 32))
    nmask = falses(length(seq))
    i = 1
    for j in 1:endof(data)
        x = UInt64(0)
        r = 0
        while r < 64 && i ≤ endof(seq)
            nt = seq[i]
            if nt ≤ DNA_T
                x |= convert(UInt64, nt) << r
            elseif nt == DNA_N
                nmask[i] = true
            else
                throw(ArgumentError("found a invalid symbol $(seq[i]) ∉ {A,C,G,T,N} at $i"))
            end
            i += 1
            r += 2
        end
        data[j] = x
    end
    return ReferenceSequence(data, NMask(nmask), 1:length(seq))
end

function Base.convert(::Type{ReferenceSequence}, seq::AbstractString)
    return ReferenceSequence(DNASequence(seq))
end

function Base.convert(::Type{DNASequence}, seq::ReferenceSequence)
    bioseq = DNASequence(length(seq))
    for i in 1:endof(seq)
        bioseq[i] = seq[i]
    end
    return bioseq
end

function Base.convert{S<:AbstractString}(::Type{S}, seq::ReferenceSequence)
    return S([Char(nt) for nt in seq])
end

Base.endof(seq::ReferenceSequence) = length(seq.part)
Base.length(seq::ReferenceSequence) = length(seq.part)
Base.isempty(seq::ReferenceSequence) = length(seq) == 0

function Base.checkbounds(seq::ReferenceSequence, i::Integer)
    if 1 ≤ i ≤ endof(seq)
        return true
    end
    throw(BoundsError(i))
end

function Base.checkbounds(seq::ReferenceSequence, part::UnitRange)
    if 1 ≤ part.start && part.stop ≤ endof(seq)
        return true
    end
    throw(BoundsError(part))
end

@inline function Base.getindex(seq::ReferenceSequence, i::Integer)
    checkbounds(seq, i)
    return unsafe_getindex(seq, i)
end

@inline function unsafe_getindex(seq::ReferenceSequence, i::Integer)
    j = Int(i) + seq.part.start - 2
    d, r = divrem32(j)
    @inbounds begin
        if seq.nmask[j+1]
            return DNA_N
        else
            return DNANucleotide((seq.data[d+1] >> 2r) & 0b11)
        end
    end
end

function Base.getindex{T<:Integer}(seq::ReferenceSequence, part::UnitRange{T})
    checkbounds(seq, part)
    return ReferenceSequence(seq, part)
end

divrem32(i) = i >> 5, i & 0b11111

Base.eltype(::Type{ReferenceSequence}) = DNANucleotide

# iterator
Base.start(seq::ReferenceSequence) = 1
Base.done(seq::ReferenceSequence, i) = i > endof(seq)
Base.next(seq::ReferenceSequence, i) = unsafe_getindex(seq, i), i + 1


# Comparison
# ----------

function Base.(:(==))(seq1::ReferenceSequence, seq2::ReferenceSequence)
    if length(seq1) != length(seq2)
        return false
    end
    for (x, y) in zip(seq1, seq2)
        if x != y
            return false
        end
    end
    return true
end

function Base.(:(==)){A<:DNAAlphabet}(seq1::ReferenceSequence, seq2::BioSequence{A})
    if length(seq1) != length(seq2)
        return false
    end
    for (x, y) in zip(seq1, seq2)
        if x != y
            return false
        end
    end
    return true
end

Base.(:(==)){A<:DNAAlphabet}(seq1::BioSequence{A}, seq2::ReferenceSequence) =
    seq2 == seq1

function Base.cmp(seq1::ReferenceSequence, seq2::ReferenceSequence)
    for (x, y) in zip(seq1, seq2)
        c = cmp(x, y)
        if c != 0
            return c
        end
    end
    return cmp(length(seq1), length(seq2))
end

function Base.cmp{A<:DNAAlphabet}(seq1::ReferenceSequence, seq2::BioSequence{A})
    for (x, y) in zip(seq1, seq2)
        c = cmp(x, y)
        if c != 0
            return c
        end
    end
    return cmp(length(seq1), length(seq2))
end

Base.cmp{A<:DNAAlphabet}(seq1::BioSequence{A}, seq2::ReferenceSequence) =
    -cmp(seq2, seq1)

Base.isless(seq1::ReferenceSequence, seq2::ReferenceSequence) =
    cmp(seq1, seq2) == -1
Base.isless{A<:DNAAlphabet}(seq1::ReferenceSequence, seq2::BioSequence{A}) =
    cmp(seq1, seq2) == -1
Base.isless{A<:DNAAlphabet}(seq1::BioSequence{A}, seq2::ReferenceSequence) =
    cmp(seq1, seq2) == -1


# Finders
# -------

function Base.findnext(seq::ReferenceSequence, val, start::Integer)
    checkbounds(seq, start)
    v = convert(DNANucleotide, val)
    if v == DNA_N
        return findnextn(seq.nmask, start)
    else
        for i in Int(start):endof(seq)
            x = unsafe_getindex(seq, i)
            if x == val
                return i
            end
        end
    end
    return 0
end

function Base.findprev(seq::ReferenceSequence, val, start::Integer)
    checkbounds(seq, start)
    v = convert(DNANucleotide, val)
    if v == DNA_N
        return findprevn(seq.nmask, start)
    else
        for i in Int(start):-1:1
            x = unsafe_getindex(seq, i)
            if x == val
                return i
            end
        end
    end
    return 0
end

Base.findfirst(seq::ReferenceSequence, val) = findnext(seq, val, 1)
Base.findlast(seq::ReferenceSequence, val)  = findprev(seq, val, endof(seq))

# Printers
# --------

Base.summary(seq::ReferenceSequence) = string(length(seq), "nt Reference Sequence")

function Base.show(io::IO, seq::ReferenceSequence)
    println(io, summary(seq), ':')
    showcompact(io, seq)
end

function Base.showcompact(io::IO, seq::ReferenceSequence)
    width = Base.tty_size()[2]
    if length(seq) > width
        h = div(width, 2)
        for i in 1:h-2
            print(io, seq[i])
        end
        print(io, '…')
        for i in endof(seq)-h+1:endof(seq)
            print(io, seq[i])
        end
    else
        for i in 1:endof(seq)
            print(io, seq[i])
        end
    end
end

function Base.print(io::IO, seq::ReferenceSequence)
    for nt in seq
        print(io, nt)
    end
end

# TODO: implement following functions for incremantal construction
# Base.push!(seq::ReferenceSequence, x)
# Base.append!(seq::ReferenceSequence, other)
