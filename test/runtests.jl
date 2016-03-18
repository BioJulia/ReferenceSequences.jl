using ReferenceSequences
using Bio.Seq

if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

@testset "Conversion" begin
    @test isa(convert(ReferenceSequence, dna""), ReferenceSequence)
    @test isa(convert(ReferenceSequence, dna"ACGTN"), ReferenceSequence)
    @test isa(convert(ReferenceSequence, dna"ACGTN"^50), ReferenceSequence)
    @test_throws ArgumentError convert(ReferenceSequence, dna"ACGRT")

    @test isa(convert(ReferenceSequence, ""), ReferenceSequence)
    @test isa(convert(ReferenceSequence, "ACGTN"), ReferenceSequence)

    refseq = ReferenceSequence(dna"ACGTN")
    @test isa(convert(DNASequence, refseq), DNASequence)
    @test isa(convert(ASCIIString, refseq), ASCIIString)
end

@testset "Basic Operations" begin
    seq = ReferenceSequence(dna"")
    @test length(seq) === endof(seq) === 0
    @test isempty(seq)
    @test seq == dna""
    @test cmp(seq, seq) == 0
    @test cmp(seq, dna"") == 0
    @test cmp(dna"", seq) == 0
    @test !isless(seq, seq)
    @test isless(seq, dna"A")
    @test_throws BoundsError seq[-1]
    @test_throws BoundsError seq[0]
    @test_throws BoundsError seq[1]
    @test_throws BoundsError seq[2:3]
    @test ReferenceSequence("") == seq

    seq = ReferenceSequence(dna"ACGTN")
    @test length(seq) === endof(seq) === 5
    @test !isempty(seq)
    @test seq[1] === DNA_A
    @test seq[2] === DNA_C
    @test seq[3] === DNA_G
    @test seq[4] === DNA_T
    @test seq[5] === DNA_N
    @test seq[2:3] == dna"CG"
    @test seq[3:5] == dna"GTN"
    @test seq[5:4] == dna""
    @test_throws BoundsError seq[0]
    @test_throws BoundsError seq[6]
    @test_throws BoundsError seq[0:2]
    @test_throws BoundsError seq[5:6]
    @test_throws MethodError seq[1] = DNA_A
    @test collect(seq) == [DNA_A, DNA_C, DNA_G, DNA_T, DNA_N]
    @test seq == seq
    @test seq == dna"ACGTN"
    @test dna"ACGTN" == seq
    @test cmp(seq, seq) == 0
    @test ReferenceSequence("ACGTN") == seq
end

@testset "Long Sequences" begin
    seq = ReferenceSequence(dna"ACGT"^10000)
    @test length(seq) == 4 * 10000
    @test seq == dna"ACGT"^10000
    @test_throws BoundsError seq[0]
    @test_throws BoundsError seq[4 * 10000 + 1]
    @test ReferenceSequence("ACGT"^10000) == seq

    seq = ReferenceSequence(
        dna"NNNN"^1000 *
        dna"ACGT"^1000 *
        dna"NNNN"^1000 *
        dna"TGCA"^1000 *
        dna"NNNN"^1000
    )
    @test length(seq) == 4 * 1000 * 5
    @test seq == (
        dna"NNNN"^1000 *
        dna"ACGT"^1000 *
        dna"NNNN"^1000 *
        dna"TGCA"^1000 *
        dna"NNNN"^1000)
    @test ReferenceSequence(
        "NNNN"^1000 *
        "ACGT"^1000 *
        "NNNN"^1000 *
        "TGCA"^1000 *
        "NNNN"^1000) == seq
end
