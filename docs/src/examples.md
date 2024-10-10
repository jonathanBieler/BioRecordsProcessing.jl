```@meta
CurrentModule = BioRecordsProcessing

DocTestSetup = quote
    using BioRecordsProcessing, FASTX, BioSequences
    
    dir = mktempdir()
    filepath = joinpath(dir, "test_1.fa")

    seq = [
        "CTTGGCATACTCAAACTCTT",
        "TGGCATACTCACTAACTCTT",
    ]
    writer = open(FASTA.Writer, filepath)
    for i=1:2
        write(writer, FASTA.Record("seq$i", seq[i]))
    end
    close(writer)

    filepath2 = joinpath(dir, "test_2.fa")
    seq = [
        "GCAAACTCTTCTTGGCATAC",
        "ATACTCAAACTCTTCTTGGC",
    ]
    writer = open(FASTA.Writer, filepath2)
    for i=1:2
        write(writer, FASTA.Record("seq$i", seq[i]))
    end
    close(writer)
end

# filter out temporary folder
DocTestFilters = [
    r"Process\(.*",
    r"/var/folders.*",
    r" \"/var/folders.*",
    r"/tmp/.*",
    r" \"/tmp/.*",
]
```

# Examples

All examples use `FASTA` files but should work with `FASTQ`, compressed `.gz` files, `VCF` from 
`VariantCallFormat.jl` and `XAM.jl` types. More examples can be found in the tests.

## Reading a FASTA file into memory

```jldoctest
using BioRecordsProcessing, FASTX, BioSequences

# the file contains two 20bp reads
p = Pipeline(
    Reader(FASTX.FASTA, File(filepath)),
    record -> begin
        sequence(LongDNA{4}, record)
    end,
    Collect(LongDNA{4}),
)
run(p)

# output
2-element Vector{LongSequence{DNAAlphabet{4}}}:
 CTTGGCATACTCAAACTCTT
 TGGCATACTCACTAACTCTT
```

## Transforming a FASTA file

```jldoctest
using BioRecordsProcessing, FASTX, BioSequences

# the file contains two 20bp reads, trim first 10bp
p = Pipeline(
    Reader(FASTX.FASTA, File(filepath)),
    record -> begin
        seq = sequence(LongDNA{4}, record)
        FASTA.Record(FASTA.identifier(record), seq[10:end])
    end,
    Writer(FASTX.FASTA, dir; suffix = ".trimmed"),
)
out = run(p)
run(`head $out`)

# output
>seq1
CTCAAACTCTT
>seq2
CACTAACTCTT
Process(`head /var/folders/8g/xj7pzy251n53px06l17vr0_00000gr/T/jl_mL4pM7/test_1.trimmed.fa`, ProcessExited(0))
```
## Reading a pair of FASTA file

```jldoctest
using BioRecordsProcessing, FASTX, BioSequences

# first in pair is named "_1.fasta", second "_2.fasta"
p = Pipeline(
    Reader(FASTX.FASTA, File(filepath; second_in_pair = x -> replace(x, "_1" => "_2"))),
    (r1, r2) -> begin
        sequence(LongDNA{4}, r1), sequence(LongDNA{4}, r2)
    end,
    Collect(LongDNA{4}; paired = true),
)
run(p)

# output
2-element Vector{Tuple{LongSequence{DNAAlphabet{4}}, LongSequence{DNAAlphabet{4}}}}:
 (CTTGGCATACTCAAACTCTT, GCAAACTCTTCTTGGCATAC)
 (TGGCATACTCACTAACTCTT, ATACTCAAACTCTTCTTGGC)
```

## Processing all files in a directory

```jldoctest
using BioRecordsProcessing, FASTX, BioSequences


p = Pipeline(
    Reader(FASTX.FASTA, Directory(dir, "*.fa")),
    record -> begin
        seq = sequence(LongDNA{4}, record)
        FASTA.Record(FASTA.identifier(record), seq[10:end])
    end,
    Writer(FASTX.FASTA, dir; suffix = ".trimmed"),
)
out = run(p; verbose = false)
basename.(out)# run returns the path to output files

# output
2-element Vector{String}:
 "test_1.trimmed.fa"
 "test_2.trimmed.fa"
```

## Write sequences in memory into a file

```jldoctest
using BioRecordsProcessing, FASTX, BioSequences

data = [FASTA.Record("seq1", dna"ATGC")]

p = Pipeline(
    Buffer(data; filename = "test.fa"),
    Writer(FASTX.FASTA, dir),
)
out = run(p)
run(`head $out`)

# output
>seq1
ATGC
Process(`head /var/folders/8g/xj7pzy251n53px06l17vr0_00000gr/T/jl_NSdfEq/test.fa`, ProcessExited(0))
```

# Cookbook

## Loading paired-end reads from a BAM file in a genomic interval

In general records can be grouped using a `RecordGrouper`. For pair-end BAMs `BAMPairedReadGrouper`
is provided. It will group the two first reads (primary alignments only) with the same
read name it encounters and release them for processing as a pair:

```julia
using BioRecordsProcessing, XAM, GenomicFeatures

region = Interval("9", 22331023, 24542023)

p = Pipeline(
    Reader(XAM.BAM, File(bamfile; interval = region)),
    BioRecordsProcessing.BAMPairedReadGrouper(),
    (r1,r2) -> begin 
        (r1,r2)
    end,
    Collect(Tuple{XAM.BAM.Record, XAM.BAM.Record}),
)
out = run(p)
```

## Generate artificial FASTQ reads

A generator can be passed as input in a Buffer : 

```jldoctest
using BioRecordsProcessing, FASTX, BioSequences

generate_read(i) = FASTA.Record("seq$i", randdnaseq(150))

p = Pipeline(
    Buffer(generate_read(i) for i in 1:10; filename = "test.fa"),
    Writer(FASTX.FASTA, dir),
)
out = run(p)

# output
"/tmp/jl_0mSqQJ/test.fastq"
```

Paired-end reads can also be generated :

```jldoctest
generate_read(i) = FASTQ.Record("seq$i", randdnaseq(150), fill(UInt8(40), 150))
generate_reads(i) = (generate_read(i), generate_read(i))

p = Pipeline(
    Buffer(generate_reads(i) for i in 1:10; filename = "test_R1.fastq.gz"),
    Writer(FASTX.FASTQ, dir; paired = true, second_in_pair = x -> replace(x, "_R1" => "_R2")),
)
out = run(p)

# output
2-element Vector{String}:
 "/tmp/jl_srnqiA/test_R1.fastq.gz"
 "/tmp/jl_srnqiA/test_R2.fastq.gz"
```

## Reading a VCF in a DataFrame

Note : This probably wont' work because of compatibilty issues, see :
https://github.com/rasmushenningsson/VariantCallFormat.jl/issues/5

```julia
using VariantCallFormat

get_depth(r) = parse(Int, VCF.genotype(r, 1, "DP"))
get_vaf(r) = parse(Float64, VCF.genotype(r, 1, "VAF"))

T = NamedTuple{(:chr, :pos, :ref, :alt, :depth, :vaf, :quality), Tuple{String, Int64, String, String, Int64, String, Float64}}

p = Pipeline(
    Reader(VCF, File(filepath)),
    r -> begin
        VCF.filter(r) !=  ["PASS"] && return nothing #filter variants that are not PASS

        (chr = VCF.chrom(r), pos = VCF.pos(r), ref = VCF.ref(r), alt = VCF.alt(r)[1], 
        depth=get_depth(r), vaf=VCF.genotype(r, 1, "VAF"), quality = VCF.qual(r)
        )
    end,
    Collect(T),
)

df = run(p) |> DataFrame
```

## Aligning a FASTQ file to the reference genome :

See this blog post : https://jonathanbieler.github.io/blog/fastq2cnv/
