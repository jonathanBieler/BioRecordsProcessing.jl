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