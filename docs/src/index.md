```@meta
CurrentModule = BioRecordsProcessing

DocTestSetup = quote
    using BioRecordsProcessing, FASTX, BioSequences
    
    dir = mktempdir()
    filepath = joinpath(dir, "test.fa")

    seq = [
        "CTTGGCATACTCAAACTCTT",
        "CTTGGCATACTCAAACTCTT",
    ]
    writer = open(FASTA.Writer, filepath)
    for i=1:2
        write(writer, FASTA.Record("seq$i", seq[i]))
    end
    close(writer)
end
```

# BioRecordsProcessing


In BioRecordsProcessing records are processed using a `Pipeline` that is constructed by taking a source (producing records),
a user-defined function to process the records and a sink that will store the output of the processing function. The pipeline can then be run.

In this example a FASTA file is read from the disk, the sequence is extracted from the records and collected in an array :

```jldoctest
using BioRecordsProcessing, FASTX, BioSequences

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
 CTTGGCATACTCAAACTCTT
```

By using different combinations of source and sink, and with user defined processing function, this allows
to handle many common cases of biological records processing.

### Conventions

- If the processing function returns `nothing` the record will not be written to the sink, allowing to filter out records.
- When writing a file to the disk the sink will get the filename from the source, so a source need to have a filename provided in this case. 

### Sources

```@docs
BioRecordsProcessing.Reader
```

```@docs
BioRecordsProcessing.Buffer
```

### Sinks

```@docs
Writer
```

```@docs
Collect
```

### Examples
