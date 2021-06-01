# BioRecordsProcessing

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jonathanBieler.github.io/BioRecordsProcessing.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jonathanBieler.github.io/BioRecordsProcessing.jl/dev)
[![Build Status](https://travis-ci.com/jonathanBieler/BioRecordsProcessing.jl.svg?branch=master)](https://travis-ci.com/jonathanBieler/BioRecordsProcessing.jl)
[![Coverage](https://codecov.io/gh/jonathanBieler/BioRecordsProcessing.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jonathanBieler/BioRecordsProcessing.jl)

Easily process files containing biological records. If several files are found in a directory
they will be processed in parallel (one thread per file).

## Examples

### Filter fasta reads shorter than 50bp :

```julia
using BioRecordsProcessing, FASTX

#test only on first 100 records
BioRecordsProcessing.process_directory(FASTX.FASTA, input_directory, "*.fa", output_directory; max_records=100) do record
    return FASTX.FASTA.seqlen(record) < 50 ? nothing : record
end
```

### Trim 10bp at the start of paired-end reads :

```julia
using BioRecordsProcessing, FASTX

get_f2 = f1 -> replace(f1, "_1" => "_2")#function to get the name of the file containing second reads

trim_record(r, trim) =  begin 
    s = FASTQ.sequence(r)
    sel = min(length(s), trim):length(s)
    return FASTQ.Record(FASTQ.identifier(r), FASTQ.description(r), s[sel], FASTQ.quality(r)[sel])
end

BioRecordsProcessing.process_directory_paired(FASTX.FASTQ, input_directory, "*_1.fastq", get_f2, output_directory) do r1,r2
    r1, r2 = trim_record(r1, 10), trim_record(r2, 10)
    return r1, r2
end
```

### Read a BAM file and groups reads by template name

```julia
using BioRecordsProcessing, XAM

grouper = BioRecordsProcessing.RecordGrouper{XAM.BAM.Record, String}(
    r -> XAM.BAM.tempname(r),           # key used to group reads 
    records -> length(records) == 2     # a group is released when this condition is met
)

BioRecordsProcessing.process(XAM.BAM, grouper, bam_file, output_directory) do r1,r2
    # ... 
    return (r1, r2)
end
```