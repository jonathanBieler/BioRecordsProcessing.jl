This contains documentation for `process_` API, it is still supported but is more limited than the new `Pipeline` type, as it can only do file-to-file.

## Usage

To process a directory containing records files the user has to provide an `input_directory`, an `output_directory`,
a [glob](https://github.com/vtjnash/Glob.jl) `pattern` that will select the files to process in the input directory, a module containting the record
type of interrest (e.g. `FASTX.FASTA`, `XAM.BAM`, `VariantCallFormat.VCF`) and finally a `process_function` taking a record as input and returning either a record or `nothing` if the record is meant to be filtered out. If several files are found in a directory they will be processed in parallel using one thread per file.

    process_directory(process_function, ReadType, input_directory, pattern, output_directory)
   
The `process_function` can be provided using a `do` block (see examples bellow). For testing the function it is convinient to limit the number of records to be processed using the optional `max_records` argument.
   
Paired files can be processed together using `process_directory_paired`, an additional function giving the filename of the second file in the pair given the filename of the first has to be provided. A single file can be processed using `process`.

## Examples

### Filter out fasta reads shorter than 50bp :

```julia
using BioRecordsProcessing, FASTX

#test only on first 100 records
BioRecordsProcessing.process_directory(FASTX.FASTA, input_directory, "*.fa", output_directory; max_records=100) do record
    return length(sequence(record)) < 50 ? nothing : record
end
```

### Trim 10bp at the start of paired-end reads :

```julia
using BioRecordsProcessing, FASTX

get_f2 = f1 -> replace(f1, "_1" => "_2")#function to get the name of the second file of the pair

# takes a read and trim it
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

### Read a BAM file and groups reads by template name (experimental)

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