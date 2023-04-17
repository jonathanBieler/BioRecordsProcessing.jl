"""
    AbstractSink
"""
abstract type AbstractSink end

close(sink::T) where T<:AbstractSink = nothing
return_value(sink::T) where T<:AbstractSink = nothing

"""
```julia
Collect(T::DataType)
```

Write the output of the processing function into an `Array` in memory. The type of output has to be provided.
"""
struct Collect{T} <: AbstractSink 
    data::Vector{T}
    Collect(T::DataType) = new{T}(T[])
end

open_writer(sink::Collect, filepath, filename, extension) = begin empty!(sink.data); sink end
write(sink::Collect, record) = push!(sink.data, ismutable(record) ? copy(record) : record) #need to copy because record can be updated in-place
return_value(sink::Collect) = sink.data

"""
```julia
Writer(record_module::Module, output_directory::String; suffix = "")
```

Write the output of the processing function into a file, the first
argument is the module that owns the `Record` type (e.g `FASTX.FASTA`, `VCF`, ...),
and the second the ouput directory. The filename is determined by the source, to which optional suffix can be added.
To avoid overwriting existing files, the pipeline will check that the output file is different from the input file.
"""
mutable struct Writer <: AbstractSink
    record_module::Module
    output_directory::String
    suffix::String
    return_value::String
end
Writer(record_module::Module, output_directory::String; suffix = "") = Writer(record_module, output_directory, suffix, "")

record_type(sink::Writer) = sink.record_module

write(sink::Writer, record) = write(sink, record)
return_value(sink::Writer) = sink.return_value

# for gz files, e.g. .fastq.gz put suffix before the fastq
# name.fastq .gz
# name .fastq 
function insert_suffix(filename, extension, suffix)
    if occursin('.', filename)
        parts = split(filename, '.')
        parts[1] = parts[1] * suffix
        filename = join(parts, '.')
    else
        filename = filename * suffix
    end
    return filename * extension
end

function open_writer(sink::Writer, filepath, filename, extension)

    if filename == ""
        error("Empty filename, make sure the source is able to provide a valid filename.")
    end

    out_file = joinpath(sink.output_directory, insert_suffix(filename, extension, sink.suffix))
    sink.return_value = out_file
    @assert out_file != filepath

    RecordType = record_type(sink)

    if extension == ".gz"
        writer = RecordType.Writer(GzipCompressorStream(open(out_file, "w")))

    elseif extension == ".bam"
        #writer = ReadType.Writer(BGZFStream(out_file, "w"), reader.header)
        error("Need to find a way to handle getting a header for the bam")
    else
        writer = RecordType.Writer(open(out_file, "w"))
    end
    writer
end