"""
    AbstractSink
"""
abstract type AbstractSink end

close(sink::T) where T<:AbstractSink = nothing
return_value(sink::T) where T<:AbstractSink = nothing
return_type(sink::T) where T<:AbstractSink = error("Need to implement this method for type $T")

"""
```julia
Collect(T::DataType; paired=false)
```

Write the output of the processing function into an vector in memory. The type of output has to be provided.
For paired files the option paired need to be set to `true`, the output will then consists of a vector of tuples.
"""
struct Collect{T} <: AbstractSink 
    data::Vector{T}
    is_paired::Bool
    Collect(T::DataType; paired=false) = begin
        T = paired ? Tuple{T,T} : T
        new{T}(T[], paired)
    end
end

#need to copy because record can be updated in-place
write(sink::Collect, record) = begin
    if sink.is_paired
        @assert length(record) == 2
        push!(sink.data, ismutable(record[1]) ? (copy(record[1]), copy(record[2])) : record)
    else
        push!(sink.data, ismutable(record) ? copy(record) : record) 
    end
end
return_value(sink::Collect) = sink.data
return_type(sink::Collect{T}) where T = Vector{T}

open_writer(sink::Collect, filepath, filename, extension) = begin empty!(sink.data); sink end
open_writer_paired(sink::Collect, filepath1, filename1, filepath2, filename2, extension) = open_writer(sink, filepath1, filename1, extension)

Base.copy(sink::Collect) = deepcopy(sink)

"""
```julia
Writer(record_module::Module, output_directory::String; suffix = "")
```

Write the output of the processing function into a file, the first
argument is the module that owns the `Record` type (e.g `FASTX.FASTA`, `VCF`, ...),
and the second the ouput directory. The filename is determined by the source, to which an optional suffix can be added.
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

# for paired records
write(sink::Tuple{W, W}, record::Tuple{T,T}) where {W, T} = begin
    write(sink[1], record[1])
    write(sink[2], record[2])
end

close(sink::Tuple{W, W}) where {W} =  begin 
    close(sink[1])
    close(sink[2]) 
end

return_value(sink::Writer) = sink.return_value
return_type(sink::Writer) = String

Base.copy(sink::Writer) = Writer(sink.record_module, sink.output_directory, sink.suffix, sink.return_value)

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

function open_writer_paired(sink::Writer, filepath1, filename1, filepath2, filename2, extension)

    writer1 = open_writer(sink, filepath1, filename1, extension)
    writer2 = open_writer(sink, filepath2, filename2, extension)

    writer1, writer2
end