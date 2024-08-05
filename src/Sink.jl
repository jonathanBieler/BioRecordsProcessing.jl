"""
    AbstractSink
"""
abstract type AbstractSink end

close(sink::T) where T<:AbstractSink = nothing
return_value(sink::T) where T<:AbstractSink = nothing
return_type(sink::T) where T<:AbstractSink = error("Need to implement this method for type $T")
is_paired(sink::T) where T<:AbstractSink = sink.is_paired

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

function Base.show(io::IO, sink::Collect{T}) where T
    print(io, "Collect{$(T)}\n")
end

return_value(sink::Collect) = sink.data
return_type(sink::Collect{T}) where T = Vector{T}

open_writer(sink::Collect, filepath, filename, extension) = begin empty!(sink.data); sink end
open_writer_paired(sink::Collect, filepath1, filename1, filepath2, filename2, extension) = open_writer(sink, filepath1, filename1, extension)

Base.copy(sink::Collect) = deepcopy(sink)

"""
```julia
Writer(record_module::Module, output_directory::String; 
    suffix = "", 
    paired = false, 
    second_in_pair = nothing, 
    extension = nothing, 
    header = nothing
)
```

Write the output of the processing function into a file, the first
argument is the module that owns the `Record` type (e.g `FASTX.FASTA`, `VCF`, ...),
and the second the ouput directory. The filename is determined by the source, to which an optional suffix can be added. 
If the type ouput is different from the type of the output (e.g. SAM to BAM), the extension (".bam") should be specified.
For SAM & BAM a SAM.Header should be provided.

To avoid overwriting existing files, the pipeline will check that the output file is different from the input file.
"""
mutable struct Writer <: AbstractSink
    record_module::Module
    output_directory::String
    suffix::String
    return_value::Union{String, Vector{String}}
    is_paired::Bool
    second_in_pair::Function
    extension::Union{Nothing,String}
    header::Union{Nothing, SAM.Header}
end

function Writer(record_module::Module, output_directory::String;
    suffix = "", paired = false, second_in_pair = nothing, extension = nothing, header = nothing) 

    return_value = (paired || !isnothing(second_in_pair)) ? String[] : ""

    Writer(
        record_module, output_directory, suffix, return_value, paired, 
        isnothing(second_in_pair) ? identity : second_in_pair, extension, header
    )
end

function Base.show(io::IO, sink::Writer)
    print(io, "Writer($(sink.record_module), \"$(sink.output_directory)\")")
    if sink.is_paired
        print(io, "#Paired")
    end
end

record_type(sink::Writer) = sink.record_module

write(sink::Writer, record) = write(sink, record)

# for paired records, single writer
write(sink, record::Tuple{T,T}) where T = begin
    write(sink, record[1])
    write(sink, record[2])
end

# Need to copy because record can be updated in-place
write(sink::Collect, records::Tuple{T,T}) where T = begin
    if sink.is_paired 
        @assert length(records) == 2
    end
    if any(ismutable(r) for r in records)
        push!(sink.data, copy.(records))
    else
        push!(sink.data, records)
    end
end

write(sink::Collect, record::T) where T = push!(sink.data, ismutable(record) ? copy(record) : record)

# for paired records, W can be e.g. BAM.Writer
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

Base.copy(sink::Writer) = Writer(
    sink.record_module, sink.output_directory, sink.suffix, sink.return_value,
    sink.is_paired, sink.second_in_pair, sink.extension, sink.header
)

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

    # override extension when specified
    if !isnothing(sink.extension)
        extension = sink.extension
    end

    out_file = joinpath(sink.output_directory, insert_suffix(filename, extension, sink.suffix))
    
    if sink.return_value isa String
        sink.return_value = out_file
    else
        push!(sink.return_value, out_file)
    end
    if out_file == filepath
        error("Output file $(out_file) is identical to input $(filepath), change output directory or specify a suffix.")
    end

    RecordType = record_type(sink)

    if length(extension) >= 3 && extension[end-2:end] == ".gz"
        writer = RecordType.Writer(GzipCompressorStream(open(out_file, "w")))
    elseif extension ∈ (".bam", ".sam")
        if isnothing(sink.header)
            error("A header should be provided for BAM type")
        end
        writer = RecordType.Writer(BGZFStream(out_file, "w"), sink.header)
    else
        writer = RecordType.Writer(open(out_file, "w"))
    end
    writer
end

function open_writer_paired(sink::Writer, filepath1, filename1, filepath2, filename2, extension)

    if is_paired(sink)
        writer1 = open_writer(sink, filepath1, filename1, extension)
        writer2 = open_writer(sink, filepath2, filename2, extension)
        return writer1, writer2
    else
        writer = open_writer(sink, filepath1, filename1, extension)
        return writer
    end
end