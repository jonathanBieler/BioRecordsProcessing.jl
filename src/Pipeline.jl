"""
    AbstractFileProvider
"""
abstract type AbstractFileProvider end

struct File <: AbstractFileProvider
    filename::String
end

struct Directory <: AbstractFileProvider
    directory::String
    glob_pattern::String
end

"""
    AbstractSource
"""
abstract type AbstractSource end

struct Reader{F} <: AbstractSource where {F <: AbstractFileProvider}
    record_module::Module
    file_provider::F
end
record_type(reader::Reader{F}) where {F} = reader.record_module

function open_reader(s::AbstractSource, filepath, filename, extension)

    RecordType = record_type(s)

    if extension == ".gz"
        reader = RecordType.Reader(GzipDecompressorStream(open(filepath)))

    elseif extension == ".bam"
        index_file = filepath * ".bai"
        if !isfile(index_file)
            @warn "Index file not found : $index_file"
            reader = RecordType.Reader(open(filepath))
        else
            reader = RecordType.Reader(open(filepath); index = index_file)
        end

    else
        reader = RecordType.Reader(open(filepath))
    end
    reader, RecordType
end


"""
    AbstractSink
"""
abstract type AbstractSink end

close(sink::T) where T<:AbstractSink = nothing
return_value(sink::T) where T<:AbstractSink = nothing

struct Collect{T} <: AbstractSink 
    data::Vector{T}
    Collect(T::DataType) = new{T}(T[])
end

open_writer(sink::Collect, filepath, filename, extension) = sink
write(sink::Collect, record) = push!(sink.data, record)
return_value(sink::Collect) = sink.data

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

"""
    AbstractPipeline
"""
abstract type AbstractPipeline end

mutable struct Pipeline{So, Si} <: AbstractPipeline where {So <: AbstractSource, Si <: AbstractSink}
    source::So
    process_function::Function
    sink::Si
end

##

function run(p::Pipeline{<:Reader, Si}; max_records = Inf) where {Si <: AbstractSink}

    filepath = p.source.file_provider.filename
    filename, extension = splitext(basename(filepath))

    reader, RecordType = open_reader(p.source, filepath, filename, extension)
    writer = open_writer(p.sink, filepath, filename, extension)
    
    record = RecordType.Record()
    k = 0 
    while !eof(reader)

        read!(reader, record)
        out_record = p.process_function(record)
        !isnothing(out_record) && write(writer, out_record)

        k += 1
        k > max_records && break
        if mod(k, 100_000) == 0
            @info "$(Threads.threadid()), $(basename(file)) : Processed $(div(k, 1000))k records..."
        end     
    end
    close(reader)
    close(writer)

    return_value(p.sink)
end




##