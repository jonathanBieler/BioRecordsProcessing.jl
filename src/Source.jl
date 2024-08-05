"""
    AbstractFileProvider
"""
abstract type AbstractFileProvider end
is_paired(file_provider::AbstractFileProvider) = file_provider.paired
interval(file_provider::AbstractFileProvider) = nothing

"""
```julia
File(filename; second_in_pair = nothing)
```

For paired files a function taking as argument the filename of the first file in pair and
returning the filename of the second file can be provided. For example one can use `replace`
 or a dictionnary, e.g. `second_in_pair = f1 -> replace(f1, "_1" => "_2")`.
"""
struct File <: AbstractFileProvider
    filename::String
    paired::Bool
    second_in_pair::Function
    interval::Union{Interval, Nothing}
    File(filename; second_in_pair = nothing, interval = nothing) = 
        new(filename, !isnothing(second_in_pair), isnothing(second_in_pair) ? identity : second_in_pair, interval)
end

_filename(file::File) = file.filename
interval(file::File) = file.interval

function Base.show(io::IO, source::File) 
    print(io, "File(\"$(source.filename)\")")
    if source.paired
        print(io, "#Paired")
    end
end


"""
```julia
Directory(directory::String, glob_pattern::String; second_in_pair = nothing)
```

List all files matching the `glob_pattern` (See [Glob.jl](https://github.com/vtjnash/Glob.jl)) in `directory`.
For paired files a function taking as argument the filename of the first file in pair and
returning the filename of the second file can be provided.

```@example
Directory(input_directory, "*.fastq")
```
"""
struct Directory <: AbstractFileProvider
    directory::String
    glob_pattern::String
    files::Vector{String}
    paired::Bool
    second_in_pair::Function
    interval::Union{Interval, Nothing}
    Directory(directory::String, glob_pattern::String; second_in_pair = nothing, interval = nothing) = begin
        files = glob(glob_pattern, directory)
        paired = !isnothing(second_in_pair)
        new(directory, glob_pattern, files, paired, paired ? second_in_pair : identity, interval)
    end
end

function Base.show(io::IO, source::Directory) 
    print(io, "Directory(\"$(source.directory)\")", "\"$(source.glob_pattern)\"")
    if source.paired
        print(io, "#Paired")
    end
end

interval(directory::Directory) = directory.interval

"""
    AbstractSource
"""
abstract type AbstractSource end

# by default source is always in reading interval
is_outside_interval(source::AbstractSource, record) = false

"""
```julia
Reader(record_module::Module, file_provider::F) where {F <: AbstractFileProvider}
```

Read a file or a directory on the disk and produce records of type `record_module.Record`. 
The second argument can be a `File` or a `Directory`.

If a string is passed the second argment will default to `File`.

```@example
Reader(FASTX.FASTA, "test.fa")
Reader(FASTX.FASTA, File("test.fa"))
Reader(FASTX.FASTQ, Directory("data/", "*.fastq"))
```
"""
struct Reader{F,I} <: AbstractSource where {F <: AbstractFileProvider, I}
    record_module::Module
    file_provider::F
    index::I

    Reader(record_module::Module, file_provider::F, index::I) where {F <: AbstractFileProvider, I} = new{F, I}(record_module, file_provider, index)
end
Reader(record_module::Module, file_provider::F; index=nothing) where {F <: AbstractFileProvider} =  Reader(record_module, file_provider, index)
Reader(record_module::Module, filename::String; index=nothing) = Reader(record_module, File(filename); index=index)

record_type(reader::Reader{F}) where {F} = reader.record_module
_filename(reader::Reader) = _filename(reader.file_provider)
is_paired(reader::Reader) = is_paired(reader.file_provider)

function open_reader(source::Reader, filepath, filename, extension)

    RecordType = record_type(source)

    if extension == ".gz"
        reader = RecordType.Reader(GzipDecompressorStream(open(filepath)))

    elseif extension == ".bam"
        
        if !isnothing(source.index)
            reader = RecordType.Reader(open(filepath); index = source.index)
        else
            index_file = filepath * ".bai"
            if !isfile(index_file)
                @warn "Index file not found : $index_file"
                reader = RecordType.Reader(open(filepath))
            else
                reader = RecordType.Reader(open(filepath); index = index_file)
            end
        end

        if !isnothing(interval(source.file_provider))
            seek_region!(reader, interval(source.file_provider))
        end

    else
        reader = RecordType.Reader(open(filepath))
    end
    reader, RecordType
end


function Base.show(io::IO, source::Reader) 
    print(io, "Reader($(source.record_module), $(source.file_provider))")
end


"""
```julia
Buffer(data::T; filename = "")
```

Use the collection `data` as a source of records. An optional filename can be provided when a `Writer`
is used as a sink.
"""
struct Buffer{T} <: AbstractSource 
    data::T
    filename::String
    Buffer(data::T; filename = "") where T = new{T}(data, filename)
end

function Base.show(io::IO, source::Buffer{T}) where T
    print(io, "Buffer{$T}(; filename = \"$(source.filename)\")")
end

_filename(buffer::Buffer) = buffer.filename


## methods for interval selection

function seek_region!(reader::BAM.Reader, region)

    refindex = findfirst(isequal(region.seqname), reader.refseqnames)
    refindex == nothing && throw(ArgumentError("sequence name $(iter.refname) is not found in the header"))
    
    chunks = XAM.BAM.Indexes.overlapchunks(reader.index.index, refindex, region.first:region.last)
    if !isempty(chunks)
        seek(reader, first(chunks).start)
    end
end

function is_outside_interval(source::T, record::BAM.Record) where T<:AbstractSource
    region = interval(source.file_provider)
    if !isnothing(region)
        BAM.position(record) > region.last && return true
    end
    false
end