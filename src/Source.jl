"""
    AbstractFileProvider
"""
abstract type AbstractFileProvider end


"""
    File <: AbstractFileProvider

File(filename)
"""
struct File <: AbstractFileProvider
    filename::String
end

"""
    Directory <: AbstractFileProvider

Directory(directory, glob_pattern)
"""
struct Directory <: AbstractFileProvider
    directory::String
    glob_pattern::String
    files::Vector{String}
    Directory(directory::String, glob_pattern::String) = begin
        files = glob(glob_pattern, directory)        
        new(directory, glob_pattern, files)
    end
end

"""
    AbstractSource
"""
abstract type AbstractSource end

"""
```julia
Reader(record_module::Module, file_provider::F) where {F <: AbstractFileProvider}
```

Read a file or a directory on the disk and produce records of type `record_module.Record`. 
The second argument can be a `File` or a `Directory`.

If a string is passed the second argment will default to `File`.

```@example
Reader(FASTX.FASTQ, "test.fa")
```
"""
struct Reader{F} <: AbstractSource where {F <: AbstractFileProvider}
    record_module::Module
    file_provider::F

    Reader(record_module::Module, file_provider::F) where {F <: AbstractFileProvider} = new{F}(record_module, file_provider)
end
Reader(record_module::Module, filename::String) = Reader(record_module, File(filename))

record_type(reader::Reader{F}) where {F} = reader.record_module
_filename(reader::Reader) = reader.file_provider.filename

function open_reader(s::Reader, filepath, filename, extension)

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
```julia
Buffer(data::Vector{T}; filename = "")
```

Use the array `data` as a source of records. An optional filename can be provided when a `Writer`
is used as a source.
"""
struct Buffer{T} <: AbstractSource 
    data::Vector{T}
    filename::String
    Buffer(data::Vector{T}; filename = "") where T = new{T}(data, filename)
end

_filename(buffer::Buffer) = buffer.filename
