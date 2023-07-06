"""
    AbstractFileProvider
"""
abstract type AbstractFileProvider end
is_paired(file_provider::AbstractFileProvider) = file_provider.paired

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
Reader(FASTX.FASTA, "test.fa")
Reader(FASTX.FASTA, File("test.fa"))
Reader(FASTX.FASTQ, Directory("data/", "*.fastq"))
```
"""
struct Reader{F} <: AbstractSource where {F <: AbstractFileProvider}
    record_module::Module
    file_provider::F

    Reader(record_module::Module, file_provider::F) where {F <: AbstractFileProvider} = new{F}(record_module, file_provider)
end
Reader(record_module::Module, filename::String) = Reader(record_module, File(filename))

record_type(reader::Reader{F}) where {F} = reader.record_module
_filename(reader::Reader) = _filename(reader.file_provider)
is_paired(reader::Reader) = is_paired(reader.file_provider)

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
