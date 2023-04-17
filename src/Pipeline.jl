"""
    AbstractPipeline
"""
abstract type AbstractPipeline end

"""
    Pipeline{So, Si} <: AbstractPipeline where {So <: AbstractSource, Si <: AbstractSink}


"""
mutable struct Pipeline{So, Si} <: AbstractPipeline where {So <: AbstractSource, Si <: AbstractSink}
    source::So
    process_function::Function
    sink::Si
end
Pipeline(source::So, sink::Si) where {So <: AbstractSource, Si <: AbstractSink} = Pipeline(source, identity, sink)

##

# pipeline for a Bio reader
function run(p::Pipeline{<:Reader{File}, Si}; max_records = Inf, verbose = true) where {Si <: AbstractSink}

    filepath = _filename(p.source)
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
        if verbose && (mod(k, 100_000) == 0)
            @info "$(Threads.threadid()), $(basename(filepath)) : Processed $(div(k, 1000))k records..."
        end     
    end
    close(reader)
    close(writer)

    return_value(p.sink)
end

# pipeline for a in memory array
function run(p::Pipeline{<:Buffer, Si}; max_records = Inf, verbose = true) where {Si <: AbstractSink}

    filepath = _filename(p.source)
    filename, extension = splitext(basename(filepath))

    writer = open_writer(p.sink, filepath, filename, extension)
    
    k = 0 
    for record in p.source.data

        out_record = p.process_function(record)
        !isnothing(out_record) && write(writer, out_record)

        k += 1
        k > max_records && break
        if verbose && (mod(k, 100_000) == 0)
            @info "$(Threads.threadid()), $(basename(filepath)) : Processed $(div(k, 1000))k records..."
        end     
    end
    close(writer)

    return_value(p.sink)
end


# pipeline for a folder
function run(p::Pipeline{<:Reader{Directory}, Si}; max_records = Inf, verbose = true) where {Si <: AbstractSink}

    files = p.source.files
    # create Pipelines's and run them in parallel
    
end


##