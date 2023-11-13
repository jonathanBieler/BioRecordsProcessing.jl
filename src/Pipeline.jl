"""
    AbstractPipeline
"""
abstract type AbstractPipeline end


const OptionalSink = Union{AbstractSink, Nothing}
const OptionalGrouper = Union{RecordGrouper, Nothing}
const FileOrSource = Union{AbstractFileProvider, AbstractSource}

"""
```julia
Pipeline(source, processor, sink)
Pipeline(source, sink)
```

Build a Pipeline, if `processor` is omitted it will default to `identity`.
"""
mutable struct Pipeline{So, G, P, Si} <: AbstractPipeline where {So <: FileOrSource, G <:OptionalGrouper, P <:AbstractProcessor, Si <: OptionalSink}
    source::So
    grouper::G
    processor::P
    sink::Si
end

# no processor
Pipeline(source::So, sink::Si) where {So <: AbstractSource, Si <: OptionalSink} = Pipeline(source, identity, sink)

# no sink
Pipeline(source::So, processor::P) where {So <: AbstractSource, P <: Union{AbstractProcessor, Function}} = Pipeline(source, processor, nothing)

# function as processor without grouper
Pipeline(source::So, process_function::Function, sink::Si) where {So <: AbstractSource, Si <: OptionalSink} =
    Pipeline(source, Processor(process_function), sink)

# function as processor with grouper
Pipeline(source::So, grouper::G, process_function::Function, sink::Si) where {So <: AbstractSource, G <: OptionalGrouper, Si <: OptionalSink} =
    Pipeline(source, grouper, Processor(process_function), sink)

# No grouper with Processor
Pipeline(source::So, processor::P, sink::Si) where {So <: AbstractSource, P <: AbstractProcessor, Si <: OptionalSink} =
    Pipeline(source, nothing, processor, sink)

# ExternalTool take a AbstractFileProvider as source
Pipeline(source::So, processor::ExternalTool) where {So <: AbstractFileProvider} = Pipeline(source, nothing, processor, nothing)

# pipeline for a Bio reader
"""
```julia
run(p::Pipeline; max_records = Inf, verbose = true)
```

Run the pipeline, the processing will stop after `max_records` have been read. Depending on the
sink it will return a path to the output file or an array.
"""
function run(p::Pipeline{<:Reader{File}, <:OptionalGrouper, P, Si}; max_records = Inf, verbose = true) where {Si <: AbstractSink, P <: AbstractProcessor}
    if is_paired(p.source)
        run_paired(p, p.source.file_provider.second_in_pair; max_records=max_records, verbose=verbose)
    else
        run_single(p; max_records=max_records, verbose=verbose)
    end
end

# run single file with no grouper
function run_single(p::Pipeline{<:Reader{File}, <:Nothing, P, Si}; max_records = Inf, verbose = true) where {Si <: AbstractSink, P <: AbstractProcessor}

    filepath = _filename(p.source)
    filename, extension = splitext(basename(filepath))

    reader, RecordType = open_reader(p.source, filepath, filename, extension)
    writer = open_writer(p.sink, filepath, filename, extension)
    
    record = RecordType.Record()
    k = 0 
    while !eof(reader)

        read!(reader, record)
        is_outside_interval(p.source, record) && break
    
        out_record = p.processor(record)
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

# run single file with grouper
function run_single(p::Pipeline{<:Reader{File}, <:RecordGrouper, P, Si}; max_records = Inf, verbose = true) where {Si <: AbstractSink, P <: AbstractProcessor}

    filepath = _filename(p.source)
    filename, extension = splitext(basename(filepath))

    reader, RecordType = open_reader(p.source, filepath, filename, extension)
    writer = open_writer(p.sink, filepath, filename, extension)
    
    record = RecordType.Record()
    k = 0 
    while !eof(reader)

        record, idx = get_record(p.grouper)
        read!(reader, record)
        is_outside_interval(p.source, record) && break
        
        key, isdone = group_record!(p.grouper, record, idx)

        if isdone
            records = (p.grouper.records[i] for i in p.grouper.groups[key])
            out_records = p.processor(records...)
            #for out_record in out_records
            !isnothing(out_records) && write(writer, out_records)

            for i in p.grouper.groups[key]
                free_idx!(p.grouper, i)
            end
            #end
        end
        
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

# pipeline for a Bio reader, paired
function run_paired(p::Pipeline{<:Reader{File}, <:OptionalGrouper, P, Si}, second_in_pair; max_records = Inf, verbose = true) where {Si <: AbstractSink, P <: AbstractProcessor}

    filepath1 = _filename(p.source)
    filename1, extension = splitext(basename(filepath1))

    filename2 = second_in_pair(filename1)
    filepath2 = joinpath(dirname(filepath1), filename2 * extension)

    reader1, RecordType = open_reader(p.source, filepath1, filename1, extension)
    reader2, _          = open_reader(p.source, filepath2, filename2, extension)

    writer = open_writer_paired(p.sink, filepath1, filename1, filepath2, filename2, extension)
    
    record1, record2 = RecordType.Record(), RecordType.Record()
    k = 0 
    while !eof(reader1)

        read!(reader1, record1)
        read!(reader2, record2)

        is_outside_interval(p.source, record1) && is_outside_interval(p.source, record2) && break

        out_record = p.processor(record1, record2)
        !isnothing(out_record) && write(writer, out_record)

        k += 1
        k > max_records && break
        if verbose && (mod(k, 100_000) == 0)
            @info "$(Threads.threadid()), $(basename(filepath1)) : Processed $(div(k, 1000))k records..."
        end     
    end
    close(reader1)
    close(reader2)
    close(writer)

    return_value(p.sink)
end

# pipeline for a in memory array
function run(p::Pipeline{<:Buffer, <:OptionalGrouper, P, Si}; max_records = Inf, verbose = true) where {Si <: AbstractSink, P <: AbstractProcessor}

    filepath = _filename(p.source)
    filename, extension = splitext(basename(filepath))

    writer = open_writer(p.sink, filepath, filename, extension)
    
    k = 0 
    for record in p.source.data

        out_record = p.processor(record)
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
function run(p::Pipeline{<:Reader{Directory}, <:OptionalGrouper, P, Si}; max_records = Inf, verbose = true) where {Si <: AbstractSink, P <: AbstractProcessor}

    files = p.source.file_provider.files
    # create Pipelines's and run them in parallel
    if verbose 
        @info "Processing files:"
        println.(files)
    end

    pipelines = [Pipeline(
        Reader(record_type(p.source), File(
            # File and Directory constructors put identiy as the function when not paired, so we need to check here
            filename; second_in_pair = is_paired(p.source) ? p.source.file_provider.second_in_pair : nothing)
        ),
        deepcopy(p.grouper),
        p.processor,
        isnothing(p.sink) ? p.sink : copy(p.sink)
    ) for filename in files]

    out = Vector{return_type(p.sink)}(undef, length(files))
    Threads.@threads for i in eachindex(pipelines)
        out[i] = run(pipelines[i]; max_records = max_records, verbose = verbose)      
    end
    out
end

## ExternalTool

function run(p::Pipeline{<:File, <:OptionalGrouper, <:ExternalTool, Si}; max_records = Inf, verbose = true) where {Si <: OptionalSink}
    if is_paired(p.source)
        run_paired(p, p.source.second_in_pair; max_records=max_records, verbose=verbose)
    else
        run_single(p; max_records=max_records, verbose=verbose)
    end
end

function run_single(p::Pipeline{<:File, <:OptionalGrouper, <:ExternalTool, Si}; max_records = Inf, verbose = true) where {Si <: OptionalSink}

    filepath = _filename(p.source)
    p.processor(filepath)
end

# pipeline for a Bio reader, paired
function run_paired(p::Pipeline{<:File, <:OptionalGrouper, <:ExternalTool, Si}, second_in_pair; max_records = Inf, verbose = true) where {Si <: OptionalSink}

    filepath1 = _filename(p.source)
    filename1, extension = splitext(basename(filepath1))

    filename2 = second_in_pair(filename1)
    filepath2 = joinpath(dirname(filepath1), filename2 * extension)

    p.processor(filepath1, filepath2)
end

function run(p::Pipeline{<:Directory, <:OptionalGrouper, <:ExternalTool, Si}; max_records = Inf, verbose = true) where {Si <: OptionalSink}

    files = p.source.files
    # create Pipelines's and run them in parallel
    if verbose 
        @info "Processing files:"
        println.(files)
    end

    pipelines = [Pipeline(
        File(
            # File and Directory constructors put identiy as the function when not paired, so we need to check here
            filename; second_in_pair = is_paired(p.source) ? p.source.second_in_pair : nothing
        ),
        deepcopy(p.grouper),
        p.processor,
        isnothing(p.sink) ? p.sink : copy(p.sink)
    ) for filename in files]

    out = Vector{Any}(undef, length(files))
    Threads.@threads for i in eachindex(pipelines)
        out[i] = run(pipelines[i]; max_records = max_records, verbose = verbose)      
    end
    out
end



##
