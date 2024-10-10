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
    hasbeenrun::Bool
end


Pipeline(source::So, grouper::G, processor::P, sink::Si) where {So <: FileOrSource, G <:OptionalGrouper, P <:AbstractProcessor, Si <: OptionalSink} =
    Pipeline(source, grouper, processor, sink, false)

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

function Base.show(io::IO, p::Pipeline) 
    print(io, "Pipeline:\n")
    print(io, "  $(p.source)\n")
    print(io, "     ↓ \n")
    print(io, "  $(p.sink)\n")
end

function log_progress(k, verbose, filepath)
    if verbose && (mod(k, 100_000) == 0)
        @info "$(Threads.threadid()), $(basename(filepath)) : Processed $(div(k, 1000))k records..."
    end
end

# pipeline for a Bio reader
"""
```julia
run(p::Pipeline; max_records = Inf, verbose = true)
```

Run the pipeline, the processing will stop after `max_records` have been read. Depending on the
sink it will return a path to the output file or an array.
"""
function run(p::Pipeline{<:Reader{File{I}}, <:OptionalGrouper, P, Si}; max_records = Inf, verbose = true) where {Si <: AbstractSink, P <: AbstractProcessor, I}
    p.hasbeenrun && error("Running the same pipeline twice is currently not supported.")
    p.hasbeenrun
    if is_paired(p.source)
        run_paired(p, p.source.file_provider.second_in_pair; max_records=max_records, verbose=verbose)
    else
        run_single(p; max_records=max_records, verbose=verbose)
    end
end

# run single file with no grouper
function run_single(p::Pipeline{<:Reader{File{I}}, <:Nothing, P, Si}; max_records = Inf, verbose = true) where {Si <: AbstractSink, P <: AbstractProcessor, I}

    filepath = _filename(p.source)
    filename, extension = splitext(basename(filepath))

    reader, RecordType = open_reader(p.source, filepath, filename, extension)
    
    if is_paired(p.sink)
        filename2 = p.sink.second_in_pair(filename)
        filepath2 = joinpath(dirname(filepath), filename2 * extension)
        writer = open_writer_paired(p.sink, filepath, filename, filepath2, filename2, extension)
    else
        writer = open_writer(p.sink, filepath, filename, extension)
    end

    record = RecordType.Record()
    k = 0 
    interval_found = false

    while !eof(reader)

        read!(reader, record)

        # if record are outside interval, skip 
        !interval_found && !is_in_interval(p.source, record) && continue
        interval_found = true

        # if we are after interval, stop
        is_after_interval(p.source, record) && break
    
        out_record = p.processor(record)
        !isnothing(out_record) && write(writer, out_record)

        k += 1
        k > max_records && break
        log_progress(k, verbose, filepath)
    end
    close(reader)
    close(writer)

    return_value(p.sink)
end

# run single file with grouper
function run_single(p::Pipeline{<:Reader{File{I}}, <:RecordGrouper, P, Si}; max_records = Inf, verbose = true) where {Si <: AbstractSink, P <: AbstractProcessor, I}

    filepath = _filename(p.source)
    filename, extension = splitext(basename(filepath))

    reader, RecordType = open_reader(p.source, filepath, filename, extension)

    if is_paired(p.sink)
        filename2 = p.sink.second_in_pair(filename)
        filepath2 = joinpath(dirname(filepath), filename2 * extension)
        
        writer = open_writer_paired(p.sink, filepath, filename, filepath2, filename2, extension)
    else
        writer = open_writer(p.sink, filepath, filename, extension)
    end
    record = RecordType.Record()
    k = 0 
    interval_found = false
    while !eof(reader)

        record, idx = get_record(p.grouper)
        read!(reader, record)

        # if record are outside interval, skip 
        if !interval_found && !is_in_interval(p.source, record) 
            free_idx!(p.grouper, idx)
            continue
        end
        interval_found = true

        # if we are after interval, stop
        is_after_interval(p.source, record) && break
        
        key, isdone = group_record!(p.grouper, record, idx)

        if isdone
            records = (p.grouper.records[i] for i in p.grouper.groups[key])
            out_records = p.processor(records...)
            
            !isnothing(out_records) && write(writer, out_records)

            for i in p.grouper.groups[key]
                free_idx!(p.grouper, i)
            end
            free_group!(p.grouper, key)
        end
        
        k += 1
        k > max_records && break
        log_progress(k, verbose, filepath)
    end
    close(reader)
    close(writer)

    return_value(p.sink)
end

# pipeline for a Bio reader, paired, with RecordGrouper
function run_paired(p::Pipeline{<:Reader{File{I}}, <:RecordGrouper, P, Si}, second_in_pair; max_records = Inf, verbose = true) where {Si <: AbstractSink, P <: AbstractProcessor, I}

    filepath1 = _filename(p.source)
    filename1, extension = splitext(basename(filepath1))

    filename2 = second_in_pair(filename1)
    filepath2 = joinpath(dirname(filepath1), filename2 * extension)

    reader1, RecordType = open_reader(p.source, filepath1, filename1, extension)
    reader2, _          = open_reader(p.source, filepath2, filename2, extension)

    writer = open_writer_paired(p.sink, filepath1, filename1, filepath2, filename2, extension)
    
    record1, record2 = RecordType.Record(), RecordType.Record()
    k = 0 
    interval_found = false
    while !eof(reader1)

        (record1, record2), idx = get_record(p.grouper)

        read!(reader1, record1)
        read!(reader2, record2)

        # if records are outside interval, skip them
        !interval_found && !is_in_interval(p.source, record1) && !is_in_interval(p.source, record2) && continue
        interval_found = true

        # if we are after interval, stop
        is_after_interval(p.source, record1) && is_after_interval(p.source, record2) && break

        key, isdone = group_record!(p.grouper, (record1, record2), idx)

        if isdone
            records = [p.grouper.records[i] for i in p.grouper.groups[key]]
            out_records = p.processor(records)
            
            if !isnothing(out_records)
                for out_record in out_records
                    !isnothing(out_record) && write(writer, out_record)
                end
            end

            for i in p.grouper.groups[key]
                free_idx!(p.grouper, i)
            end
            free_group!(p.grouper, key) # so we can reuse the same key
            
        end

        k += 1
        k > max_records && break
        log_progress(k, verbose, filepath1) 
    end
    close(reader1)
    close(reader2)
    close(writer)

    return_value(p.sink)
end

function run_paired(p::Pipeline{<:Reader{File{I}}, <:OptionalGrouper, P, Si}, second_in_pair; max_records = Inf, verbose = true) where {Si <: AbstractSink, P <: AbstractProcessor, I}

    filepath1 = _filename(p.source)
    filename1, extension = splitext(basename(filepath1))

    filename2 = second_in_pair(filename1)
    filepath2 = joinpath(dirname(filepath1), filename2 * extension)

    reader1, RecordType = open_reader(p.source, filepath1, filename1, extension)
    reader2, _          = open_reader(p.source, filepath2, filename2, extension)

    writer = open_writer_paired(p.sink, filepath1, filename1, filepath2, filename2, extension)
    
    record1, record2 = RecordType.Record(), RecordType.Record()
    k = 0 
    interval_found = false
    while !eof(reader1)

        read!(reader1, record1)
        read!(reader2, record2)

        # if records are outside interval, skip them
        !interval_found && !is_in_interval(p.source, record1) && !is_in_interval(p.source, record2) && continue
        interval_found = true

        # if we are after interval, stop
        is_after_interval(p.source, record1) && is_after_interval(p.source, record2) && break

        out_record = p.processor(record1, record2)
        !isnothing(out_record) && write(writer, out_record)

        k += 1
        k > max_records && break
        log_progress(k, verbose, filepath1)
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

    if is_paired(p.sink)
        filename2 = p.sink.second_in_pair(filename)
        filepath2 = joinpath(dirname(filepath), filename2 * extension)
        writer = open_writer_paired(p.sink, filepath, filename, filepath2, filename2, extension)
    else
        writer = open_writer(p.sink, filepath, filename, extension)
    end
    
    k = 0 
    for record in p.source.data

        out_record = p.processor(record)
        !isnothing(out_record) && write(writer, out_record)

        k += 1
        k > max_records && break
        log_progress(k, verbose, filepath)
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
            filename; second_in_pair = is_paired(p.source) ? p.source.file_provider.second_in_pair : nothing,
            interval = interval(p.source.file_provider)
            ),
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

function run(p::Pipeline{<:File{I}, <:OptionalGrouper, <:ExternalTool, Si}; max_records = Inf, verbose = true) where {Si <: OptionalSink, I}
    if is_paired(p.source)
        run_paired(p, p.source.second_in_pair; max_records=max_records, verbose=verbose)
    else
        run_single(p; max_records=max_records, verbose=verbose)
    end
end

function run_single(p::Pipeline{<:File{I}, <:OptionalGrouper, <:ExternalTool, Si}; max_records = Inf, verbose = true) where {Si <: OptionalSink, I}

    filepath = _filename(p.source)
    p.processor(filepath)
end

# pipeline for a Bio reader, paired
function run_paired(p::Pipeline{<:File{I}, <:OptionalGrouper, <:ExternalTool, Si}, second_in_pair; max_records = Inf, verbose = true) where {Si <: OptionalSink, I}

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
            filename; second_in_pair = is_paired(p.source) ? p.source.second_in_pair : nothing,
            interval = interval(p.source)
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
