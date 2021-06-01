module BioRecordsProcessing

    using Glob, CodecZlib, BGZFStreams

    function process_directory(process_function, ReadType, input_directory, pattern, output_directory;
            prefix = "",
            max_records = Inf,
            paired = false,
        )

        files = glob(pattern, input_directory)
        println.(files)
        mkpath(output_directory)

        Threads.@threads for file in files
            process(process_function, ReadType, file, output_directory;
                prefix = prefix, max_records = max_records
            )
        end
    end

    function get_writer_reader(ReadType, file, extension, out_file)
        if extension == ".gz"
            reader = ReadType.Reader(GzipDecompressorStream(open(file)))
            writer = ReadType.Writer(GzipCompressorStream(open(out_file, "w")))

        elseif extension == ".bam"
            index_file = file * ".bai"
            if !isfile(index_file)
                @warn "Index file not found : $index_file"
                reader = ReadType.Reader(open(file))
            else
                reader = ReadType.Reader(open(file); index = index_file)
            end
            
            writer = ReadType.Writer(BGZFStream(out_file, "w"), reader.header)
        else
            reader = ReadType.Reader(open(file))
            writer = ReadType.Writer(open(out_file, "w"))
        end
        reader, writer
    end

    function process(process_function, ReadType, file, output_directory;
        prefix = "",
        max_records = Inf,
    )

        filename, extension = splitext(basename(file))
        out_file = joinpath(output_directory, string(filename, prefix, extension))
        @assert out_file != file

        reader, writer = get_writer_reader(ReadType, file, extension, out_file)

        record = ReadType.Record()
        k = 0 
        while !eof(reader)

            read!(reader, record)
            out_record = process_function(record)
            !isnothing(out_record) && write(writer, out_record)

            k += 1
            k > max_records && break
            if mod(k, 100_000) == 0
                @info "$(Threads.threadid()), $(basename(file)) : Processed $(div(k, 1000))k records..."
            end     
        end
        close(reader)
        close(writer)
    end

    # paired files

    function process_directory_paired(process_function, ReadType, input_directory, pattern, to_file2, output_directory;
        prefix = "",
        max_records = Inf,
        paired = false,
    )

        files1 = glob(pattern, input_directory)
        files2 = to_file2.(files1)
        println.(zip(files1, files2))
        mkpath(output_directory)

        Threads.@threads for i in eachindex(files1)
            process_paired(process_function, ReadType, files1[i], files2[i], output_directory;
                prefix = prefix, max_records = max_records
            )
        end
    end

    function process_paired(process_function, ReadType, file1, file2, output_directory;
        prefix = "",
        max_records = Inf,
    )

        filename1, extension1 = splitext(basename(file1))
        filename2, extension2 = splitext(basename(file2))
        out_file1 = joinpath(output_directory, string(filename1, prefix, extension1))
        out_file2 = joinpath(output_directory, string(filename2, prefix, extension2))
        @assert out_file1 != file1
        @assert out_file2 != file2

        reader1, writer1 = get_writer_reader(ReadType, file1, extension1, out_file1)
        reader2, writer2 = get_writer_reader(ReadType, file2, extension2, out_file2)

        record1 = ReadType.Record()
        record2 = ReadType.Record()
        k = 0 
        while !(eof(reader1) || eof(reader2))

            read!(reader1, record1)
            read!(reader2, record2)
            out_record1, out_record2 = process_function(record1, record2)
            !isnothing(out_record1) && write(writer1, out_record1)
            !isnothing(out_record2) && write(writer2, out_record2)

            k += 2
            k > max_records && break
            if mod(k, 100_000) == 0
                @info "$(Threads.threadid()), $(basename(file1)) : Processed $(div(k, 1000))k records..."
            end     
        end
        close(reader1)
        close(reader2)
        close(writer1)
        close(writer2)
    end

    # Grouped things

    struct RecordGrouper{R, G, F <: Function, FD <: Function}
        records::Vector{R}
        free_idx::Set{Int}
        groups::Dict{G, Vector{Int}}
        get_key::F
        isdone::FD
        
        RecordGrouper{R,G}(get_key::Function, isdone::Function) where {R, G} = 
            new{R, G, typeof(get_key), typeof(isdone)}(
                [R() for i=1:1000],
                Set{Int}(1:1000),
                Dict{G, Vector{Int}}(),
                get_key,
                isdone
            )
    end

    function check_free_idx!(rg::RecordGrouper{R,G}) where {R,G}
        if isempty(rg.free_idx)
            N = length(rg.records)
            for i = N:2N
                push!(rg.records, R())
                push!(rg.free_idx, i)
            end 
        end
    end

    function get_record(rg::RecordGrouper{R,G}) where {R,G}
        check_free_idx!(rg)
        idx = pop!(rg.free_idx)
        rg.records[idx], idx
    end

    function group_record!(rg::RecordGrouper, r, idx)
        key = rg.get_key(r)
        if haskey(rg.groups, key)
            push!(rg.groups[key], idx)
        else
            rg.groups[key] = [idx]
        end
        records = (rg.records[i] for i in rg.groups[key])
        key, rg.isdone(records)
    end

    function process(process_function, ReadType, rg::RecordGrouper, file, output_directory;
        prefix = "",
        max_records = Inf,
    )

        filename, extension = splitext(basename(file))
        out_file = joinpath(output_directory, string(filename,'.', prefix, extension))
        @assert out_file != file

        reader, writer = get_writer_reader(ReadType, file, extension, out_file)

        k = 0 
        while !eof(reader)

            record, idx = get_record(rg)
            read!(reader, record)

            key, isdone = group_record!(rg, record, idx)

            if isdone
                records = (fg.records[i] for i in fg.groups[key])
                out_records = process_function(records...)
                for r in out_records
                    write(writer, r)
                end
            end

            k += 1
            k > max_records && break
            if mod(k, 100_000) == 0
                @info "$(Threads.threadid()), $(basename(f)) : Processed $(div(k, 1000))k records..."
            end     
        end
        close(reader)
        close(writer)
    end

end
