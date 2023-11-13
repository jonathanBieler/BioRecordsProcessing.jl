mutable struct RecordGrouper{R, G, F <: Function, FD <: Function}
    records::Vector{R}
    free_idx::Set{Int}
    groups::Dict{G, Vector{Int}}
    get_key::F
    isdone::FD
    
    RecordGrouper{R,G}(get_key::Function, isdone::Function) where {R, G} = 
        new{R, G, typeof(get_key), typeof(isdone)}(
            [R() for i=1:10_000],
            Set{Int}(1:10_000),
            Dict{G, Vector{Int}}(),
            get_key,
            isdone
        )
end

function check_free_idx!(rg::RecordGrouper{R,G}) where {R,G}
    if isempty(rg.free_idx)
        N = length(rg.records)
        @warn "No space left in RecordGrouper, increasing buffer size to $(2N) Records"
        for i = N:2N
            push!(rg.records, R())
            push!(rg.free_idx, i)
        end 
    end
end

function free_idx!(rg::RecordGrouper{R,G}, idx) where {R,G}
    push!(rg.free_idx, idx)
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

BAMPairedReadGrouper() = RecordGrouper{BAM.Record, String}(
    r -> BAM.tempname(r),
    records -> length(records) >= 2
)