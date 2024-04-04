module BioRecordsProcessing

    using Glob, CodecZlib, BGZFStreams, XAM, GenomicFeatures

    import Base: close, run, write

    export Pipeline, Reader, Writer, File, Directory, Collect, Buffer, ExternalTool
    export RecordGrouper, BAMPairedReadGrouper

    include("Source.jl")
    include("Sink.jl")
    include("Processor.jl")
    include("RecordGrouper.jl")
    include("Pipeline.jl")
    
end
