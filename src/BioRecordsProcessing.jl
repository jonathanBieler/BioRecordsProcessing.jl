module BioRecordsProcessing

    using Glob, CodecZlib, BGZFStreams

    import Base: close, run, write

    export Pipeline, Reader, Writer, File, Collect

    include("process.jl")
    include("Pipeline.jl")

end
