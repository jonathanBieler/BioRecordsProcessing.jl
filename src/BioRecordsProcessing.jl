module BioRecordsProcessing

    using Glob, CodecZlib, BGZFStreams

    import Base: close, run, write

    export Pipeline, Reader, Writer, File, Collect, Buffer

    include("process.jl")
    include("Source.jl")
    include("Sink.jl")
    include("Pipeline.jl")

end
