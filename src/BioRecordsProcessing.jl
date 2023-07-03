module BioRecordsProcessing

    using Glob, CodecZlib, BGZFStreams

    import Base: close, run, write

    export Pipeline, Reader, Writer, File, Directory, Collect, Buffer, ExternalTool

    include("process.jl")
    include("Source.jl")
    include("Sink.jl")
    include("Processor.jl")
    include("Pipeline.jl")

end
