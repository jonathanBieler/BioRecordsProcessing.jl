using BioRecordsProcessing
using Test, FASTX, XAM, BioSequences, FormatSpecimens
#VariantCallFormat

@testset "Internals" begin
    @test BioRecordsProcessing.insert_suffix("name", ".fastq", ".processed") == "name.processed.fastq"
    @test BioRecordsProcessing.insert_suffix("name.fastq", ".gz", ".processed") == "name.processed.fastq.gz"
end

include("test_external_tool.jl")
include("test_pipeline.jl")