
using BioRecordsProcessing, FASTX, BioSequences
filename = "/Users/jbieler/.julia/dev/BioRecordsProcessing/test/data/illumina_full_range_as_illumina.fastq"

@testset "Pipeline" begin

    @testset "FASTA Collect" begin
        mktempdir() do dir
            
            filepath = joinpath(dir, "test.fa")
            writer = open(FASTA.Writer, filepath)
            for i=1:10
                write(writer, FASTA.Record("seq1", randdnaseq(i)))
            end
            close(writer)

            p = Pipeline(
                Reader(FASTX.FASTA, File(filepath)),
                record -> begin
                    length(sequence(record))
                end,
                Collect(Int),
            )
            lengths = run(p)
            @test lengths == collect(1:10)

        end
    end
    
    @testset "FASTQ.gz Writer" begin
        mktempdir() do dir
                
            filepath = joinpath(@__DIR__, "data", "illumina_full_range_as_illumina.fastq.gz")
            
            p = Pipeline(
                Reader(FASTX.FASTQ, File(filepath)),
                identity,
                Writer(
                    FASTX.FASTQ, dir;
                    suffix = ".processed"
                ),
            )
            outfile = joinpath(dir, "illumina_full_range_as_illumina.processed.fastq.gz")

            @test  outfile == run(p; max_records = 100)
            @test isfile(outfile)

        end
    end


    @testset "VCF" begin
        mktempdir() do dir
            filepath = joinpath(path_of_format("VCF"), "adeno_virus.vcf")
        
            p = Pipeline(
                Reader(VCF, File(filepath)),
                identity,
                Collect(VCF.Record),
            )
            records = run(p)
            
            @test length(unique(records)) == 14
            
        end
    end

end