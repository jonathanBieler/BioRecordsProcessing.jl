
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

    @testset "VCF overwrite" begin
        mktempdir() do dir
            filepath = joinpath(path_of_format("VCF"), "adeno_virus.vcf")
        
            p = Pipeline(
                Reader(VCF, File(filepath)),
                Writer(VCF, dirname(filepath)),
            )
            @test_throws AssertionError run(p)    
        end
    end

    @testset "Buffer + Collect" begin
        
        input = rand(10)
        p = Pipeline(
            Buffer(input),
            x -> 2x,
            Collect(Float64),
        )
        output = run(p)
        @test all(output .≈ 2*input)
    
    end

    @testset "Buffer + Writer" begin
        mktempdir() do dir
            input = [FASTX.FASTA.Record("test$i",randdnaseq(100)) for i in 1:3]
            p = Pipeline(
                Buffer(input; filename = "test.fa"),
                Writer(FASTX.FASTA, dir),
            )
            filename = run(p)
            p = Pipeline(
                Reader(FASTX.FASTA, filename),
                Collect(FASTX.FASTA.Record),
            )
            output = run(p)
            @test output == input
        end
    end

end