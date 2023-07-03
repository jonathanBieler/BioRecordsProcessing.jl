
using BioRecordsProcessing, FASTX, BioSequences
filename = "/Users/jbieler/.julia/dev/BioRecordsProcessing/test/data/illumina_full_range_as_illumina.fastq"

@testset "ExternalTool + File" begin
    mktempdir() do dir
        spec = list_valid_specimens("SAM")
        bam = joinpath(path_of_format("SAM"), "ce#1.sam")

        p = Pipeline(
            File(bam),
            ExternalTool(filepath ->
                read(`head -1 $filepath`, String)
            ),
        )
        out = run(p)
        @test out == "@SQ\tSN:CHROMOSOME_I\tLN:1009800\n"
    end
end

@testset "ExternalTool + Paired File" begin
    
    p = Pipeline(
        File("test1", second_in_pair = x -> "test2"),
        ExternalTool(
            (filepath1, filepath2) ->
            filepath1 * filepath2
        ),
    )
    out = run(p)
    @test out == "test1test2"
end

@testset "ExternalTool + Directory" begin
    mktempdir() do dir
        samdir = path_of_format("SAM")
        
        p = Pipeline(
            Directory(samdir, "xx#*"),
            ExternalTool(filepath ->
                read(`head -1 $filepath`, String)
            ),
        )
        out = run(p)
        @test out == ["", "@SQ\tSN:xx\tLN:20\n"]
    end
end

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
    @testset "FASTA + Directory" begin
        mktempdir() do dir
                
            input_directory = joinpath(@__DIR__, "data")
            directory = Directory(input_directory, "*.fastq")
            p = Pipeline(
                Reader(FASTX.FASTQ, directory),
                identity,
                Writer(
                    FASTX.FASTQ, dir;
                    suffix = ".processed"
                ),
            )
            out = run(p)
            @test all(isfile.(out))

            read_file = file -> Pipeline(Reader(FASTX.FASTQ, file), Collect(FASTX.FASTQ.Record))
            @test run(read_file(directory.files[1])) == run(read_file(out[1]))
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
            input = [FASTX.FASTA.Record("test$i", randdnaseq(100)) for i in 1:3]
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

    @testset "Paired FASTA + Collect" begin
        mktempdir() do dir

            indir = path_of_format("FASTA")
            second_in_pair = f1 -> replace(f1, "_1" => "_2")
            
            tot_length1 = tot_length2 = 0
            
            p = Pipeline(
                Reader(FASTX.FASTA, Directory(indir, "multi_1.fasta"; second_in_pair = second_in_pair)),
                (r1, r2) -> begin
                    tot_length1 += length(sequence(r1))
                    tot_length2 += length(sequence(r2))
                    return r1, r2
                end,
                Collect(FASTX.FASTA.Record, paired=true)
            )
            out = run(p)
            # the two files are the same
            @test sequence(out[1][2][1]) == "MPPPETPSEGRQPSPSPSPTERAPASEEEFQFLRCQQCQAEAKCPKLLPCLHTLCSGCLEASGMQCPICQ"
            @test sequence(out[1][2][2]) == "MPPPETPSEGRQPSPSPSPTERAPASEEEFQFLRCQQCQAEAKCPKLLPCLHTLCSGCLEASGMQCPICQ"
            
            @test tot_length1 == 6*70
            @test tot_length2 == 6*70
        end
    end

    @testset "Paired FASTA + Writer" begin
        mktempdir() do dir

            indir = path_of_format("FASTA")
            second_in_pair = f1 -> replace(f1, "_1" => "_2")
            
            p = Pipeline(
                Reader(FASTX.FASTA, File(joinpath(indir, "multi_1.fasta"); second_in_pair = second_in_pair)),
                (r1, r2) -> begin
                    return r1, r2
                end,
                Writer(FASTX.FASTA, dir)
            )
            out = run(p)
            read_file = file -> Pipeline(Reader(FASTX.FASTA, file), Collect(FASTX.FASTA.Record)) |> run

            @test read_file(out) == read_file(joinpath(indir, "multi_1.fasta"))
        end
    end

    @testset "BAM + Collect" begin
        mktempdir() do dir
            spec = list_valid_specimens("BAM")
            bam = joinpath(path_of_format("BAM"), "R_12h_D06.uniq.q40.bam")

            p = Pipeline(
                Reader(XAM.BAM, File(joinpath(dir, bam))),
                Collect(XAM.BAM.Record)
            )
            out = run(p)
            @test XAM.BAM.sequence(out[end]) == dna"GGACTTGGCGGTACTTTATATCCATCTAGAGGAGCCTGTTCTATAATCGATAAACCCCGCTCTACCTCACC"
        end
    end

end
