
# using BioRecordsProcessing, FASTX, BioSequences
# filename = "/Users/jbieler/.julia/dev/BioRecordsProcessing/test/data/illumina_full_range_as_illumina.fastq"

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
            @show p
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
            @show p
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
            @show p
            out = run(p)
            @test all(isfile.(out))

            read_file = file -> Pipeline(Reader(FASTX.FASTQ, file), Collect(FASTX.FASTQ.Record))
            @test run(read_file(directory.files[1])) == run(read_file(out[1]))
            @test run(read_file(directory.files[1])) != run(read_file(out[2]))
        end
    end

    
    @testset "Directory + Collect" begin
        input_directory = path_of_format("FASTQ")
        directory = Directory(input_directory, "solexa*.fastq")
        p = Pipeline(
            Reader(FASTX.FASTQ, directory),
            r -> FASTQ.identifier(r),
            Collect(String)
        )
        @show p
        out = run(p)
        @test out[1] != out[2] # make sure Collect is copied
    end

    @testset "VCF" begin
        mktempdir() do dir
            filepath = joinpath(path_of_format("VCF"), "adeno_virus.vcf")
        
            p = Pipeline(
                Reader(VCF, File(filepath)),
                identity,
                Collect(VCF.Record),
            )
            @show p
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
            @show p
            @test_throws ErrorException run(p)    
        end
    end

    @testset "Buffer + Collect" begin
        
        input = rand(10)
        p = Pipeline(
            Buffer(input),
            x -> 2x,
            Collect(Float64),
        )
        @show p
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
            # make sure we don't have alias records
            @test length(unique(out[1])) == length(out[1])
            # the two files are the same
            @test sequence(out[1][2][1]) == "MPPPETPSEGRQPSPSPSPTERAPASEEEFQFLRCQQCQAEAKCPKLLPCLHTLCSGCLEASGMQCPICQ"
            @test sequence(out[1][2][2]) == "MPPPETPSEGRQPSPSPSPTERAPASEEEFQFLRCQQCQAEAKCPKLLPCLHTLCSGCLEASGMQCPICQ"
            
            @test tot_length1 == 6*70
            @test tot_length2 == 6*70
        end
    end

    @testset "Paired FASTA + Paired Writer" begin
        mktempdir() do dir

            indir = path_of_format("FASTA")
            second_in_pair = f1 -> replace(f1, "_1" => "_2")
            
            p = Pipeline(
                Reader(FASTX.FASTA, File(joinpath(indir, "multi_1.fasta"); second_in_pair = second_in_pair)),
                (r1, r2) -> begin
                    return r1, r2
                end,
                Writer(FASTX.FASTA, dir; paired = true)
            )
            out = run(p)
            
            read_file = file -> Pipeline(Reader(FASTX.FASTA, file), Collect(FASTX.FASTA.Record)) |> run

            @test read_file(out[1]) == read_file(joinpath(indir, "multi_1.fasta"))
        end
    end

    @testset "Paired FASTA + Single Writer" begin
        mktempdir() do dir

            indir = path_of_format("FASTA")
            second_in_pair = f1 -> replace(f1, "_1" => "_2")
            
            # read paired and write single file
            p = Pipeline(
                Reader(FASTX.FASTA, File(joinpath(indir, "multi_1.fasta"); second_in_pair = second_in_pair)),
                (r1, r2) -> begin
                    return r1, r2
                end,
                Writer(FASTX.FASTA, dir; paired = false)
            )
            single_file = run(p)
            
            # read single and write paired
            single_2_paired(single_file) = begin
                odd_read = [FASTX.FASTA.Record()]
                k = 0
                p = Pipeline(
                    Reader(FASTX.FASTA, single_file),
                    r -> begin
                        k += 1
                        
                        if isodd(k)
                            odd_read[1] = copy(r)
                            return nothing
                        else
                            return (odd_read[1], r)
                        end
                    end,
                    Writer(FASTX.FASTA, dir; paired = true, second_in_pair = second_in_pair, suffix = ".test")
                )
                out = run(p)
            end
            out = single_2_paired(single_file)
            
            # compare original and paired files
            read_file = file -> Pipeline(Reader(FASTX.FASTA, file), Collect(FASTX.FASTA.Record)) |> run

            @test read_file(out[1]) == read_file(joinpath(indir, "multi_1.fasta"))
            @test read_file(out[2]) == read_file(joinpath(indir, "multi_2.fasta"))
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
            @test length(unique(out)) == length(out)
            @test XAM.BAM.sequence(out[end]) == dna"GGACTTGGCGGTACTTTATATCCATCTAGAGGAGCCTGTTCTATAATCGATAAACCCCGCTCTACCTCACC"
        end
    end

    # Note : BWA can write additional alignements, in that case you can get more than two
    # reads with the same name
    @testset "BAM + Collect + Group by read name" begin
        mktempdir() do dir
            spec = list_valid_specimens("BAM")
            bam = joinpath(path_of_format("BAM"), "bam1.bam")

            p = Pipeline(
                Reader(BAM, File(joinpath(dir, bam))),
                BAMPairedReadGrouper(),
                (r1,r2) -> BAM.tempname(r1) == BAM.tempname(r2),
                Collect(Bool)
            )
            @show p
            out = run(p)
            @test length(out) == 100 # 200 reads in the bam
            @test all(out)

        end
    end

    @testset "BAM to paired FASTQ" begin
        mktempdir() do dir
            spec = list_valid_specimens("BAM")
            bam = joinpath(path_of_format("BAM"), "bam1.bam")

            p = Pipeline(
                Reader(XAM.BAM, File(joinpath(dir, bam))),
                BAMPairedReadGrouper(),
                (r1,r2) ->  begin
                    @assert BAM.tempname(r1) == BAM.tempname(r2)
                    return (FASTA.Record(BAM.tempname(r1), BAM.sequence(r1)),
                            FASTA.Record(BAM.tempname(r2), BAM.sequence(r2)))
                end,
                Writer(
                    FASTX.FASTA, dir; 
                    paired = true, 
                    second_in_pair = x -> "bam1_2",
                    extension = ".fasta.gz"
                )
            )
            out = run(p)
            read_file = file -> run(Pipeline(Reader(FASTX.FASTA, file), Collect(FASTX.FASTA.Record)))
            
            r1 = read_file(out[1])
            r2 = read_file(out[2])
            
            @test length(r1) == 100
            @test length(r2) == 100

            @test FASTA.identifier(r1[1]) == "HWI-1KL120:88:D0LRBACXX:1:1101:2852:2134"
            @test FASTA.identifier(r2[1]) == "HWI-1KL120:88:D0LRBACXX:1:1101:2852:2134"
            
            @test FASTA.sequence(r1[1]) == "GAGAGGTCAGCGTGAGCCCCTTGCCTCACACCGGCCCCTCTCACGCCGAGAGAGGTCAGCGTGAGCCCCTTGCCTCACACCGGCCCCTCCCACGCCGAGAG"
            @test FASTA.sequence(r2[1]) == "TCACGGTGGCCTGTTGAGGCAGGGGCTCACGCTGACCTCTCTCGGCGTGGGAGGGGCCGGTGTGAGGCAAGGGCTCACGCTGACCTCTCTCGGCGTGGGAG"
            
        end
    end

end
