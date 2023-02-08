using BioRecordsProcessing
using Test, FASTX, XAM, VariantCallFormat, BioSequences, FormatSpecimens

@testset "BioRecordsProcessing.jl" begin
    
    @testset "FASTA" begin
        mktempdir() do dir
            
            filepath = joinpath(dir, "test.fa")
            writer = open(FASTA.Writer, filepath)
            for i=1:10
                write(writer, FASTA.Record("seq1", randdnaseq(i)))
            end
            close(writer)

            N_larger_than_5 = 0
            BioRecordsProcessing.process_directory(FASTX.FASTA, dir, "*.fa", dir; prefix = "out") do record
                N_larger_than_5 += length(sequence(record)) > 5
            
                return record
            end
            @test N_larger_than_5 == 5

        end
    end

    @testset "Paired FASTA" begin
        mktempdir() do dir

            indir = path_of_format("FASTA")
            f2 = f1 -> replace(f1, "_1" => "_2")
            
            BioRecordsProcessing.process_directory_paired(FASTX.FASTA, indir, "multi_1.fasta", f2, dir; prefix = "out") do r1, r2
                return r1, r2
            end
        end
    end

    @testset "BAM" begin
        mktempdir() do dir
            spec = list_valid_specimens("BAM")
            bam = joinpath(path_of_format("BAM"), "R_12h_D06.uniq.q40.bam")

            BioRecordsProcessing.process(XAM.BAM, bam, dir; prefix = "out") do record
                return record
            end
        end
    end

    @testset "Paired BAM" begin
        mktempdir() do dir
            rg = BioRecordsProcessing.RecordGrouper{XAM.BAM.Record, String}(
                XAM.BAM.tempname, 
                records -> length(records) == 2
            )
            spec = list_valid_specimens("BAM")
            bam = joinpath(path_of_format("BAM"), "R_12h_D06.uniq.q40.bam")

            BioRecordsProcessing.process(XAM.BAM, rg, bam, dir; prefix = "out") do r1,r2
                (r1, r2)
            end
        end
    end
    
    @testset "VCF" begin
        mktempdir() do dir
            vcf = joinpath(path_of_format("VCF"), "adeno_virus.vcf")

            BioRecordsProcessing.process(VCF, vcf, dir; prefix = "out") do record
                return VCF.filter(record) == ["PASS"] ? record : nothing 
            end
        end
    end

end


