
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