module BcubeTutorialsTests
using BcubeTutorials
using Test
using SHA: sha1
using DelimitedFiles

# This dir will be removed at the end of the tests
tempdir = mktempdir()

function test_files(dir, files2check)
    # Reading sha1 checksums
    f = readdlm(joinpath(@__DIR__, "checksums.sha1"), String)
    fname2sum = Dict(r[2] => r[1] for r in eachrow(f))
    for f in files2check
        printstyled("   ➡ testing output file: ", f, "\n"; color = :light_black)
        @test fname2sum[f] == bytes2hex(open(sha1, joinpath(dir, f)))
    end
end

"""
Custom way to "include" a file to print infos.
"""
function custom_include(path)
    filename = split(path, "/")[end]
    basename = first(splitext(filename))
    # println("Running file " * filename * "...")
    # @show basename
    @testset "$basename" begin
        include(path)
    end
    println("done.")
end

ENV["TestMode"] = "true"
@testset "BcubeTutorials.jl" begin
    custom_include("../src/tutorial/helmholtz.jl")
end

end
