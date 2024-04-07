module BcubeTutorialsTests
using BcubeTutorials
using Test
using SHA: sha1
using DelimitedFiles
using Serialization

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

function get_ref_checksum(key)
    # Reading sha1 checksums
    f = readdlm(joinpath(@__DIR__, "checksums.sha1"), String)
    fname2sum = Dict(r[2] => r[1] for r in eachrow(f))
    return fname2sum[key]
end

function compute_checksum(value)
    path, io = mktemp()
    serialize(path, value)
    return bytes2hex(open(sha1, path))
end

"""
    check_value(value, key)

Compare the checksum of `value` with the reference checksum
store for `key` in file "checksums.sha1"
"""
function check_value(value, key)
    compute_checksum(value) == get_ref_checksum(key)
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
