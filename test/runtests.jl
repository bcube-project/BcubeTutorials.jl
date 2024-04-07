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

function get_ref_checksum(key)
    # Reading sha1 checksums
    f = readdlm(joinpath(@__DIR__, "checksums.sha1"), String)
    fname2sum = Dict(r[2] => r[1] for r in eachrow(f))
    return fname2sum[key]
end

function compute_checksum(value::AbstractArray{<:Number}; digits = 10)
    return bytes2hex(sha1(reinterpret(UInt8, vec(round.(value; digits = digits)))))
end

"""
    check_value(value, key; digits = 10)

Compare the checksum of `value` with the reference checksum
store for `key` in file "checksums.sha1".
Keyword argument `digits` can be use to control rounding of `value`
before computing the checksum. This option allows to avoid test failure
due to approximative floating-point precision.
Note that the reference checksum must be generated with the number of `digits`.
"""
function check_value(value, key; digits = 10)
    _ref = get_ref_checksum(key)
    _cur = compute_checksum(value; digits = digits)
    printstyled("   ➡ key = ", key, " \n"; color = :light_black)
    printstyled("      ↪ checksum ref = ", _ref, " \n"; color = :light_black)
    printstyled("      ↪ checksum     = ", _cur, " \n"; color = :light_black)
    return _cur == _ref
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
    @testset "check_value" begin
        a = [1.0, 2.0, 3.0]
        b = a .+ 1.0e-14 * rand(Float64, 3)
        @test check_value(a, "checkvalue_vector")
        @test check_value(b, "checkvalue_vector")
    end
    custom_include("../src/tutorial/helmholtz.jl")
end

end
