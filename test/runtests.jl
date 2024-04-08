module BcubeTutorialsTests
using BcubeTutorials
using Test
using SHA: sha256
using DelimitedFiles
using Printf: Format, format
using FileIO, JLD2
using ReferenceTests

# This dir will be removed at the end of the tests
tempdir = mktempdir()

# function test_files(dir, files2check)
#     # Reading sha256 checksums
#     f = readdlm(joinpath(@__DIR__, "checksums.sha256"), String)
#     fname2sum = Dict(r[2] => r[1] for r in eachrow(f))
#     for f in files2check
#         printstyled("   ➡ testing output file: ", f, "\n"; color = :light_black)
#         @test fname2sum[f] == bytes2hex(open(sha256, joinpath(dir, f)))
#     end
# end

function get_ref_checksum(key)
    # Reading sha256 checksums
    f = readdlm(joinpath(@__DIR__, "checksums.sha256"), String)
    fname2sum = Dict(r[2] => r[1] for r in eachrow(f))
    return fname2sum[key]
end

scientific_format(x, digits = 10) = format(Format("%10.$(digits)f"), x)

function compute_checksum(value::AbstractArray{<:Number}; digits = 10)
    #return bytes2hex(sha256(reinterpret(UInt8, vec(round.(value; digits = digits)))))
    path, _ = mktemp()
    str_value = map(Base.Fix2(scientific_format, digits), round.(value, digits = digits))
    writedlm(path, str_value)
    checksum = bytes2hex(open(sha256, path))
    return checksum
end

"""
    check_value(value, key; digits = 10)

Compare the checksum of `value` with the reference checksum
store for `key` in file "checksums.sha256".
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

function comp(
    d1::Dict{String, <:AbstractArray},
    d2::Dict{String, <:AbstractArray},
    rtol,
    atol,
)
    keys(d1) != keys(d2) && return false
    for (v1, v2) in zip(values(d1), values(d2))
        any((!isapprox).(v1, v2; atol = atol, rtol = rtol)) && return false
    end
    return true
end
comp(atol::Real, rtol::Real) = (a, b) -> comp(a, b, atol, rtol)
comp() = comp(1.0e-12, 1.0e-12)

refpath(filename) = joinpath(@__DIR__, "references/", filename)

function test_ref(filename_ref::String, dict::Dict, comp::Function = comp())
    @test_reference refpath(filename_ref) dict by = comp
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
