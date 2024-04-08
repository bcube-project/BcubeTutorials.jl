module BcubeTutorialsTests
using BcubeTutorials
using Test
using SHA: sha256
using DelimitedFiles
using Printf: Format, format
using FileIO, JLD2
using ReferenceTests
using SparseArrays

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
    return get(fname2sum, key, nothing)
end

scientific_format(x, digits = 10) = format(Format("%10.$(digits)f"), x)

function compute_checksum(value::AbstractArray{<:Number}; digits = 10)
    path, _ = mktemp()
    str_value = map(Base.Fix2(scientific_format, digits), round.(value, digits = digits))
    writedlm(path, str_value)
    checksum = bytes2hex(open(sha256, path))
    return checksum
end

"""
    compare_checksum(key, value; digits = 10)

Compare the checksum of `value` with the reference checksum
store for `key` in file "checksums.sha256".
Keyword argument `digits` can be use to control rounding of `value`
before computing the checksum. This option avoids test failure
due to approximative floating-point precision.
Note that the reference checksum must be generated with the same number of `digits`.
"""
function compare_checksum(key, value; digits = 10, verbose = true)
    _ref = get_ref_checksum(key)
    _cur = compute_checksum(value; digits = digits)
    result = _cur == _ref
    if (!result) && verbose
        printstyled("   ➡ TEST FAIL : key = ", key, " \n"; color = :light_black)
        printstyled("      ↪ checksum ref = ", _ref, " \n"; color = :light_black)
        printstyled("      ↪ checksum     = ", _cur, " \n"; color = :light_black)
    end
    return result
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
    # Reminder from the official doc for `Dict`:
    # "When the values are stored internally in a hash table,
    # as is the case for Dict, the order in which they are returned may vary"
    length(keys(d1)) != length(keys(d2)) && return false
    any(x -> !haskey(d1, x), keys(d2)) && return false
    for key in keys(d1)
        any((!isapprox).(d1[key], d2[key]; atol = atol, rtol = rtol)) && return false
    end
    return true
end
comp(atol::Real, rtol::Real) = (a, b) -> comp(a, b, atol, rtol)
comp() = comp(1.0e-12, 1.0e-12)

refpath(filename) = joinpath(@__DIR__, "references/", filename)

function test_ref(filename_ref::String, dict::Dict, comp::Function = comp())
    @test_reference refpath(filename_ref) dict by = comp
end
function test_ref(filename_ref::String, data, args...; kwargs...)
    test_ref(filename_ref, _as_dict(data), args...; kwargs...)
end
_as_dict(a::AbstractArray) = Dict("array" => a)
_as_dict(a::AbstractSparseMatrix) = Dict(zip("sparse_" .* ("I", "J", "V"), findnz(a)))

ENV["TestMode"] = "true"
@testset "BcubeTutorials.jl" begin
    @testset "compare_checksum" begin
        a = [1.0, 2.0, 3.0]
        b = a .+ 2.0e-11
        c = a .+ 2.0e-10
        @test compare_checksum("checkvalue_vector", a; digits = 10)
        @test compare_checksum("checkvalue_vector", b; digits = 10)
        @test !compare_checksum("checkvalue_vector", c; digits = 10, verbose = false)
    end
    custom_include("../src/tutorial/helmholtz.jl")
end

end
