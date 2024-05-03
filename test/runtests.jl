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

include("./utils.jl")

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
    custom_include("../src/tutorial/heat_equation.jl")
    custom_include("../src/tutorial/helmholtz.jl")
    custom_include("../src/tutorial/linear_transport.jl")

    custom_include("../src/example/constrained_poisson/constrained_poisson.jl")
end

end
