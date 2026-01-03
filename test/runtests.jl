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
    @testset "Tutorials" begin
        custom_include("../src/tutorial/heat_equation.jl")
        custom_include("../src/tutorial/helmholtz.jl")
        custom_include("../src/tutorial/linear_transport.jl")
    end
    @testset "Examples" begin
        custom_include("../src/example/heat_equation_sphere/heat_equation_sphere.jl")
        custom_include("../src/example/constrained_poisson/constrained_poisson.jl")
        custom_include(
            "../src/example/heat_equation_two_layers/heat_equation_two_layers.jl",
        )
        custom_include("../src/example/covo/covo.jl")
    end
end

end
