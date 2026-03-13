module BcubeTutorialsTests
using Pkg
using Test
using ReferenceTests
using SparseArrays

# This dir will be removed at the end of the tests
tempdir = mktempdir()

include("./utils.jl")

ENV["TestMode"] = "true"
has_custom_bcube_branch = haskey(ENV, "BCUBE_BRANCH")
bcubeSpec = if has_custom_bcube_branch
    branch = get(ENV, "BCUBE_BRANCH")
    @info "Running tests with custom Bcube branch '$branch'"
    Pkg.PackageSpec(; name = "Bcube", rev = branch)
else
    nothing
end

SRC_DIR = joinpath(@__DIR__, "..", "src")
names = (
    "constrained_poisson",
    "covo",
    "linear_transport",
    "heat_equation",
    "heat_equation_two_layers",
    "heat_equation_sphere",
    "helmholtz",
)

@testset "BcubeTutorials.jl" begin
    for name in names
        filename = name * ".jl"
        dir = joinpath(SRC_DIR, name)
        filepath = joinpath(dir, filename)
        @testset "$filename" begin
            Pkg.activate(dir)
            has_custom_bcube_branch && Pkg.add(bcubeSpec)
            Pkg.instantiate()
            include(filepath)
        end
    end
end

end
