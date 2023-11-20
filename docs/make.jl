push!(LOAD_PATH, "../src/")

using BcubeTutorials
using Documenter
using Literate

# Alias for `Literate.markdown`
function gen_markdown(src, name, dir)
    Literate.markdown(joinpath(src, name), dir; documenter = false, execute = false)
end

"""
Build a markdown file with just the content of the julia file in it.
"""
function julia_to_markdown(src_dir, target_dir, filename, title)
    open(joinpath(target_dir, split(filename, ".")[1] * ".md"), "w") do io
        println(io, "# " * title)
        println(io, "```julia")
        f = open(joinpath(src_dir, filename), "r")
        lines = readlines(f)
        close(f)
        map(line -> println(io, line), lines)
        println(io, "```")
    end
end

# Generate tutorials
# `documenter = false` to avoid Documenter to execute cells
tutorial_names =
    ["helmholtz", "heat_equation", "linear_transport", "phase_field_supercooled"]
tutorial_src = joinpath(@__DIR__, "..", "src", "tutorial")
tutorial_dir = joinpath(@__DIR__, "src", "tutorial")
Sys.rm(tutorial_dir; recursive = true, force = true)
map(filename -> gen_markdown(tutorial_src, "$(filename).jl", tutorial_dir), tutorial_names)

# Generate "commented" examples
# `documenter = false` to avoid Documenter to execute cells
example_src = joinpath(@__DIR__, "..", "src", "example")
example_dir = joinpath(@__DIR__, "src", "example")
Sys.rm(example_dir; recursive = true, force = true)
mkdir(example_dir)
# gen_markdown(example_src, "euler_naca_steady.jl", example_dir)
# gen_markdown(example_src, "covo.jl", example_dir)
# gen_markdown(example_src, "linear_elasticity.jl", example_dir)

# Generate "uncommented" examples
julia_to_markdown(
    example_src,
    example_dir,
    "euler_naca_steady.jl",
    "Euler equations on a NACA0012",
)
julia_to_markdown(example_src, example_dir, "covo.jl", "Euler equations - covo")
julia_to_markdown(example_src, example_dir, "linear_elasticity.jl", "Linear elasticity")
julia_to_markdown(
    example_src,
    example_dir,
    "linear_thermoelasticity.jl",
    "Linear thermo-elasticity",
)

makedocs(;
    modules = [BcubeTutorials],
    authors = "Ghislain Blanchard, Lokman Bennani and Maxime Bouyges",
    sitename = "BcubeTutorials",
    clean = true,
    doctest = false,
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://bcube-project.github.io/BcubeTutorials.jl",
        assets = String[],
    ),
    checkdocs = :none,
    pages = [
        "Home" => "index.md",
        "Tutorials" => ["tutorial/$(filename).md" for filename in tutorial_names],
        "Advanced examples" => Any[
            "example/covo.md",
            "example/euler_naca_steady.md",
            "example/linear_elasticity.md",
            "example/linear_thermoelasticity.md",
        ],
    ],
)

deploydocs(; repo = "github.com/bcube-project/BcubeTutorials.jl.git", push_preview = true)
