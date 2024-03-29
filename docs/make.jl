push!(LOAD_PATH, "../src/")

using BcubeTutorials
using Documenter
using Literate

# Alias for `Literate.markdown`
function gen_markdown_with_literate(src, name, dir)
    Literate.markdown(joinpath(src, name), dir; documenter = false, execute = false)
end

"""
Build a markdown file with just the content of the julia file in it.
"""
function julia_to_markdown(src_dir, target_dir, filename, title)
    filename_noext = splitext(filename)[1]
    open(joinpath(target_dir, "$(filename_noext).md"), "w") do io
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
map(
    filename -> gen_markdown_with_literate(tutorial_src, "$(filename).jl", tutorial_dir),
    tutorial_names,
)

# Generate examples
# `documenter = false` to avoid Documenter to execute cells
example_src = joinpath(@__DIR__, "..", "src", "example")
example_dir = joinpath(@__DIR__, "src", "example")
Sys.rm(example_dir; recursive = true, force = true)
mkdir(example_dir)
# gen_markdown(example_src, "euler_naca_steady.jl", example_dir)
# gen_markdown(example_src, "covo.jl", example_dir)
# gen_markdown(example_src, "linear_elasticity.jl", example_dir)

# Generate "uncommented" examples (= without `Literate`)
for (script_name, name) in (
    ("linear_elasticity.jl", "Linear elasticity"),
    ("linear_thermoelasticity.jl", "Linear thermo-elasticity"),
    ("covo.jl", "Euler equations - covo"),
    ("euler_naca_steady.jl", "Euler equations on a NACA0012"),
    ("shallow_water.jl", "Shallow water"),
    ("poisson_dg.jl", "Poisson equation (DG)"),
)
    julia_to_markdown(
        joinpath(example_src, splitext(script_name)[1]),
        example_dir,
        script_name,
        name,
    )
end

# Generator markdown with `Literate`
gen_markdown_with_literate(
    joinpath(example_src, "linear_thermoelasticity"),
    "linear_thermoelasticity.jl",
    example_dir,
)
gen_markdown_with_literate(
    joinpath(example_src, "constrained_poisson"),
    "constrained_poisson.jl",
    example_dir,
)
gen_markdown_with_literate(
    joinpath(example_src, "transport_supg"),
    "transport_supg.jl",
    example_dir,
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
            "example/constrained_poisson.md",
            "example/transport_supg.md",
            "example/poisson_dg.md",
        ],
    ],
    # remotes = nothing, # tmp fix for bmxam windows
)

deploydocs(; repo = "github.com/bcube-project/BcubeTutorials.jl.git", push_preview = true)
