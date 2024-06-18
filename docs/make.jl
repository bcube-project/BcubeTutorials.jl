push!(LOAD_PATH, "../src/")

using BcubeTutorials
using Documenter
using Literate

# Alias for `Literate.markdown`
function gen_markdown_with_literate(src, name, dir)
    Literate.markdown(
        joinpath(src, name),
        dir;
        documenter = false,
        execute = false,
        preprocess = replace_plain_code,
    )
end

# replace line "# @__PLAIN_CODE__" in `content` with the uncommented content version
function replace_plain_code(content::String)
    tempdir = mktempdir()
    filename = joinpath(tempdir, "tmp")
    write(filename, content)
    Literate.script(filename, tempdir; credit = false)
    plaincode = read(filename * ".jl", String)
    content = replace(content, "# @__PLAIN_CODE__" => plaincode)
    return content
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
    ("covo.jl", "Euler equations - covo"),
    ("linear_elasticity.jl", "Linear elasticity"),
    ("linear_thermoelasticity.jl", "Linear thermo-elasticity"),
    ("euler_naca_steady.jl", "Euler equations on a NACA0012"),
    ("shallow_water.jl", "Shallow water"),
    ("poisson_dg.jl", "Poisson equation (DG)"),
    ("heat_equation_two_layers.jl", "Heat equation with two layers"),
)
    julia_to_markdown(
        joinpath(example_src, first(splitext(script_name))),
        example_dir,
        script_name,
        name,
    )
end

# Generator markdown with `Literate`
for name in (
    "linear_thermoelasticity",
    "constrained_poisson",
    "transport_supg",
    "heat_equation_sphere",
    "transport_hypersurface",
)
    gen_markdown_with_literate(joinpath(example_src, name), "$(name).jl", example_dir)
end

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
            "example/heat_equation_sphere.md",
            "example/heat_equation_two_layers.md",
            "example/transport_hypersurface.md",
        ],
    ],
    # remotes = nothing, # tmp fix for bmxam windows
)

deploydocs(; repo = "github.com/bcube-project/BcubeTutorials.jl.git", push_preview = true)
