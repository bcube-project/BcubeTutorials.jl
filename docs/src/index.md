```@meta
CurrentModule = BcubeTutorials
```

# BcubeTutorials

This website gather several commented tutorials as well as various examples of application of the library [Bcube](https://github.com/bcube-project/Bcube.jl).

Bcube is a Julia library providing tools for the spatial discretization of partial differential equation(s) (PDE). The main objective is to provide a set of tools to quickly assemble an algorithm solving partial differential equation(s) efficiently.

## Tutorials vs examples

Tutorials will always follow the `Bcube` development : they will be kept updated. Examples are more advanced use-cases of `Bcube` and could be compatible only with a specific version of `Bcube` (and some other dependencies such as `DifferentialEquations.jl`).

## Run the scripts locally

All the tutorials (and most of the examples) can be ran locally with the following steps.

First, clone the repository

```bash
$ git clone https://bcube-project.github.io/BcubeTutorials.jl
$ cd BcubeTutorials.jl/
```

Then, set up the environnement

```julia-repl
julia> using Pkg
julia> Pkg.activate(".")
julia> Pkg.add(PackageSpec(url="https://github.com/bcube-project/Bcube.jl"))
julia> Pkg.instantiate()
julia> include("src/tutorial/helmholtz.jl")
```

## Authors

Ghislain Blanchard, Lokman Bennani and Maxime Bouyges.
