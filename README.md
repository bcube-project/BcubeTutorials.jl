# BcubeTutorials.jl

[![](https://img.shields.io/badge/docs-release-blue.svg)](https://bcube-project.github.io/BcubeTutorials.jl)

Documented tutorials and various examples for the [Bcube.jl](https://github.com/bcube-project/Bcube.jl) project. Browse the [online documentation](https://bcube-project.github.io/BcubeTutorials.jl).

## Run the scripts locally

All the **tutorials** can be ran locally with the following steps.

First, clone the repository

```bash
$ git clone https://bcube-project.github.io/BcubeTutorials.jl
$ cd BcubeTutorials.jl/
```

Then, set up the environnement and run the script.

```julia-repl
julia> using Pkg
julia> Pkg.activate(".")
julia> Pkg.add(PackageSpec(url="https://github.com/bcube-project/Bcube.jl"))
julia> Pkg.instantiate()
julia> include("src/tutorial/helmholtz.jl")
```

Regarding the **examples**, some of them require additionnal dependencies. Hence each example is associated to a specific environment:
```julia-repl
julia> cd("src/example/covo")
julia> using Pkg
julia> Pkg.activate(".")
julia> Pkg.add(PackageSpec(url="https://github.com/bcube-project/Bcube.jl"))
julia> Pkg.instantiate()
julia> include("covo.jl")
```

## Authors

Ghislain Blanchard, Lokman Bennani and Maxime Bouyges.
