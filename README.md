# BcubeTutorials.jl

[![](https://img.shields.io/badge/docs-release-blue.svg)](https://bcube-project.github.io/BcubeTutorials.jl) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://bcube-project.github.io/BcubeTutorials.jl/dev)

Documented tutorials and various examples for the [Bcube.jl](https://github.com/bcube-project/Bcube.jl) project. Browse the [online documentation](https://bcube-project.github.io/BcubeTutorials.jl).

## Run the scripts locally

All the examples can be ran locally with the following steps.

First, clone the repository

```bash
$ git clone https://github.com/bcube-project/BcubeTutorials.jl.git
$ cd BcubeTutorials.jl/
```

Then, choose your tutorial, set up the environnement and run the script.

```julia-repl
pkg> activate ./src/helmholtz
(hemlholtz) pkg> instantiate
julia> include("src/helmholtz.jl")
```

## Build the documentation

To browse the online documentation, simply click on the blue badge at the top of this README. If you would like to build the documentation yourself, for an offline access for instance, you can do it with the following commands.

```bash
$ git clone https://github.com/bcube-project/BcubeTutorials.jl.git
$ cd BcubeTutorials.jl
```

```julia-repl
pkg> activate .
(BcubeTutorials) pkg> instantiate
julia> include("docs/make.jl")
```

You can then browse `BcubeTutorials.jl/docs/build/index.html`.

## Authors

Ghislain Blanchard, Lokman Bennani and Maxime Bouyges.
