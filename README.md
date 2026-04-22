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

## Developer guide

### Running the tests

The test suite can be executed locally with:

```julia-repl
pkg> activate .
(BcubeTutorials) pkg> test
```

### Testing with a custom Bcube branch

By default, the tests use the version of Bcube.jl specified in the `Project.toml` files. However, you may want to run the tests against a specific branch of Bcube.jl (e.g., for development or debugging purposes).

**Locally:** Set the `BCUBE_BRANCH` environment variable before running the tests:

```bash
export BCUBE_BRANCH=dev
julia --project=test test/runtests.jl
```

**Via GitHub Actions (manual trigger):** Go to the Actions tab, select the "CI" workflow, click "Run workflow", and enter the branch name in the "Bcube branch to use for tests" field.

**Via PR:** On a pull request, include the following syntax in the PR description or in the latest commit message to trigger CI with a custom Bcube branch:

```
/bcube-branch <branch_name>
```

For example: `/bcube-branch dev`

The same mechanism works for `BcubeVTK` and `BcubeGmsh`, and all can be combined. For instance with the commit message:
```
/bcube-branch <name_1>
/bcubegmsh-branch <name_2>
/bcubevtk-branch <name_3>
```

## Authors

Ghislain Blanchard, Lokman Bennani and Maxime Bouyges.
