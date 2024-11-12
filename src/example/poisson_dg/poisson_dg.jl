module Poisson_DG #hide
println("Running poisson DG example...") #hide

# # Poisson equation (DG)
# This example is based on a [Gridap tutorial](https://gridap.github.io/Tutorials/stable/pages/t006_dg_discretization/)

# import necessary packages
using Bcube
using BcubeVTK
using LinearAlgebra
using SparseArrays
using StaticArrays

const outputpath = joinpath(@__DIR__, "..", "..", "..", "myout", "poisson_dg")
isdir(outputpath) || mkpath(outputpath)
const degree = 3
const degree_quad = 2 * degree + 1
const γ = degree * (degree + 1)
const n = 4
const Lx = 1.0
const h = Lx / n

const uₐ = PhysicalFunction(x -> 3 * x[1] + x[2]^2 + 2 * x[1]^3)
const f = PhysicalFunction(x -> -2 - 12 * x[1])
const g = uₐ

avg(u) = 0.5 * (side⁺(u) + side⁻(u))

function main()

    # Build mesh
    meshParam = (nx = n + 1, ny = n + 1, lx = Lx, ly = Lx, xc = 0.0, yc = 0.0)
    tmp_path = joinpath(tempdir(), "tmp.msh")
    gen_rectangle_mesh(tmp_path, :quad; meshParam...)
    mesh = read_msh(tmp_path)

    # Choose degree and define function space, trial space and test space
    fs = FunctionSpace(:Lagrange, degree)
    U = TrialFESpace(fs, mesh, :discontinuous)
    V = TestFESpace(U)

    # Define volume and boundary measures
    dΩ = Measure(CellDomain(mesh), degree_quad)
    dΓ = Measure(InteriorFaceDomain(mesh), degree_quad)
    dΓb = Measure(BoundaryFaceDomain(mesh), degree_quad)
    nΓ = get_face_normals(dΓ)
    nΓb = get_face_normals(dΓb)

    a_Ω(u, v) = ∫(∇(v) ⋅ ∇(u))dΩ
    l_Ω(v) = ∫(v * f)dΩ

    function a_Γ(u, v)
        ∫(
            -jump(v, nΓ) ⋅ avg(∇(u)) - avg(∇(v)) ⋅ jump(u, nΓ) +
            γ / h * jump(v, nΓ) ⋅ jump(u, nΓ),
        )dΓ
    end

    fa_Γb(u, ∇u, v, ∇v, n) = -v * (∇u ⋅ n) - (∇v ⋅ n) * u + (γ / h) * v * u
    a_Γb(u, v) = ∫(fa_Γb ∘ map(side⁻, (u, ∇(u), v, ∇(v), nΓb)))dΓb

    fl_Γb(v, ∇v, n, g) = -(∇v ⋅ n) * g + (γ / h) * v * g
    l_Γb(v) = ∫(fl_Γb ∘ map(side⁻, (v, ∇(v), nΓb, g)))dΓb

    a(u, v) = a_Ω(u, v) + a_Γ(u, v) + a_Γb(u, v)
    l(v) = l_Ω(v) + l_Γb(v)

    sys = Bcube.AffineFESystem(a, l, U, V)
    uh = Bcube.solve(sys)

    l2(u) = sqrt(sum(Bcube.compute(∫(u ⋅ u)dΩ)))
    h1(u) = sqrt(sum(Bcube.compute(∫(u ⋅ u + ∇(u) ⋅ ∇(u))dΩ)))
    e = uₐ - uh

    vars = Dict("uh" => uh, "u_ref" => uₐ, "error" => e)
    write_file(joinpath(outputpath, "output.pvd"), mesh, U, vars)

    el2 = l2(e)
    eh1 = h1(e)
    tol = 1.e-10
    @show el2, eh1
    @assert el2 < tol
    @assert eh1 < tol
end

main()

end #hide
