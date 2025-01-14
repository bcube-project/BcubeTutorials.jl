module HeatEquationSphere #hide
# # Heat equation on a sphere
# Adapted from https://www.chebfun.org/examples/sphere/SphereHeatConduction.html.
# The equation solved on the sphere is
# ```math
#   \partial_t T = \alpha \Delta_\Gamma T
# ```
# The temperature is initialized by
# the following sum of spherical harmonic:
# ```math
#   u(\lambda, \theta) = Y^0_6(\lambda, \theta) + \sqrt{14/11}Y^6_6(\lambda,\theta)
# ```
# where $\lambda$ and $\theta$ are two angles parametrizing the sphere.
# Hence an analytical solution can be found : $u(\lambda,\theta, t) = e^{-42 \alpha t} u_0(\lambda, \theta)$.
#
# The following animation is obtained after running the simulation:
#
# ![](../assets/heat_equation_sphere.gif)
#
using Bcube
using BcubeVTK
using BcubeGmsh
using LinearAlgebra
using StaticArrays
using FastTransforms
using Random
using ProgressMeter
using Test #src

"""
From wikipedia (physics)

Other mathematical valid possibilites:
* CoordinateTransformations: https://github.com/JuliaGeometry/CoordinateTransformations.jl/blob/b97c74a35c6835f6b015440cb4af6298824dbe1d/src/coordinatesystems.jl#L161
* Custom : θ = atan(z, sqrt(x^2 + y^2)), ϕ = atan(y, x)
"""
function cart2sphere(xyz)
    x, y, z = xyz

    r = norm(xyz)

    θ = acos(z / r)
    ϕ = sign(y) * acos(x / sqrt(x^2 + y^2))

    return SA[r, θ, ϕ]
end

function rotMat(θx, θy, θz)
    Rx = @SMatrix[
        1.0 0.0 0.0
        0.0 cos(θx) sin(θx)
        0.0 (-sin(θx)) cos(θx)
    ]
    Ry = @SMatrix[
        cos(θy) 0.0 (-sin(θy))
        0.0 1.0 0.0
        sin(θy) 0.0 cos(θy)
    ]
    Rz = @SMatrix[
        cos(θz) sin(θz) 0.0
        -sin(θz) cos(θz) 0.0
        0.0 0.0 1.0
    ]
    return Rx * Ry * Rz
end

mutable struct VtkHandler
    basename::Any
    ite::Any
    mesh::Any
    θ::Any
    θ_centers::Any
    ϕ::Any
    ϕ_centers::Any
    ν::Any
    ν_centers::Any
    function VtkHandler(basename, mesh)
        θ = PhysicalFunction(x -> cart2sphere(x)[2])
        θ_centers = MeshCellData(var_on_centers(θ, mesh))

        ϕ = PhysicalFunction(x -> cart2sphere(x)[3])
        ϕ_centers = MeshCellData(var_on_centers(ϕ, mesh))

        ν = get_cell_normals(CellDomain(mesh))
        ν_centers = var_on_centers(ν, mesh)
        _ν_centers = MeshCellData([SA[ν_centers[i, :]...] for i in 1:ncells(mesh)])

        new(basename, 0, mesh, θ, θ_centers, ϕ, ϕ_centers, ν, _ν_centers)
    end
end

function append_vtk(vtk, u::Bcube.AbstractFEFunction, t)
    values_centers = MeshCellData(var_on_centers(u, vtk.mesh))

    ## Write
    write_file(
        vtk.basename,
        vtk.mesh,
        Dict(
            "u_centers" => values_centers,
            "u_vertices" => u,
            "θ_centers" => vtk.θ_centers,
            "θ_vertices" => vtk.θ,
            "ϕ_centers" => vtk.ϕ_centers,
            "ϕ_vertices" => vtk.ϕ,
            "ν_centers" => vtk.ν_centers,
            "ν_vertices" => vtk.ν,
        ),
        vtk.ite,
        t;
        collection_append = vtk.ite > 0,
    )

    ## Update counter
    vtk.ite += 1
end

function run(;
    α = 1.0 / 42,
    degree = 1,
    tfinal = 1,
    nout = 100,
    lc = 1e-1,
    CFL = 0.5,
    vtk_output = true,
)

    ## Settings
    out_dir = joinpath(@__DIR__, "..", "..", "..", "myout", "heat_eqn_sphere")
    mkpath(out_dir)

    ## Mesh
    mesh_path = joinpath(out_dir, "mesh.msh")
    BcubeGmsh.gen_sphere_mesh(mesh_path; radius = 1.0, lc = lc)
    mesh = read_mesh(mesh_path)
    rng = MersenneTwister(0)
    R = rotMat(rand(rng, 3)...)
    transform!(mesh, x -> R * x) # rotate to avoid being "aligned" with an axis

    ## Domain
    dΩ = Measure(CellDomain(mesh), 2 * degree + 1)

    ## Prepare output
    vtk = VtkHandler(joinpath(out_dir, "result_d$(degree).pvd"), mesh)

    ## Discretization
    U = TrialFESpace(FunctionSpace(:Lagrange, degree), mesh)
    V = TestFESpace(U)
    u = FEFunction(U)

    ## Display infos
    println("degree = $degree, CFL = $CFL, lc = $lc, ndofs = $(get_ndofs(U))")

    ## Init
    u0 = PhysicalFunction(
        x -> begin
            _r, _θ, _ϕ = cart2sphere(x)
            return sphevaluate(_θ, _ϕ, 6, 0) + √(14 / 11) * sphevaluate(_θ, _ϕ, 6, 5)
        end,
    )
    projection_l2!(u, u0, mesh)

    ## Bilinear forms
    m(u, v) = ∫(u ⋅ v)dΩ
    a(u, v) = ∫(∇ₛ(u) ⋅ ∇ₛ(v))dΩ

    ## Assemble
    M = assemble_bilinear(m, U, V)
    K = assemble_bilinear(a, U, V)

    ## Time loop
    Δt = CFL * lc^2 / α
    nite = floor(Int, tfinal / Δt)
    _nout = min(nite, nout)
    t = 0.0
    b = Bcube.allocate_dofs(U)
    vtk_output && append_vtk(vtk, u, t)
    progress = Progress(nite)
    for ite in 1:nite
        b .= (M + Δt * α * K) \ (M * get_dof_values(u))
        set_dof_values!(u, b)

        t += Δt
        next!(progress)

        ## Output results
        if ite % (nite ÷ _nout) == 0
            vtk_output && append_vtk(vtk, u, t)
        end
    end

    ## Compute L2 error with the analytical solution
    utrue = exp(-42 * α * t) * u0
    errL2 = sum(Bcube.compute(∫((utrue - u)^2)dΩ))
    return get_ndofs(U), errL2
end

ndofs, errL2 = run(;
    degree = 1,
    α = 1.0 / 42,
    tfinal = 1,
    CFL = 0.5,
    nout = 100,
    lc = 4e-2,
    vtk_output = true,
)
@show errL2

if get(ENV, "TestMode", "false") == "true" #src
    @test errL2 < 4.83e-5                  #src
end                                        #src

end #hide
