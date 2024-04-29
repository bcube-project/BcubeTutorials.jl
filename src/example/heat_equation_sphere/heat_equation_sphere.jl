module HeatEquationSphere
# Adapted from https://www.chebfun.org/examples/sphere/SphereHeatConduction.html
using Bcube
using WriteVTK
using LinearAlgebra
using StaticArrays
using FastTransforms
using Random
using Plots
using ProgressMeter

"""
From wikipedia

or

From CoordinateTransformations
https://github.com/JuliaGeometry/CoordinateTransformations.jl/blob/b97c74a35c6835f6b015440cb4af6298824dbe1d/src/coordinatesystems.jl#L161
"""
function cart2sphere(xyz)
    x, y, z = xyz

    r = norm(xyz)

    # From Wikipedia ("physics")
    θ = acos(z / r)
    ϕ = sign(y) * acos(x / sqrt(x^2 + y^2))

    # From "me" (adapted from wiki) -> physics
    # θ = atan(z, sqrt(x^2 + y^2))
    # ϕ = atan(y, x)

    # Formulae from CoordinateTransformations
    # θ = atan(y, x)
    # ϕ = atan(z / sqrt(x^2 + y^2))
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
    θ_vertices::Any
    ϕ::Any
    ϕ_centers::Any
    ϕ_vertices::Any
    ν::Any
    ν_centers::Any
    ν_vertices::Any
    function VtkHandler(basename, mesh)
        θ = PhysicalFunction(x -> cart2sphere(x)[2])
        θ_centers = var_on_centers(θ, mesh)
        θ_vertices = var_on_vertices(θ, mesh)

        ϕ = PhysicalFunction(x -> cart2sphere(x)[3])
        ϕ_centers = var_on_centers(ϕ, mesh)
        ϕ_vertices = var_on_vertices(ϕ, mesh)

        ν = Bcube.CellNormal(mesh)
        ν_centers = transpose(var_on_centers(ν, mesh))
        ν_vertices = transpose(var_on_vertices(ν, mesh))

        new(
            basename,
            0,
            mesh,
            θ,
            θ_centers,
            θ_vertices,
            ϕ,
            ϕ_centers,
            ϕ_vertices,
            ν,
            ν_centers,
            ν_vertices,
        )
    end
end

function append_vtk(vtk, u::Bcube.AbstractFEFunction, t)
    values_vertices = var_on_vertices(u, vtk.mesh)
    values_centers = var_on_centers(u, vtk.mesh)

    ## Write
    Bcube.write_vtk(
        vtk.basename,
        vtk.ite,
        t,
        vtk.mesh,
        Dict(
            "u_centers" => (values_centers, VTKCellData()),
            "u_vertices" => (values_vertices, VTKPointData()),
            "θ_centers" => (vtk.θ_centers, VTKCellData()),
            "θ_vertices" => (vtk.θ_vertices, VTKPointData()),
            "ϕ_centers" => (vtk.ϕ_centers, VTKCellData()),
            "ϕ_vertices" => (vtk.ϕ_vertices, VTKPointData()),
            "ν_centers" => (vtk.ν_centers, VTKCellData()),
            "ν_vertices" => (vtk.ν_vertices, VTKPointData()),
        ),
        ;
        append = vtk.ite > 0,
    )

    ## Update counter
    vtk.ite += 1
end

function run(;
    α = 1.0 / 42,
    degree = 1,
    tfinal = 1,
    Δt = 1e-2,
    nout = 100,
    lc = 1e-1,
    vtk_output = true,
)

    # Settings
    out_dir = joinpath(@__DIR__, "../../../myout/heat_eqn_sphere")
    mkpath(out_dir)

    # Mesh
    mesh_path = joinpath(out_dir, "mesh.msh")
    Bcube.gen_sphere_mesh(mesh_path; radius = 1.0, lc = lc) # 1e-1
    mesh = read_msh(mesh_path)
    rng = MersenneTwister(0)
    R = rotMat(rand(rng, 3)...)
    transform!(mesh, x -> R * x) # rotate to avoid being "aligned" with an axis

    # Domain
    dΩ = Measure(CellDomain(mesh), 2 * degree + 1)

    # Prepare output
    vtk = VtkHandler(joinpath(out_dir, "result_d$(degree)"), mesh)

    # Discretization
    U = TrialFESpace(FunctionSpace(:Lagrange, degree), mesh)
    V = TestFESpace(U)
    u = FEFunction(U)
    @show get_ndofs(U)

    # Init
    u0 = PhysicalFunction(
        x -> begin
            _r, _θ, _ϕ = cart2sphere(x)
            return sphevaluate(_θ, _ϕ, 6, 0) + √(14 / 11) * sphevaluate(_θ, _ϕ, 6, 5)
        end,
    )
    projection_l2!(u, u0, mesh)

    # Bilinear forms
    m(u, v) = ∫(u ⋅ v)dΩ
    a(u, v) = ∫(∇ₛ(u) ⋅ ∇ₛ(v))dΩ

    # Assemble
    M = assemble_bilinear(m, U, V)
    K = assemble_bilinear(a, U, V)

    # Time loop
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

        # Output results
        if ite % (nite ÷ _nout) == 0
            vtk_output && append_vtk(vtk, u, t)
        end
    end

    utrue = exp(-42 * α * t) * u0
    errL2 = norm(var_on_vertices(utrue, mesh) .- var_on_vertices(u, mesh))
    return errL2
end

tfinal = 1
α = 1.0 / 42
nout = 100
for lc in [3e-2, 4e-2, 5e-2, 1e-1, 5e-1]
    @show run(; degree = 1, α, tfinal, Δt = 1e-2, nout, lc, vtk_output = false)
end
end