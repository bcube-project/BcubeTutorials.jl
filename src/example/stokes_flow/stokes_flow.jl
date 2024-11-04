module stokes_flow_API #hide
println("Running Stokes flow API example...") #hide

# # Steady and unsteady Stokes flow example (FE with Taylor-Hood elements)
# In this example, the equations of motion for Stokes flow are solved using $\mathbb{P2}-\mathbb{P1}$ (Taylor-Hood) elements.
#
# # Theory

# # Commented code

# import necessary packages
using Bcube
using LinearAlgebra
using SparseArrays
using WriteVTK
using StaticArrays
using SpecialFunctions
using Test #src

const fspace = :Lagrange
const degree_u = 2
const degquad = 2 * degree_u + 1
const degree_p = 1

const ρ = 1400.0
const μ = 14.0
const ν = μ / ρ

const f = 1.0
const k = sqrt(π * f / ν)

const Δt = 1.0e-3
const finalTime = 2.0

const outputpath = joinpath(@__DIR__, "..", "..", "..", "myout", "stokes")
mkpath(outputpath)

# Function that runs the steady case:
function run_steady()
    # Read 2D mesh
    mesh_path =
        joinpath(@__DIR__, "..", "..", "..", "input", "mesh", "domainTriangle_tri.msh")
    mesh = read_msh(mesh_path)

    fsu = FunctionSpace(fspace, degree_u)
    U_vel = TrialFESpace(
        fsu,
        mesh,
        Dict(
            "edge_right" => SA[0.0, 0.0],
            "edge_left" => SA[0.0, 0.0],
            "edge_top" => SA[1.0, 0.0],
        );
        size = 2,
    )
    V_vel = TestFESpace(U_vel)

    fsp = FunctionSpace(fspace, degree_p)
    U_pre = TrialFESpace(fsp, mesh, Dict("edge_top" => SA[0.0]); size = 1)
    V_pre = TestFESpace(U_pre)

    U = MultiFESpace(U_vel, U_pre)
    V = MultiFESpace(V_vel, V_pre)

    velocity = FEFunction(U_vel, 0.0)
    pressure = FEFunction(U_pre, 0.0)

    #################################################################

    # Define measures for cell
    dΩ = Measure(CellDomain(mesh), degquad)

    f = PhysicalFunction(x -> SA[0.0, 0.0])

    # definition of bilinear and linear forms
    a((u, p), (v, q)) = ∫(μ * ∇(u) ⊡ ∇(v) - tr(∇(v)) * p + tr(∇(u)) * q)dΩ
    l((v, q)) = ∫(f ⋅ v)dΩ

    A = assemble_bilinear(a, U, V)
    L = assemble_linear(l, V)

    Wd = Bcube.assemble_dirichlet_vector(U, V, mesh)

    ## Apply lift
    L = L - A * Wd

    ## Apply homogeneous dirichlet condition
    Bcube.apply_homogeneous_dirichlet_to_vector!(L, U, V, mesh)
    Bcube.apply_dirichlet_to_matrix!(A, U, V, mesh)

    sol = A \ L

    sol = sol .+ Wd

    ϕ = FEFunction(V)

    # Write solution and compare to analytical solution
    set_dof_values!(ϕ, sol)
    velocity, pressure = ϕ

    vars_1 = Dict("Velocity" => velocity)
    Bcube.write_vtk_lagrange(
        joinpath(outputpath, "output_velocity_steady"),
        vars_1,
        mesh,
        U_vel,
    )
    vars_2 = Dict("Pressure" => pressure)
    Bcube.write_vtk_lagrange(
        joinpath(outputpath, "output_pressure_steady"),
        vars_2,
        mesh,
        U_pre,
    )
end

function reference_solution(t, x)
    Uₛ = exp(-k * x[2]) * sin(2.0 * π * f * t - k * x[2])
    T = 2.0 * π * f * t
    Y = x[2] / sqrt(ν / (2.0 * π * f))
    C = 1.0 - im
    Uₜ = imag(
        0.5 *
        exp(-im * T - C * Y / sqrt(2)) *
        erfc(sqrt(0.5 * T) * (C - Y / (T * sqrt(2)))) -
        0.5 *
        exp(-im * T + C * Y / sqrt(2)) *
        erfc(sqrt(0.5 * T) * (C + Y / (T * sqrt(2)))),
    )
    return Uₛ + Uₜ
end

# Function that runs the unsteady case:
function run_unsteady()
    # Read 2D mesh
    mesh_path = joinpath(outputpath, "mesh.msh")
    gen_rectangle_mesh(
        mesh_path,
        :quad;
        nx = 11,
        ny = 51,
        lx = 8.0e-1,
        ly = 5.0e-1,
        xc = 4.0e-1,
        yc = 2.5e-1,
    )
    mesh = read_msh(mesh_path)

    fsu = FunctionSpace(fspace, degree_u)
    U_vel = TrialFESpace(
        fsu,
        mesh,
        Dict("South" => (x, t) -> [sin(2.0 * π * f * t), 0.0]);
        size = 2,
    )
    V_vel = TestFESpace(U_vel)

    fsp = FunctionSpace(fspace, degree_p)
    U_pre = TrialFESpace(fsp, mesh, Dict("North" => SA[0.0]); size = 1)
    V_pre = TestFESpace(U_pre)

    U = MultiFESpace(U_vel, U_pre)
    V = MultiFESpace(V_vel, V_pre)

    velocity = FEFunction(U_vel, 0.0)
    pressure = FEFunction(U_pre, 0.0)

    #################################################################

    # Define measures for cell
    dΩ = Measure(CellDomain(mesh), degquad)

    fv = PhysicalFunction(x -> SA[0.0, 0.0])

    # definition of bilinear and linear forms
    m((u, p), (v, q)) = ∫(ρ * u ⋅ v)dΩ
    a((u, p), (v, q)) = ∫(μ * ∇(u) ⊡ ∇(v) - tr(∇(v)) * p + tr(∇(u)) * q)dΩ
    l((v, q)) = ∫(fv ⋅ v)dΩ

    M = assemble_bilinear(m, U, V)
    A = assemble_bilinear(a, U, V)
    L = assemble_linear(l, V)

    M0 = copy(M)
    A0 = copy(A)

    Bcube.apply_dirichlet_to_matrix!(A, U, V, mesh)
    Bcube.apply_dirichlet_to_matrix!(M, U, V, mesh)

    ϕ = FEFunction(V)

    time = 0.0
    itime = 0
    sol = zero(L)

    while time <= finalTime
        time += Δt
        itime += 1

        println("Time stepping : time = ", time, " / ", finalTime)

        Wd = Bcube.assemble_dirichlet_vector(U, V, mesh, time)

        ## Apply lift
        L = -A0 * Wd

        ## Apply homogeneous dirichlet condition
        Bcube.apply_homogeneous_dirichlet_to_vector!(L, U, V, mesh)

        sol = (M .+ Δt * A) \ (Δt * L + M * sol)

        if itime % 10 == 0
            set_dof_values!(ϕ, sol .+ Wd)
            velocity, pressure = ϕ

            velocity_ref = PhysicalFunction(x -> SA[reference_solution(time, x), 0.0])
            error = velocity - velocity_ref

            vars_1 = Dict(
                "Velocity" => velocity,
                "Velocity_ref" => velocity_ref,
                "error" => error,
            )
            Bcube.write_vtk_lagrange(
                joinpath(outputpath, "output_velocity_unsteady"),
                vars_1,
                mesh,
                U_vel,
                itime,
                time,
            )
            vars_2 = Dict("Pressure" => pressure)
            Bcube.write_vtk_lagrange(
                joinpath(outputpath, "output_pressure_unsteady"),
                vars_2,
                mesh,
                U_pre,
                itime,
                time,
            )
        end
    end
end

#run_steady()
run_unsteady()

end #hide
