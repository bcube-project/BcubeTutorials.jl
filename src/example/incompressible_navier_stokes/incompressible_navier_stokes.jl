module incompressible_navier_stokes #hide
println("Running incompressible Navier-Stokes example...") #hide

# # Incompressible Navier-Stokes (FEM) - flow around a cylinder
# In this example, the laminar flow around a cylinder is simulated by solving the incompressible Navier-Stokes equations.
# The equations are discretized and solved using two methods, the projection method and the "mixed form" method. This tutorial is based
# on a [DFG Benchmark](https://wwwold.mathematik.tu-dortmund.de/~featflow/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark2_re100.html) and is also an
# example in [Ferrite.jl](https://ferrite-fem.github.io/Ferrite.jl/stable/tutorials/ns_vs_diffeq/).
#
# # Description of the case
# The simulation domain consists in a channel bounded above and below by two walls and within which a cylinder is placed. The geometrical parameters of the case
# are shown in the Figure below.
#
# ![](../assets/navier_stokes_cylindre.png)
#
# The flow around the cylinder is governed by the incompressible Navier-Stokes equations (here $\rho=1$):
# ```math
#   \partial_t u + u \cdot \nabla u = - \nabla p + \nu \Delta u
# ```
# ```math
#   \nabla \cdot u = 0
# ```
#
# To discretize the equations $\mathbb{P}_2$-$\mathbb{P}_1$ elements are used (Taylor-Hood).
#
# No-slip boundary conditions are applied on all walls. The inlet is given by an imposed parabolic velocity profile which
# is ramped up in time:
# ```math
#   u(t,x,y) = u_{in}(t) \times \left( 4  y (0.41-y)/0.41^2 , 0 \right)
# ```
# where $u_{in}(t) = \textrm{clamp}(t,0.0,1.5)$.
#
# At the outlet, $p=0$ is imposed when using the projection method and a "do-nothing" ($-pn + \nu \nabla u \cdot n = 0$) condition is applied when using the "mixed form" method.
#
# # Code
# Load the necessary packages
const dir = string(@__DIR__, "/") # bcube/example dir
using Bcube
using LinearAlgebra
using BcubeVTK
using BcubeGmsh
using StaticArrays

# Function space (here we shall use Taylor-Hood P2-P1 elements) and quadrature degree.
const fspace = :Lagrange
const degree_u = 2
const degquad = 2 * degree_u + 1
const degree_p = 1

# Input and output paths
const outputpath = joinpath(dir, "..", "..", "..", "myout", "navier_stokes/")
const meshpath = joinpath(dir, "../../../input/mesh/cylinder_navier_stokes_tri.msh")
mkpath(outputpath)

# Kinematic viscosity
const ν = 0.001

# Time step and simulation time
const Δt = 1.0e-3
const finalTime = 3.0

# Function that defines the inlet velocity profile
function inlet_velocity(x, t)
    SA[4.0 * clamp(t, 0.0, 1.5) * x[2] * (0.41 - x[2]) / (0.41 * 0.41), 0.0]
end

# Function that solves the problem using a mixed formalism
function run_unsteady_mixed()
    # Read mesh
    mesh = read_mesh(meshpath)

    # Definition of trial and test function spaces (with associated Dirichlet boundary conditions)
    fsu = FunctionSpace(fspace, degree_u)
    U_vel = TrialFESpace(
        fsu,
        mesh,
        Dict(
            "left" => inlet_velocity,
            "top" => SA[0.0, 0.0],
            "bottom" => SA[0.0, 0.0],
            "cylinder" => SA[0.0, 0.0],
        );
        size = 2,
    )
    V_vel = TestFESpace(U_vel)

    fsp = FunctionSpace(fspace, degree_p)
    U_pre = TrialFESpace(fsp, mesh; size = 1)
    V_pre = TestFESpace(U_pre)

    U = MultiFESpace(U_vel, U_pre)
    V = MultiFESpace(V_vel, V_pre)

    # Define measures for cell
    dΩ = Measure(CellDomain(mesh), degquad)

    # Definition of bilinear and linear forms
    m((u, p), (v, q)) = ∫(u ⋅ v)dΩ
    a((u, p), (v, q)) = ∫(ν * ∇(u) ⊡ ∇(v) - tr(∇(v)) * p + tr(∇(u)) * q)dΩ

    # Assemble and factorize matrices
    M = assemble_bilinear(m, U, V)
    A = assemble_bilinear(a, U, V)

    f_Mtime = M .+ Δt * A
    Bcube.apply_dirichlet_to_matrix!(f_Mtime, U, V, mesh)
    f_Mtime = factorize(f_Mtime)

    ϕ = FEFunction(U)

    # Initial output
    velocity, pressure = ϕ
    vars = Dict("Velocity" => velocity, "Pressure" => pressure)
    filepath = joinpath(outputpath, "output_mixed.pvd")
    write_file(filepath, mesh, vars, 0, 0.0; collection_append = false)

    # Time stepping
    t = 0.0
    itime = 0

    while t <= finalTime
        t += Δt
        itime += 1

        println("Time stepping : time = ", t, " / ", finalTime)

        velocity, pressure = ϕ
        l((v, q)) = -∫((∇(velocity) * velocity) ⋅ v)dΩ
        L = assemble_linear(l, V)

        b = Δt * L + M * get_dof_values(ϕ)
        Bcube.apply_dirichlet_to_vector!(b, U, V, mesh, t)

        ## Compute solution
        sol = f_Mtime \ b

        set_dof_values!(ϕ, sol)

        ## Write outputs
        if itime % 100 == 0
            velocity, pressure = ϕ
            vars = Dict("Velocity" => velocity, "Pressure" => pressure)
            write_file(filepath, mesh, vars, itime, t; collection_append = true)
        end
    end
end

function run_tc_channel()
    # Read mesh
    mesh = rectangle_mesh(20, 20; xmin = 0.0, xmax = 8, ymin = -1.0, ymax = 1.0)

    # Definition of trial and test function spaces (with associated Dirichlet boundary conditions)
    fsu = FunctionSpace(fspace, degree_u)
    U_vel = TrialFESpace(
        fsu,
        mesh,
        Dict("ymin" => SA[0.0, 1.0], "ymax" => SA[0.0, -1.0], "xmin" => SA[0.0, 0.0]);
        size = 2,
    )
    V_vel = TestFESpace(U_vel)

    fsp = FunctionSpace(fspace, degree_p)
    U_pre = TrialFESpace(fsp, mesh; size = 1)
    V_pre = TestFESpace(U_pre)

    U = MultiFESpace(U_vel, U_pre)
    V = MultiFESpace(V_vel, V_pre)

    # Define measures for cell
    dΩ = Measure(CellDomain(mesh), degquad)

    # Definition of bilinear and linear forms
    m((u, p), (v, q)) = ∫(u ⋅ v)dΩ
    a((u, p), (v, q)) = ∫(ν * ∇(u) ⊡ ∇(v) - tr(∇(v)) * p + tr(∇(u)) * q)dΩ

    # Assemble and factorize matrices
    M = assemble_bilinear(m, U, V)
    A = assemble_bilinear(a, U, V)

    f_Mtime = M .+ Δt * A
    Bcube.apply_dirichlet_to_matrix!(f_Mtime, U, V, mesh)
    f_Mtime = factorize(f_Mtime)

    ϕ = FEFunction(U)

    # Initial output
    velocity, pressure = ϕ
    vars = Dict("Velocity" => velocity, "Pressure" => pressure)
    filepath = joinpath(outputpath, "tc_channel.pvd")
    write_file(filepath, mesh, vars, 0, 0.0; collection_append = false)

    # Time stepping
    t = 0.0
    itime = 0

    while t <= finalTime
        t += Δt
        itime += 1

        println("Time stepping : time = ", t, " / ", finalTime)

        velocity, pressure = ϕ
        l((v, q)) = -∫((∇(velocity) * velocity) ⋅ v)dΩ
        L = assemble_linear(l, V)

        b = Δt * L + M * get_dof_values(ϕ)
        Bcube.apply_dirichlet_to_vector!(b, U, V, mesh, t)

        ## Compute solution
        sol = f_Mtime \ b

        set_dof_values!(ϕ, sol)

        ## Write outputs
        if itime % 100 == 0
            velocity, pressure = ϕ
            vars = Dict("Velocity" => velocity, "Pressure" => pressure)
            write_file(filepath, mesh, vars, itime, t; collection_append = true)
        end
    end
end

run_tc_channel()
# run_unsteady_mixed()
end #hide
