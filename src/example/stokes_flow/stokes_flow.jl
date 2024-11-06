module stokes_flow_API #hide
println("Running Stokes flow API example...") #hide

# # Steady and unsteady Stokes flow example (FE with Taylor-Hood elements)
# This example shows how to solve the equations of motion for Stokes flow with $\mathbb{P2}-\mathbb{P1}$ (Taylor-Hood) elements using Bcube.
#
# # Theory (Steady Stokes flow - Moffat vortices)
# Let's first consider a Stokes flow within a wedge of angle $\theta = 28^\circ$ and height $H=1m$. A velocity of $1m/s$ is imposed on the top boundary (the lid) as well as zero pressure. 
# On the left and right edges a no-slip boundary condition is applied. 
#
# ![](../assets/Stokes_steady_Moffat_vortices_1.png)
#
# The set of equations are:
# ```math
# \begin{cases}
#   0 = \nabla \cdot \left( -p \mathbb{I} + \mu \nabla u \right) + f  \textrm{  in } \Omega \\
#   \nabla \cdot u = 0 \textrm{  in } \Omega  \\
#   u = (0,0)  \textrm{  on } \Gamma_L \cup \Gamma_R \\
#   u = (1,0) \textrm{ and } p=0 \textrm{  on } \Gamma_T
# \end{cases}
# ```
# Density $\rho$ is taken equal to $1400kg.m^{-3}$ and dynamic viscosity $\mu$ equal to $1.4 Pa.s$. The Reynolds number based on the height of the wedge is therefore equal to $1000$.
# The weak form of this problem is given by:
# ```math
# \begin{cases}
#   \textrm{ Find } u \in V \textrm{ and } p \in Q \textrm{ such that } \forall v \in V, \, \forall q \in Q\\
#   \int_\Omega \left( \mu \nabla u : \nabla v - p \nabla \cdot v \right) \, dx  = \int_\Omega f \cdot v \, dx  \\
#   \int_\Omega q \nabla \cdot u\, dx  = 0  
# \end{cases}
# ```
# where $V$ and $Q$ are suitable function spaces. This can be recast in "mixed" formalism:
# ```math
# \begin{cases}
#   \textrm{ Find } (u,p) \in W=V \times Q  \textrm{ such that } \forall (v,q) \in W\\
#   a((u,p),(v,q)) = l((v,q))
# \end{cases}
# ``` 
# where $a((u,p),(v,q)) = \int_\Omega \left( \mu \nabla u : \nabla v - p \nabla \cdot v +q \nabla \cdot u \right)\, dx$ and
# $l((v,q))=\int_\Omega f \cdot v \, dx$
# 
# Taylor-Hood ($\mathbb{P2}-\mathbb{P1}$) elements are used to discretize the weak form of the problem which leads to the linear system:
# ```math
# \begin{bmatrix}
#  A & B^T \\
#  B & 0 
# \end{bmatrix} 
# \begin{bmatrix} 
# U \\
# P
# \end{bmatrix}
# = L
# ``` 
#
# # Commented code

# import necessary packages
using Bcube
using LinearAlgebra
using SparseArrays
using WriteVTK
using StaticArrays
using SpecialFunctions
using Test #src

# Function space: P2 for velocity and P1 for pressure
const fspace = :Lagrange
const degree_u = 2
const degquad = 2 * degree_u + 1
const degree_p = 1

# Material parameters
const ρ = 1400.0
const μ = 14.0
const ν = μ / ρ

# Parameters for the oscillating plate case
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

    # Define MultiFESpace (mixed formalism)
    U = MultiFESpace(U_vel, U_pre)
    V = MultiFESpace(V_vel, V_pre)

    velocity = FEFunction(U_vel, 0.0)
    pressure = FEFunction(U_pre, 0.0)

    # Define measures for cell
    dΩ = Measure(CellDomain(mesh), degquad)

    # volume force is equal to zero
    f = PhysicalFunction(x -> SA[0.0, 0.0])

    # definition of bilinear and linear forms
    a((u, p), (v, q)) = ∫(μ * ∇(u) ⊡ ∇(v) - tr(∇(v)) * p + tr(∇(u)) * q)dΩ
    l((v, q)) = ∫(f ⋅ v)dΩ

    # assemble matrix and RHS
    A = assemble_bilinear(a, U, V)
    L = assemble_linear(l, V)

    # build Dirichlet vector to apply lift
    Wd = Bcube.assemble_dirichlet_vector(U, V, mesh)

    # Apply lift
    L = L - A * Wd

    # Apply homogeneous dirichlet condition
    Bcube.apply_homogeneous_dirichlet_to_vector!(L, U, V, mesh)
    Bcube.apply_dirichlet_to_matrix!(A, U, V, mesh)

    # Solve linear system
    sol = A \ L

    sol = sol .+ Wd

    ϕ = FEFunction(V)

    # Write solution
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
# The obtained solution captures the Moffat vortices topology of the flow
# ![](../assets/Stokes_flow_Moffat_vortices.png)

# # Theory (Unsteady Stokes flow - oscillating plate)
# Let's now consider a Stokes flow over an oscillating flat plate of length $0.1m$ and height $0.5m$. An oscillating velocity $(\sin(2\pi f t),0)$ is imposed on the bottom boundary as well as zero pressure. 
# On all the other boundaries a "do nothing" condition is applied: $-pn + \mu \frac{\partial u}{\partial n} = 0$. 
# The set of equations are:
# ```math
# \begin{cases}
#   \rho \frac{\partial u}{\partial t} = \nabla \cdot \left( -p \mathbb{I} + \mu \nabla u \right) + f  \textrm{  in } \Omega \\
#   \nabla \cdot u = 0 \textrm{  in } \Omega  \\
#   u = (0,0)  \textrm{  on } \Gamma_L \cup \Gamma_R \\
#   u = (1,0) \textrm{ and } p=0 \textrm{  on } \Gamma_T
# \end{cases}
# ```
# The weak form of this problem is given by:
# ```math
# \begin{cases}
#   \textrm{ Find } (u,p) \in W=V \times Q  \textrm{ such that } \forall (v,q) \in W\\
#   m((u,p),(v,q))+a((u,p),(v,q)) = l((v,q))
# \end{cases}
# ``` 
# where $m((u,p),(v,q)) = \int_\Omega \rho u \cdot v \, dx$ and
# $l((v,q))=\int_\Omega f \cdot v \, dx$
# $\mathbb{P2}-\mathbb{P1}$ (Taylor-Hood) elements are used to discretize the weak form of the problem which leads to the linear system:
# ```math
# \begin{bmatrix}
#  M & 0 \\
#  0 & 0 
# \end{bmatrix} 
# \begin{bmatrix} 
# \dot{U} \\
# \dot{P}
# \end{bmatrix} +
# \begin{bmatrix}
#  A & B^T \\
#  B & 0 
# \end{bmatrix} 
# \begin{bmatrix} 
# U \\
# P
# \end{bmatrix}
# = L
# ``` 
#
# # Commented code

# Reference solution for unsteady case (Panton, R. L. (2013). Incompressible Flow, 4th edition, pages 228-236)
# This solution is valid for a semi-infinite domain.  
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
        lx = 1.0e-1,
        ly = 5.0e-1,
        xc = 0.5e-1,
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
    U_pre = TrialFESpace(fsp, mesh, Dict("South" => SA[0.0]); size = 1)
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

        # Apply lift
        L = -A0 * Wd

        # Apply homogeneous dirichlet condition
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
# The obtained solution compares well with the reference solution
# ![](../assets/Stokes_flow_unsteady.gif)

println("----- Running steady case -----")
run_steady()
println("-------------------------------")
println("")
println("----- Running unsteady case -----")
run_unsteady()
println("-------------------------------")
println("")
end #hide
