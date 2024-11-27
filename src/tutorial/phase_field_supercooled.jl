module PhaseFieldSupercooled #hide
println("Running phase field supercooled equation example...") #hide

# # Phase field model - solidification of a liquid in supercooled state
# In this tutorial, a coupled system of two unsteady equations is solved using finite elements
# and an imex time scheme. This tutorial doesn't introduce `MultiFESpace`, check the "euler" example
# for this. Warning : this file is currently quite long to run (a few minutes).
#
# # Theory
# This case is taken from: Kobayashi, R. (1993). Modeling and numerical simulations of dendritic crystal growth. Physica D: Nonlinear Phenomena, 63(3-4), 410-423.
# In particular, the variables of the problem are denoted in the same way ($p$ for the phase indicator and $T$ for temperature).
# Consider a rectangular domain $$\Omega = [0, L_x] \times [0, L_y]$$ on which we wish to solve the following equations:
# ```math
#   \tau \partial_t p = \epsilon^2 \Delta p + p (1-p)(p - \frac{1}{2} + m(T))
# ```
# ```math
#   \partial_t T = \Delta T + K \partial_t p
# ```
# where $m(T) = \frac{\alpha}{\pi} atan \left[ \gamma (T_e - T) \right]$.
# This set of equations represents the solidification of a liquid in a supercooled state. Here $T$ is a dimensionless temperature and $p$ is the solid volume fraction.
# Lagrange finite elements are used to discretize both equations. Time marching is performed with a forward Euler scheme for the first equation and a backward Euler scheme for the second one.
#
# To initiate the solidification process, a Dirichlet boundary condition ($p=1$,$T=1$) is applied at $x=0$ ("West" boundary).
#
# # Code
# Load the necessary packages
using Bcube
using BcubeVTK
using LinearAlgebra
using Random

Random.seed!(1234) # to obtain reproductible results

# Define some physical and numerical constants, as well as the `g` function
# appearing in the problem definition.
const ε = 0.01
const τ = 0.0003
const α = 0.9
const γ = 10.0
const K = 1.6
const Te = 1.0
const β = 0.0 # noise amplitude, original value : 0.01
const Δt = 0.0001 # time step
const totalTime = 1.0 # original value : 1
const nout = 50 # Number of iterations to skip before writing file
const degree = 1 # function space degree
const lx = 3.0
const ly = 1.0
const nx = 100
const ny = 20
const out_dir = joinpath(@__DIR__, "../../myout/phase_field_supercooled") # output directory
mkpath(out_dir) #hide

# Read the mesh using `gmsh`
const mesh_path = joinpath(@__DIR__, "../../input/mesh/domainPhaseField_tri.msh")

# Here we define the main algorithm in a function
# to avoid performance penalty (see [Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-tips))
function main()
    # read the mesh file
    mesh = read_msh(mesh_path)

    # Noise function : random between [-1/2,1/2]
    χ = MeshCellData(rand(ncells(mesh)) .- 0.5)

    # Build the function space and the FE Spaces. The two unknowns will share the
    # same FE spaces for this tutorial. Note the way we specify the Dirichlet condition
    # in the definition of `U`.
    fs = FunctionSpace(:Lagrange, degree)
    U = TrialFESpace(fs, mesh, Dict("West" => (x, t) -> 1.0))
    V = TestFESpace(U)

    # Build FE functions
    ϕ = FEFunction(U)
    T = FEFunction(U)

    # Define measures for cell integration
    dΩ = Measure(CellDomain(mesh), 2 * degree + 1)

    g(T) = (α / π) * atan(γ * (Te - T))

    # Define bilinear and linear forms.
    # As the linear form is assembled at each step,
    # we use the composition operator to exploit
    # the full performance of Bcube.
    a(u, v) = ∫(∇(u) ⋅ ∇(v))dΩ
    m(u, v) = ∫(u ⋅ v)dΩ
    f_l(ϕ, T, χ, v) = v * ϕ * (1.0 - ϕ) * (ϕ - 0.5 + g(T) + β * χ)
    l(v) = ∫(f_l ∘ (ϕ, T, χ, v))dΩ

    # Assemble the two constant matrices
    A = assemble_bilinear(a, U, V)
    M = assemble_bilinear(m, U, V)

    # Create iterative matrices
    C_ϕ = M + Δt / τ * ε^2 * A
    C_T = M + Δt * A

    # Apply Dirichlet conditions.
    # For this example, we don't use a lifting method to impose the Dirichlet, but `d`
    # is used to initialize the solution.
    d = assemble_dirichlet_vector(U, V, mesh)
    apply_dirichlet_to_matrix!((C_ϕ, C_T), U, V, mesh)

    # Init solution and write it to a VTK file
    set_dof_values!(ϕ, d)
    set_dof_values!(T, d)

    filepath = joinpath(out_dir, "result_phaseField_imex_1space.pvd")
    dict_vars = Dict("Temperature" => T, "Phi" => ϕ)
    write_file(filepath, mesh, dict_vars, 0, 0.0)

    # Factorize and allocate some vectors to increase performance
    C_ϕ = factorize(C_ϕ)
    C_T = factorize(C_T)
    L = zero(d)
    rhs = zero(d)
    ϕ_new = zero(d)

    # Time loop (imex time integration)
    t = 0.0
    itime = 0
    while t <= totalTime
        t += Δt
        itime += 1
        @show t, totalTime

        ## Integrate equation on ϕ
        L .= 0.0 # reset L
        assemble_linear!(L, l, V)
        rhs .= M * get_dof_values(ϕ) .+ (Δt / τ) .* L
        apply_dirichlet_to_vector!(rhs, U, V, mesh)
        ϕ_new .= C_ϕ \ rhs

        ## Integrate equation on T
        rhs .= M * (get_dof_values(T) .+ K .* (ϕ_new .- get_dof_values(ϕ)))
        apply_dirichlet_to_vector!(rhs, U, V, mesh)

        ## Update solution
        set_dof_values!(ϕ, ϕ_new)
        set_dof_values!(T, C_T \ rhs)

        ## write solution in vtk format
        if itime % nout == 0
            write_vtk(filepath, mesh, dict_vars; itime, t, collection_append = true)
        end
    end
end

# run simulation
main()

# And here is an animation of the result:
# ![](../assets/phase-field-supercooled-rectangle.gif)

end #hide
