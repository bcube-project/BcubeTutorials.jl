module LinearTransport #hide
println("Running linear transport example...") #hide
# # Linear transport (DG)
# In this tutorial, we show how to solve a linear transport equation using a discontinuous-Galerkin
# framework with Bcube.
# # Theory
# In this example, we solve the following linear transport equation using discontinuous elements:
# ```math
# \frac{\partial \phi}{\partial t} + \nabla \cdot (c \phi) = 0
# ```
# where ``c`` is a constant velocity. Using an explicit time scheme, one obtains:
# ```math
# \phi^{n+1} = \phi^n - \Delta t \nabla \cdot (c \phi^n)
# ```
# The corresponding weak form of this equation is:
# ```math
# \int_\Omega \phi^{n+1} v \mathrm{\,d}\Omega = \int_\Omega \phi^n v \mathrm{\,d}\Omega + \Delta t \left[
# \int_\Omega c \phi^n \cdot \nabla v \mathrm{\,d}\Omega - \oint_\Gamma \left( c \phi \cdot n \right) v \mathrm{\,d}\Gamma
# \right]
# ```
# where ``\Gamma = \delta \Omega``. Adopting the discontinuous Galerkin framework, this equation is written in every mesh cell
# ``\Omega_i``. The cell boundary term involves discontinuous quantities and is replaced by a "numerical flux",
# leading to the expression:
# ```math
# \int_{\Omega_i} \phi^{n+1} v \mathrm{\,d}\Omega_i = \int_{\Omega_i} \phi^n v \mathrm{\,d}\Omega_i + \Delta t \left[
# \int_{\Omega_i} c \phi^n \cdot \nabla v \mathrm{\,d}\Omega_i - \oint_{\Gamma_i} F^*(\phi) v \mathrm{\,d} \Gamma_i
# \right]
# ```
# For this example, an upwind flux will be used for ``F^*``. Using a matrix formulation, and noting that the two volumic
# terms are bilinear, the above equation can be written as:
# ```math
# M \phi^{n+1} = (M + \Delta t C) \phi^n - \Delta t f_\Gamma
# ```
# where ``M`` is the mass matrix, ``C`` is the volumic convective term and ``f_\Gamma`` the surfacic flux term.
#
# # Commented code
# Start by importing the necessary packages:
# Load the necessary packages
using Bcube
using BcubeGmsh
using BcubeVTK
using LinearAlgebra

# Define some physical and numerical constant parameters
const degree = 0 # Function-space degree (Taylor(0) = first order Finite Volume)
const c = [1.0, 0.0] # Convection velocity (must be a vector)
const nite = 100 # Number of time iteration(s)
const CFL = 1 # 0.1 for degree 1
const nx = 41 # Number of nodes in the x-direction
const ny = 41 # Number of nodes in the y-direction
const lx = 2.0 # Domain width
const ly = 2.0 # Domain height
const Δt = CFL * min(lx / nx, ly / ny) / norm(c) # Time step

# Then generate the mesh of a rectangle using Gmsh and read it
tmp_path = joinpath(@__DIR__, "..", "..", "myout", "tmp.msh")
BcubeGmsh.gen_rectangle_mesh(
    tmp_path,
    :quad;
    nx = nx,
    ny = ny,
    lx = lx,
    ly = ly,
    xc = 0.0,
    yc = 0.0,
)
mesh = read_mesh(tmp_path)
rm(tmp_path)

# As seen in the previous tutorial, the definition of trial and test spaces needs a mesh and
# a function space. Here, we select Taylor space, and build discontinuous FE spaces with it.
# Then an FEFunction, that will represent our solution, is created.
fs = FunctionSpace(:Taylor, degree)
U = TrialFESpace(fs, mesh, :discontinuous)
V = TestFESpace(U)
u = FEFunction(U)

# Define measures for cell and interior face integrations
Γ = InteriorFaceDomain(mesh)
Γ_in = BoundaryFaceDomain(mesh, "West")
Γ_out = BoundaryFaceDomain(mesh, ("North", "East", "South"))

dΩ = Measure(CellDomain(mesh), 2 * degree + 1)
dΓ = Measure(Γ, 2 * degree + 1)
dΓ_in = Measure(Γ_in, 2 * degree + 1)
dΓ_out = Measure(Γ_out, 2 * degree + 1)

# We will also need the face normals associated to the different face domains.
# Note that this operation is lazy, `nΓ` is just an abstract representation on
# face normals of `Γ`.
nΓ = get_face_normals(Γ)
nΓ_in = get_face_normals(Γ_in)
nΓ_out = get_face_normals(Γ_out)

# Let's move on to the bilinear and linear forms. First, the two easiest ones:
mass(u, v) = ∫(u ⋅ v)dΩ # Mass matrix
conv(u, v) = ∫((c * u) ⋅ ∇(v))dΩ # Volumic convective term

# For the flux term, we first need to define a numerical flux. It is convenient to define it separately
# in a dedicated function. Here is the definition of simple upwind flux.
function upwind(ui, uj, nij)
    cij = c ⋅ nij
    if cij > zero(cij)
        flux = cij * ui
    else
        flux = cij * uj
    end
    flux
end
# We then define the "flux" as the composition of the upwind function and the needed entries: namely the
# solution on the negative side of the face, the solution on the positive face, and the face normal. The
# orientation negative/positive is arbitrary, the only convention is that the face normals are oriented from
# the negative side to the positive side.
flux = upwind ∘ (side⁻(u), side⁺(u), side⁻(nΓ))
l_Γ(v) = ∫(flux * jump(v))dΓ

# Finally, we define what to perform on the "two" boundaries : inlet / oulet.
# On the inlet, we directly impose the flux with a user defined function that depends on the time
# (the input is an oscillating wave).
# On the outlet, we keep our upwind flux but we impose the ghost cell value.
bc_in = t -> PhysicalFunction(x -> c .* cos(3 * x[2]) * sin(4 * t)) # flux
l_Γ_in(v, t) = ∫(side⁻(bc_in(t)) ⋅ side⁻(nΓ_in) * side⁻(v))dΓ_in
flux_out = upwind ∘ (side⁻(u), 0.0, side⁻(nΓ_out))
l_Γ_out(v) = ∫(flux_out * side⁻(v))dΓ_out

# !!! tip
#     For this tutorial, each linear form is defined separately, but it's often more
#     handy to group them into a unique one, either by summing the integral
#     ```julia
#     l(v,t) = ∫(flux * jump(v))dΓ + ∫(side⁻(bc_in(t)) ⋅ side⁻(nΓ_in) * side⁻(v))dΓ_in + ∫(flux_out * side⁻(v))dΓ_out
#     ```
#     or by introducing a closure:
#     ```julia
#     l(v,t) = l_Γ(v) + l_Γ_in(v, t) + l_Γ_out(v)
#     ```

# Assemble the (constant) mass and "convective" matrices. The returned matrices are sparse matrices. Then factorize
# the mass matrix since it will be used in a linear system. We could also define B = M + Δt*C so speed a little bit
# up the computation, but it would require a fixed time step.
M = assemble_bilinear(mass, U, V)
C = assemble_bilinear(conv, U, V)
Mfac = factorize(M)

# Let's also create a vector to avoid allocating it at each time step
# We will then use the in-place version of `assemble_linear` : `assemble_linear!(b_Γ, l, V)`
nd = get_ndofs(V)
b_Γ = zeros(nd) # alternatively: b_Γ = Bcube.allocate_dofs(U)

# Write the initial solution to a file. Note the use of the `discontinuous` flag to write
# a true discontinuous solution with the VTK format
out_dir = joinpath(@__DIR__, "..", "..", "myout", "linear_transport")
mkpath(out_dir) #hide
out_path = joinpath(out_dir, "linear_transport.pvd")
write_file(out_path, mesh, Dict("u" => u), 0, 0.0; discontinuous = true)

# The time loop is trivial : at each time step we compute the linear forms using
# the `assemble_linear!` function, we solve the linear system, update the solution and write
# it to a file.
t = 0.0

for i in 1:nite
    global t

    ## Reset pre-allocated vector
    b_Γ .= 0.0

    ## Compute linear forms
    assemble_linear!(b_Γ, l_Γ, V)
    assemble_linear!(b_Γ, v -> l_Γ_in(v, t), V)
    assemble_linear!(b_Γ, l_Γ_out, V)

    ## Solve the linear system : use `get_dof_values(u)` to get the array of dofs
    x = Mfac\((M + Δt*C)*get_dof_values(u) - Δt*b_Γ)

    ## Update solution: `set_dof_values!(u, x)` to set the dof with an array,
    set_dof_values!(u, x)

    ## Update time
    t += Δt

    ## Append the solution to file
    write_file(
        out_path,
        mesh,
        Dict("u" => u),
        i,
        t;
        discontinuous = true,
        collection_append = true,
    )
end

# And here is an animation of the result:
# ![](../assets/linear_transport.gif)

if get(ENV, "TestMode", "false") == "true"                           #src
    import ..Tester: test_ref                           #src
    test_ref("linear_transport_sol_100ites.jld2", get_dof_values(u)) #src
end                                                                  #src

end #hide
