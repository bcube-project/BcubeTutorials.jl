module constrained_poisson_API #hide
println("Running constrained poisson API example...") #hide

# # Constrained Poisson equation (FE)
# In this example, a Poisson equation with Neumann boundary conditions is solved using a boundary integral constraint.
#
# # Theory
# Consider the following Poisson equation on the unit disk (noted $$\Omega$$ in this example, its boundary is noted $$\Gamma$$):
# ```math
#    - \Delta u = f \, \, \forall x \in \Omega
# ```
# ```math
#    \frac{\partial u}{\partial n} = 0 \, \, \forall x \in \Gamma
# ```
# Poisson's equation can be written in the form of a minimisation problem:
# ```math
#    \min_{u} J(u) = \frac{1}{2} \int_{\Omega} \nabla u . \nabla u \, dV - \int_{\Omega} f u \, dV
# ```
# There is no unique solution to this problem (adding a constant to any solution will also be a solution).
# Uniqueness can be recovered by adding a constraint to the problem. In this example the following constraint is added:
# ```math
#    \int_{\Gamma} u \, d \gamma = 2 \pi
# ```
# To solve this constrained minimisation problem, the following lagragian is introduced:
# ```math
#    \mathcal{L}(u, \lambda_u) = \frac{1}{2} \int_{\Omega} \nabla u . \nabla u \, dV - \int_{\Omega} f u \, dV + \lambda_u ( \int_{\Gamma} u \, d \gamma - 2 \pi)
# ```
# where $$\lambda_u$$ is a Lagrange multiplier.
# The first order optimality conditions translate to the problem: find (u, \lambda_u) such that for all (v, \lambda_v):
# ```math
#    \int_{\Omega} \nabla u . \nabla v \, dV + \lambda_u \int_\Gamma v \, d \gamma = \int_\Omega f v \, dV
# ```
# ```math
#    \lamda_v \int_\Gamma u \, d\gamma = 2 \pi \lambda_v
# ```
# This problem can be assembled by introducing a MultiplierFESpace and combining it with the usual FESpace using a MultiFESpace. 
# In this example, the manufactured solution $$u(x,y)=cos(4\pi(x^2 + y^2))$$ is used to test the method.

# import necessary packages
using Bcube
using LinearAlgebra
using SparseArrays
using WriteVTK

const outputpath = joinpath(@__DIR__, "../../../myout/constrained_poisson/")
mkpath(outputpath)

# Read 2D mesh
mesh_path = joinpath(outputpath, "mesh.msh")
gen_disk_mesh(mesh_path; lc = 3.2e-2)
mesh = read_msh(mesh_path)

# Choose degree and define function space, trial space and test space
const degree = 2
fs = FunctionSpace(:Lagrange, degree)
U = TrialFESpace(fs, mesh)
V = TestFESpace(U)

# Define the multiplier trial space and corresponding test space
Λᵤ = MultiplierFESpace(mesh, 1)
Λᵥ = TestFESpace(Λᵤ)

# The usual trial FE space and multiplier space are combined into a MultiFESpace
TrialMultiFESpace = MultiFESpace(U, Λᵤ; arrayOfStruct = false)
TestMultiFESpace = MultiFESpace(V, Λᵥ; arrayOfStruct = false)

# Define volume and boundary measures
dΩ = Measure(CellDomain(mesh), 2 * degree + 1)
Γ = BoundaryFaceDomain(mesh, ("BORDER",))
dΓ = Measure(Γ, 2 * degree + 1)

# Define solution FE Function
ϕ = FEFunction(U)

# Define source term function (deduced from manufactured solution)
f = PhysicalFunction(
    x ->
        64.0 * π^2 * (x[1]^2 + x[2]^2) * cos(4.0 * π * (x[1]^2 + x[2]^2)) +
        16.0 * π * sin(4.0 * π * (x[1]^2 + x[2]^2)),
)

volume = sum(Bcube.compute(∫(PhysicalFunction(x -> 1.0))dΩ))
@show volume

# Define bilinear and linear forms
function a((u, λᵤ), (v, λᵥ))
    ∫(∇(u) ⋅ ∇(v))dΩ + ∫(side⁻(λᵤ) * side⁻(v))dΓ + ∫(side⁻(λᵥ) * side⁻(u))dΓ
end

# For the time being only functionals in the form of integrals can be assembled.
# A temporary workaround is to put the multiplier in the integral and divide by the volume (the multiplier does not vary in space).
l((v, λᵥ)) = ∫(f * v + 2.0 * π * λᵥ / volume)dΩ

# Assemble to get matrices and vectors
A = assemble_bilinear(a, TrialMultiFESpace, TestMultiFESpace)
L = assemble_linear(l, TestMultiFESpace)

########################### DEBUG ###############################

# Define bilinear and linear forms
a2(u, v) = ∫(∇(u) ⋅ ∇(v))dΩ
l2(v) = ∫(f * v)dΩ
lc(v) = ∫(side⁻(v))dΓ

# Assemble to get matrices and vectors
A2 = assemble_bilinear(a2, U, V)
L2 = assemble_linear(l2, V)
Lc = assemble_linear(lc, V)

# Build augmented problem
n = size(L2)[1]

M = spzeros(n + 1, n + 1)
B = zeros(n + 1)

M[1:n, 1:n] .= A2[1:n, 1:n]
M[n + 1, 1:n] .= Lc[:]
M[1:n, n + 1] .= Lc[:]
B[1:n] .= L2[1:n]
B[n + 1] = 2.0 * π

#################################################################

display(abs.(A[(n + 1), n] .- M[(n + 1), n]))

# Solve problem
sol = A \ L

ϕ = FEFunction(TrialMultiFESpace)

# Write solution and compare to analytical solution
for i in 1:n
    d = abs(2.0 * π - L[i])
    if d < 1.0e-10
        println(i, " ", L[i], " ", n, " ", sol[i])
    end
end
set_dof_values!(ϕ, sol)
u, λ = ϕ

println(" Value of Lagrange multiplier : ", λ.dofValues[1])

ϕₑ = FEFunction(U)

projection_l2!(ϕₑ, PhysicalFunction(x -> cos(4.0 * π * (x[1]^2 + x[2]^2))), mesh)

Un = var_on_vertices(u, mesh)
Ue = var_on_vertices(ϕₑ, mesh)
mkpath(outputpath)
dict_vars = Dict(
    "Numerical Solution" => (Un, VTKPointData()),
    "Analytical solution" => (Ue, VTKPointData()),
)
write_vtk(
    outputpath * "result_constrained_poisson_equation_multiplierFESpace",
    0,
    0.0,
    mesh,
    dict_vars,
)

error = norm(Un .- Ue, Inf) / norm(Ue, Inf)
println(" Error : ", error)

end #hide
