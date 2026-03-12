module constrained_poisson_multiplierFESpace_API #hide
println("Running constrained poisson with MultiplierFESpace API example...") #hide

# # Constrained Poisson equation (FE)
# In this example, a Poisson equation with Neumann boundary conditions is solved using a boundary integral constraint.
#
# # Theory
# Consider the following Poisson equation on the unit disk (noted $\Omega$ in this example, its boundary is noted $\Gamma$):
# ```math
# \begin{cases}
# - \Delta u = f \, \, \forall x \in \Omega \\
# \frac{\partial u}{\partial n} = 0 \, \, \forall x \in \Gamma
# \end{cases}
# ```
# Poisson's equation can be written in the form of a minimisation problem:
# ```math
#    \min_{u} J(u) = \frac{1}{2} \int_{\Omega} \nabla u \cdot \nabla u \, dV - \int_{\Omega} f u \, dV
# ```
# There is no unique solution to this problem (adding a constant to any solution will also be a solution).
# Uniqueness can be recovered by adding a constraint to the problem. In this example the following constraint is added:
# ```math
#    \int_{\Gamma} u \, d \gamma = 2 \pi
# ```
# To solve this constrained minimisation problem, the following lagragian is introduced:
# ```math
#    \mathcal{L}(u, \lambda_u) = \frac{1}{2} \int_{\Omega} \nabla u \cdot \nabla u \, dV - \int_{\Omega} f u \, dV + \lambda_u ( \int_{\Gamma} u \, d \gamma - 2 \pi)
# ```
# where $\lambda_u$ is a Lagrange multiplier.
# The first order optimality conditions translate to the problem: find $(u, \lambda_u)$ such that for all $(v, \lambda_v)$:
# ```math
#    \int_{\Omega} \nabla u \cdot \nabla v \, dV + \lambda_u \int_\Gamma v \, d \gamma = \int_\Omega f v \, dV
# ```
# ```math
#    \lambda_v \int_\Gamma u \, d\gamma = 2 \pi \lambda_v
# ```
# This problem can be assembled by introducing a MultiplierFESpace and combining it with the usual FESpace using a MultiFESpace.
# In this example, the manufactured solution $u(x,y)=cos(4\pi(x^2 + y^2))$ is used to test the method.

# # Commented code

# import necessary packages
using Bcube
using BcubeGmsh
using BcubeVTK
using LinearAlgebra
using SparseArrays
using Test #src

const outputpath = joinpath(@__DIR__, "..", "..", "..", "myout", "constrained_poisson")
mkpath(outputpath)

# Read 2D mesh
mesh_path = joinpath(outputpath, "mesh.msh")
BcubeGmsh.gen_disk_mesh(mesh_path; lc = 3.2e-2)
mesh = read_mesh(mesh_path)

# Choose degree and define function space, trial space and test space
const degree = 2
fs = FunctionSpace(:Lagrange, degree)
U = TrialFESpace(fs, mesh)
V = TestFESpace(U)

# Define the multiplier trial space and corresponding test space
# The second argument of MultiplierFESpace specifies the number of
# scalar Lagrange multipliers that are to be used for the problem.
Λᵤ = MultiplierFESpace(mesh, 1)
Λᵥ = TestFESpace(Λᵤ)

# The usual trial FE space and multiplier space are combined into a MultiFESpace
P = MultiFESpace(U, Λᵤ)
Q = MultiFESpace(V, Λᵥ)

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

# Define bilinear and linear forms
function a((u, λᵤ), (v, λᵥ))
    ∫(∇(u) ⋅ ∇(v))dΩ + ∫(side⁻(λᵤ) * side⁻(v))dΓ + ∫(side⁻(λᵥ) * side⁻(u))dΓ
end

# For the time being only functionals in the form of integrals can be assembled.
# A temporary workaround is to put the multiplier in the integral and divide by the volume (the multiplier does not vary in space).
l((v, λᵥ)) = ∫(f * v + 2.0 * π * λᵥ / volume)dΩ

# Assemble to get matrices and vectors
A = assemble_bilinear(a, P, Q)
L = assemble_linear(l, Q)

# Solve problem
sol = A \ L

ϕ = FEFunction(Q)

# Write solution and compare to analytical solution
set_dof_values!(ϕ, sol)
u, λ = ϕ

println("Value of Lagrange multiplier : ", λ.dofValues[1])

u_ref = PhysicalFunction(x -> cos(4.0 * π * (x[1]^2 + x[2]^2)))
error = u_ref - u

vars = Dict("Numerical solution" => u, "Analytical solution" => u_ref, "error" => error)
filename = joinpath(outputpath, "output.pvd")
write_file(filename, mesh, U, vars)

l2(u) = sqrt(sum(Bcube.compute(∫(u ⋅ u)dΩ)))
el2 = l2(error) / l2(u_ref)
tol = 1.e-3
println("L2 relative error : ", el2)

if get(ENV, "TestMode", "false") == "true" #src
    @test el2 < tol                        #src
end                                        #src

end #hide
