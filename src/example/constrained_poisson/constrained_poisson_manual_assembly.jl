module constrained_poisson_manual_assembly #hide
println("Running constrained poisson API example...") #hide

# # Constrained Poisson equation (FE)
# In this example, a Poisson equation with Neumann boundary conditions is solved using a boundary integral constraint.
#
# # Theory
# Consider the following Poisson equation on the unit disk (noted $$\Omega$$ in this example, its boundary is noted $$\Gamma$$):
# ```math
# \begin{cases}
# - \Delta u = f \, \, \forall x \in \Omega \\
# \frac{\partial u}{\partial n} = 0 \, \, \forall x \in \Gamma
# \end{cases}
# ```
# Poisson's equation can be written in the form of a minimisation problem:
# ```math
#    \min_{u} J(u) = \frac{1}{2} \int_{\Omega} \nabla u . \nabla u \, dV + \int_{\Omega} f u \, dV
# ```
# The discrete version of this problem is:
# ```math
#    \min_{X} J_d(X) = < \frac{1}{2} A X , X >  - < L , X >
# ```
# where $$A$$ is the stiffness matrix corresponding to the bilinear form $$a(u,v) = \int_{\Omega} \nabla u \cdot \nabla v \, dV$$
# and $$L$$ is the right hand side corresponding to the linear form $$l(v) = \int_{\Omega} f v \, dV$$
# There is no unique solution to this problem (adding a constant to any solution will also be a solution).
# Uniqueness can be recovered by adding a constraint to the problem. In this example the following constraint is added:
# ```math
#    \int_{\Gamma} u \, d \gamma = 2 \pi
# ```
# The discrete version of the constraint is: $$<Lc , X > = 2 \pi$$
# where $$Lc$$ is the vector corresponding to the linear form $$l_c(v) = \int_{\Gamma} v \, dV$$.
# To solve this constrained minimisation problem, the following lagragian is introduced:
# ```math
#    L(X, \lambda) = < \frac{1}{2} A X , X >  - < L , X > + \lambda ( < Lc , X > - 2 \pi)
# ```
# where $$\lambda$$ is a Lagrange multiplier.
# The solution of this problem is given by the first order optimality conditions:
# ```math
#    AX + \lambda Lc = L
# ```
# ```math
#    Lc^T X = 2 \pi
# ```
# In this example, the manufactured solution $$u(x,y)=cos(4\pi(x^2 + y^2))$$ is used to test the method.

# # Commented code

# import necessary packages
using Bcube
using BcubeGmsh
using BcubeVTK
using LinearAlgebra
using SparseArrays

const outputpath = joinpath(@__DIR__, "../../../myout/constrained_poisson/")
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

# Define volume and boundary measures
dΩ = Measure(CellDomain(mesh), 2 * degree + 1)
Γ = BoundaryFaceDomain(mesh, ("BORDER",))
dΓ = Measure(Γ, 2 * degree + 1)

# Define solution FE Function
u = FEFunction(U)

# Define source term function (deduced from manufactured solution)
f = PhysicalFunction(
    x ->
        64.0 * π^2 * (x[1]^2 + x[2]^2) * cos(4.0 * π * (x[1]^2 + x[2]^2)) +
        16.0 * π * sin(4.0 * π * (x[1]^2 + x[2]^2)),
)

# Define bilinear and linear forms
a(u, v) = ∫(∇(u) ⋅ ∇(v))dΩ
l(v) = ∫(f * v)dΩ
lc(v) = ∫(side⁻(v))dΓ

# Assemble to get matrices and vectors
A  = assemble_bilinear(a, U, V)
L  = assemble_linear(l, V)
Lc = assemble_linear(lc, V)

# Build augmented problem
n = size(L)[1]

M = spzeros(n + 1, n + 1)
B = zeros(n + 1)

M[1:n, 1:n] .= A[1:n, 1:n]
M[n + 1, 1:n] .= Lc[:]
M[1:n, n + 1] .= Lc[:]
B[1:n] .= L[1:n]
B[n + 1] = 2.0 * π

# Solve problem
sol = M \ B

# Write solution and compare to analytical solution
set_dof_values!(u, sol[1:n])
λ = sol[n + 1]

println("Value of Lagrange multiplier : ", λ)

u_ref = PhysicalFunction(x -> cos(4.0 * π * (x[1]^2 + x[2]^2)))
error = u_ref - u

vars = Dict("Numerical solution" => u, "Analytical solution" => u_ref, "error" => error)
filename = joinpath(outputpath, "output.pvd")
write_file(filename, mesh, U, vars)

l2(u) = sqrt(sum(Bcube.compute(∫(u ⋅ u)dΩ)))
el2 = l2(error) / l2(u_ref)
tol = 1.e-3
println("L2 relative error : ", el2)
@assert el2 < tol

end #hide
