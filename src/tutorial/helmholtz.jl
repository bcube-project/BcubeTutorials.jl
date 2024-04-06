module Helmholtz #hide
println("Running Helmholtz example...") #hide
# # Helmholtz equation (FE)
# This tutorial shows how to solve the Helmholtz eigen problem with a finite-element approach using Bcube.
# # Theory
# We consider the following Helmholtz equation, representing for instance the acoustic wave propagation with Neuman boundary condition(s):
# ```math
# \begin{cases}
#   \Delta u + \omega^2 u = 0 \\
#   \dfrac{\partial u}{\partial n} = 0 \textrm{  on  } \Gamma
# \end{cases}
# ```
#
# An analytic solution of this equation can be obtained: for a rectangular domain $$\Omega = [0,L_x] \times [0,L_y]$$,
# ```math
# u(x,y) = \cos \left( \frac{k_x \pi}{L_x} x \right) \cos \left( \frac{k_y \pi}{L_y} y \right) \mathrm{~~with~~} k_x,~k_y \in \mathbb{N}
# ```
# with $$\omega^2 = \pi^2 \left( \dfrac{k_x^2}{L_x^2} + \dfrac{k_y^2}{L_y^2} \right)$$
#
# Now, both the finite-element method and the discontinuous Galerkin method requires to write the weak form of the problem:
# ```math
# \int_\Omega (\Delta u \Delta v + \omega^2u)v \mathrm{\,d}\Omega = 0
# ```
# ```math
# - \int_\Omega \nabla u \cdot \nabla v \mathrm{\,d}\Omega
# + \underbrace{\left[ (\nabla u \cdot n) v \right]_\Gamma}_{=0} + \omega^2 \int_\Omega u v \mathrm{\,d} \Omega = 0
# ```
# ```math
# \int_\Omega \nabla u \cdot \nabla v \mathrm{\,d} \Omega = \omega^2 \int_\Omega u v \mathrm{\,d} \Omega
# ```
# Introducing to bilinear forms ``a(u,v)`` and ``b(u,v)`` for respectively the left and right side terms,
# this equation can be written
# ```math
# a(u,v) = \omega^2 b(u,v)
# ```
# This is actually a generalized eigenvalue problem which can be written:
# ```math
# A u = \alpha B u
# ```
# where
# ```math
# A u = \int_\Omega \nabla u \cdot \nabla v \mathrm{\,d} \Omega,~~ B u = \int_\Omega u v \mathrm{\,d} \Omega,~~ \alpha = \omega^2
# ```

# # Uncommented code
# The code below solves the described Helmholtz eigen problem. The code with detailed comments
# is provided in the next section.

# @__PLAIN_CODE__

# # Commented code
# Load the necessary packages
using Bcube
using LinearAlgebra
using Test #src

# Mesh a 2D rectangular domain with quads.
mesh = rectangle_mesh(21, 21)

# Next, create the function space that will be used for the trial and test spaces.
# The Lagrange polynomial space is used here. The degree is set to `1`.
degree = 1
fs = FunctionSpace(:Lagrange, degree)

# The trial space is created from the function space and the mesh. By default, a scalar "continuous"
# FESpace is created. For "discontinuous" ("DG") example, check out the linear transport tutorial.
U = TrialFESpace(fs, mesh)
V = TestFESpace(U)

# Then, define the geometrical dommain on which the operators will apply. For this finite-element example,
# we only need a `CellDomain` (no `FaceDomain`).
domain = CellDomain(mesh)

# Now, integrating on a domain necessitates a "measure", basically a quadrature of given degree.
dΩ = Measure(domain, Quadrature(2 * degree + 1))

# The definition of the two bilinear forms is quite natural. Note that these definitions are lazy,
# no computation is done at this step : the computations will be triggered by the assembly.
a(u, v) = ∫(∇(u) ⋅ ∇(v))dΩ
b(u, v) = ∫(u ⋅ v)dΩ

# To obtain the two sparse matrices corresponding to the discretisation of these bilinear forms,
# simply call the `assemble_bilinear` function, providing the trial and test spaces.
A = assemble_bilinear(a, U, V)
B = assemble_bilinear(b, U, V)

# Compute eigen-values and vectors : we convert to dense matrix to avoid importing additionnal packages,
# but it is quite easy to solve it in a "sparse way".
vp, vecp = eigen(Matrix(A), Matrix(B))

# Display the "first" five eigenvalues:
@show sqrt.(abs.(vp[3:8]))

# Now we can export the solution (the eigenvectors) at nodes of the mesh for several eigenvalues.   #md
# We will restrict to the first 20 eigenvectors. To do so, we will create a `FEFunction` for each   #md
# eigenvector. This `FEFunction` can then be evaluated on the mesh centers, nodes, etc.             #md
ϕ = FEFunction(U)                               #md
nvecs = min(20, get_ndofs(U))                   #md
values = zeros(nnodes(mesh), nvecs)             #md
for ivec in 1:nvecs                             #md
    set_dof_values!(ϕ, vecp[:, ivec])           #md
    values[:, ivec] = var_on_vertices(ϕ, mesh)  #md
end                                             #md

# To write a VTK file, we need to build a dictionnary linking the variable name with its #md
# values and type #md
using WriteVTK #md
outputvtk = "../../myout/helmholtz/rectangle_mesh"                          #md
mkpath(dirname(joinpath(@__DIR__, outputvtk)))                              #src
dict_vars = Dict("u_$i" => (values[:, i], VTKPointData()) for i in 1:nvecs) #md
write_vtk(joinpath(@__DIR__, outputvtk), 0, 0.0, mesh, dict_vars)           #md

# And here is the eigenvector corresponding to the 4th eigenvalue:
# ![](../assets/helmholtz_x21_y21_vp6.png)

if get(ENV, "TestMode", "false") == "true"                  #src
    results = sqrt.(abs.(vp[3:8]))                          #src
    ref_results = [                                         #src
        3.144823462554393,                                  #src
        4.447451992013584,                                  #src
        6.309054755690625,                                  #src
        6.309054755690786,                                  #src
        7.049403274103087,                                  #src
        7.049403274103147,                                  #src
    ]                                                       #src
    @test all(results .≈ ref_results)                       #src
    files2check = outputvtk .* [".pvd", "_00000000.vtu"]    #src
    import ..BcubeTutorialsTests: test_files                #src
    test_files(@__DIR__, files2check)                       #src
end                                                         #src
end #hide
