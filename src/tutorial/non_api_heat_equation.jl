module NonAPIHeatEqn #hide
println("Running non-api heat equation example...") #hide
# In this tutorial, we manually assemble the bilinear and linear forms associated to the steady heat equation covered before.
# This helps showing what's done with the API behind the scene when writing something like
# ```julia
# a(u, v) = ∫(∇(u) ⋅ ∇(v))dΩ
# A = assemble_bilinear(a, U, V)
# ```
# However, note that the loops written in the present tutorial are not performant compared to the API implementation. It is shown
# here for pedagogic purpose.
#
# It is also important to keep in mind that there are many intermediate functions that can simplify the assembly, even
# we the highest api level is not used. For instance, one may loop on the quadrature points to perform an integral, or
# directly use the `apply_quadrature` function, or even call the `integrate_on_ref` function.
# Hence in this tutorial a few choices are made but there are many ways to proceed.
#
# Recall that the two forms we want to discretize are
# ```math
# a(u,v) = \int_\Omega \eta \nabla(u) \cdot \nabla(v) d\Omega
# l(v) = \int_\Omega qv d\Omega
# ```

# Start by importing the necessary packages.
using Bcube
using LinearAlgebra
using WriteVTK

# First we define some physical and numerical constants
const htc = 100.0 # Heat transfer coefficient (bnd cdt)
const Tr = 268.0 # Recovery temperature (bnd cdt)
const phi = 100.0
const q = 1500.0
const λ = 100.0
const η = λ
const ρCp = 100.0 * 200.0
const degree = 2
const outputpath = joinpath(@__DIR__, "../../myout/non_api_heat_equation/")
mkpath(outputpath) #hide

# Read 2D mesh
mesh_path = joinpath(mktempdir(), "mesh.msh")
gen_rectangle_mesh(mesh_path, :tri; nx = 50, ny = 50)
mesh = read_msh(mesh_path)

# Build function space and associated Trial and Test FE spaces.
# We impose a Dirichlet condition with a temperature of 260K
# on boundary "West"
fs = FunctionSpace(:Lagrange, degree)
U = TrialFESpace(fs, mesh, Dict("West" => 260.0))
V = TestFESpace(U)

# Since we want to perform a "volumic" integration over all the mesh cells,
# we create a `CellDomain` gathering all the cells, as well as a `Measure` for
# the quadrature.
Ω = CellDomain(mesh)
dΩ = Measure(Ω, 2 * degree + 1)

# Allocate the rigidity matrix and the right-hand-side vector
A = zeros(get_ndofs(U), get_ndofs(U))
b = zeros(get_ndofs(U))

# Now the big assembly loop. To sum up, we first need to loop over the mesh cells to perform the volumic 
# integrations. For each cell, we will loop on the dofs / shape functions of the TestFESpace.
# For the bilinear form, we add one layer of loop on the dofs of the TrialFESpace. For each (pair of)
# dof, we perform the integration by applying a quadrature rule (hence a loop) in the reference coordinate
# system (hence a mapping).

# The loop over the mesh cells is performed using a `DomainIterator`. The item obtained by this iterator
# is a `CellInfo`. This type contains the cell type, its nodes indices and coordinates.
for cellinfo in DomainIterator(domain)
    ctype = celltype(cellinfo) # cell entity type
    cshape = shape(ctype) # cell `Shape`
    cnodes = nodes(cellinfo) # cell nodes

    ## Reference -> physical mapping and absolute value of the determinant of the jacobian
    Fmap(ξ) = mapping(cnodes, ctype, ξ)
    J(ξ) = mapping_det_jacobian(nodes, ctype, ξ)

    ## Get the quadrature rule (nodes and weights) associated to this cell
    qrule = QuadratureRule(cshape, get_quadrature(dΩ))
    qnodes = get_nodes(qrule)
    qweights = get_weights(qrule)

    ## Get the global dofs indices of V in the cell (`g` suffix stands for global)
    idofs_g = get_dofs(V, cellindex(cellinfo))

    ## Loop over the dofs of V in the cell
    ## (...)
end

# Let's detail the loop on the dofs of the TestFESpace. There are many ways to access the shape functions of
# the TestFESpace, but remember that the shape functions are always evaluated all at once : you cannot evaluate
# alone the second shape function of `V` in cell `icell`. That's one of the reason why loop on the dofs is not
# performant.
# Below, imagine that the code is included in the previous loop. 
for (idof_l, idof_g) in enumerate(idofs_g)
    λ_V(ξ) = _scalar_shape_functions(fs, cshape, ξ)[idof_l]
    ∇λ_V(ξ) = grad_shape_functions(fs, cshape, ξ)[idof_l, :]

    ## Linear form (rhs) assembly
    ## To apply the quadrature rule, one may directly use the `apply_quadrature` using the `qrule` object
    ## retrieved in the previous loop. Here for the tutorial we choose to loop over the quadrature points.
    for (ωq, ξq) in zip(qweights, qnodes)
        b[idof_g] += q * ωq * J(ξq) * λ_V(ξq)
    end

    jdofs_g = get_dofs(U, cellindex(cellinfo))
    ## Loop over the dofs of U
    ## (...)
end

# In the loop on the dofs of U, we will again apply all the `λV` functions for each couple (idof, jdof), which is
# very bad for performance. Also, note that if the mapping to the reference coordinate system was trivial for the linear
# form, it becomes a bit more complex when it comes to gradients.
# Here again, imagine that the code below is included in the previous loop.
for (jdof_l, jdof_g) in enumerate(jdofs_g)
    λ_U(ξ) = _scalar_shape_functions(fs, cshape, ξ)[jdof_l]
    ∇λ_U(ξ) = grad_shape_functions(fs, cshape, ξ)[jdof_l, :]

    ## This time we directly use the `apply_quadrature` using the `qrule` object
    error("not implemented yet")
end

# Create an affine FE system and solve it using the `AffineFESystem` structure.
# The package `LinearSolve` is used behind the scenes, so different solver may
# be used to invert the system (ex: `solve(...; alg = IterativeSolversJL_GMRES())`)
# The result is a FEFunction (`ϕ`).
# We can interpolate it on mesh centers : the result is named `Tcn`.
sys = AffineFESystem(a, l, U, V)
ϕ = solve(sys)
Tcn = var_on_centers(ϕ, mesh)

# Compute analytical solution for comparison. Apply the analytical solution
# on mesh centers
T_analytical = x -> 260.0 + (q / λ) * x[1] * (1.0 - 0.5 * x[1])
Tca = map(T_analytical, get_cell_centers(mesh))

# Write both the obtained FE solution and the analytical solution to a vtk file.
mkpath(outputpath)
dict_vars =
    Dict("Temperature" => (Tcn, VTKCellData()), "Temperature_a" => (Tca, VTKCellData()))
write_vtk(outputpath * "result_steady_heat_equation", 0, 0.0, mesh, dict_vars)

# Compute and display the error
@show norm(Tcn .- Tca, Inf) / norm(Tca, Inf)
end #hide
