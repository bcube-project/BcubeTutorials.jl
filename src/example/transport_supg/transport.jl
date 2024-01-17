module LinearTransport #hide
using Bcube
using LinearAlgebra
using WriteVTK
using StaticArrays
using Plots

# function Bcube.idof_by_face_with_bounds(
#     ::Bcube.FunctionSpace{type, 0},
#     shape::Union{Bcube.AbstractShape{1}, Bcube.AbstractShape{2}},
# ) where {type <: Bcube.Lagrange}
#     return ntuple(i -> SA[], Bcube.nfaces(shape))
# end

# Before all, to ease to ease the solution VTK output we will write a
# structure to store the vtk filename and the number of iteration; and a function
# that exports the solution on demand. Note the use of `var_on_nodes_discontinuous`
# to export the solution on the mesh nodes, respecting the discontinuous feature of the
# solution.
mutable struct VtkHandler
    basename::Any
    ite::Any
    mesh::Any
    VtkHandler(basename, mesh) = new(basename, 0, mesh)
end

function append_vtk(vtk, u::Bcube.AbstractFEFunction, uref, t)
    ## Values on center
    values = var_on_nodes_discontinuous(u, vtk.mesh)
    values_ref = var_on_nodes_discontinuous(uref, vtk.mesh)

    ## Write
    Bcube.write_vtk_discontinuous(
        vtk.basename,
        vtk.ite,
        t,
        vtk.mesh,
        Dict("u" => (values, VTKPointData()), "uref" => (values_ref, VTKPointData())),
        1;
        append = vtk.ite > 0,
    )

    ## Update counter
    vtk.ite += 1
end

# First, we define some physical and numerical constant parameters
const degree = 1 # Function-space degree (Taylor(0) = first order Finite Volume)
const nite = 500 # Number of time iteration(s)
const CFL = 0.5 # 0.1 for degree 1
const nx = 101 # Number of nodes in the x-direction
const ny = 41 # Number of nodes in the y-direction
const lx = 1.0 # Domain width
const ly = 2.0 # Domain height
const α = CFL # upwind param

@assert degree >= 1 "Cannot apply Dirichlet when degree = 0!"

# Then generate the mesh of a rectangle using Gmsh and read it
# tmp_path = "tmp.msh"
# gen_rectangle_mesh(tmp_path, :quad; nx = nx, ny = ny, lx = lx, ly = ly, xc = 0.0, yc = 0.0)
# mesh = read_msh(tmp_path)
# rm(tmp_path)
mesh = line_mesh(nx; xmax = lx, names = ("West", "East"))

if spacedim(mesh) == 1
    const c = [1.0]
elseif spacedim(mesh) == 2
    const c = [1.0, 0]
end

const Δt = CFL * min(lx / (nx - 1), ly / (ny - 1)) / norm(c) # Time step

function run()
    # We can now init our `VtkHandler`
    out_dir = joinpath(@__DIR__, "../../../myout/linear_transport")
    mkpath(out_dir) #hide
    vtk = VtkHandler(joinpath(out_dir, "linear_transport"), mesh)

    f_bnd(x, t) = sin(10 * t) # better start with 0 at t=0

    # As seen in the previous tutorial, the definition of trial and test spaces needs a mesh and
    # a function space. Here, we select Taylor space, and build discontinuous FE spaces with it.
    # Then an FEFunction, that will represent our solution, is created.
    fs = FunctionSpace(:Lagrange, degree)
    U = TrialFESpace(fs, mesh, Dict("West" => f_bnd))
    # U = TrialFESpace(fs, mesh, Dict("West" => 1.0))
    V = TestFESpace(U)
    u = FEFunction(U)
    uref = FEFunction(U)

    # Define measures for cell and interior face integrations
    # Γ = InteriorFaceDomain(mesh)
    # Γ_in = BoundaryFaceDomain(mesh, "West")
    # Γ_out = BoundaryFaceDomain(mesh, ("North", "East", "South"))

    dΩ = Measure(CellDomain(mesh), 2 * degree + 1)

    # Compute 'h'
    vol = Bcube.compute(∫(PhysicalFunction(x -> 1))dΩ)
    if spacedim(mesh) == 1
        h = MeshCellData(vol)
    elseif spacedim(mesh) == 2
        h = MeshCellData(sqrt.(vol))
    end

    # supg(v) = v + α * h / (2 * norm(c)) * (c ⋅ ∇(v))
    supg(v) = v + Δt / 2 * c ⋅ ∇(v)
    # supg(v) = v + α * h / (2 * norm(c)) * (c ⋅ ∇(v))

    # Let's move on to the bilinear and linear forms. First, the two easiest ones:
    # a(u, v) = ∫(u ⋅ supg(v))dΩ # Mass matrix
    a(u, v) = ∫(u ⋅ v)dΩ # Mass matrix
    b(u, v) = ∫((c ⋅ ∇(u)) ⋅ supg(v))dΩ # Convection matrix

    # Assemble the (constant) mass matrix. The returned matrix is a sparse matrix. To simplify the
    # tutorial, we will directly compute the inverse mass matrix. But note that way more performant
    # strategies should be employed to solve such a problem (since we don't need the inverse, only the
    # matrix-vector product).
    A = assemble_bilinear(a, U, V)
    B = assemble_bilinear(b, U, V)
    M = I - Δt * inv(Matrix(A)) * B #WARNING : really expensive !!!

    # display(A)
    # display(B)
    # display(M)

    # The time loop is trivial : at each time step we compute the linear forms using
    # the `assemble_` methods, we complete the rhs, perform an explicit step and write
    # the solution.
    p = plot()
    t = 0.0
    apply_dirichlet_to_vector!(u.dofValues, U, V, mesh, t)

    minVal = Inf
    maxVal = -Inf

    for i in 1:nite
        ## Update time
        t += Δt

        ## Update solution
        u.dofValues .= M * u.dofValues

        ## Apply bnd condition : we can also set M[1,:] = [1, 0...]
        apply_dirichlet_to_vector!(u.dofValues, U, V, mesh, t)

        minVal = min(minVal, minimum(u.dofValues))
        maxVal = max(maxVal, maximum(u.dofValues))

        # @show u.dofValues

        projection_l2!(uref, PhysicalFunction(x -> (x[1] - c[1] * t) > 0 ? 0.0 : 1.0), mesh)

        ## Write to file
        append_vtk(vtk, u, uref, t)

        # plot!(lx / (nx - 1) * collect(0:(nx - 1)), u.dofValues)
    end
    plot!(lx / (nx - 1) * collect(0:(nx - 1)), u.dofValues)

    display(p)

    @show minVal, maxVal
end

run()

end #hide
