module TransportHypersurface #hide
# # Transport equation on hypersurfaces
# In this file, a linear transport equation is solved on several hypersurfaces. The scalar
# transport equation is
# ```math
#   \partial_t u + c \cdot \nabla_\Gamma u = 0
# ```
# where $c$ is the (tangential) transport velocity (always divergent-free in these examples); and
# $\nabla_\Gamma$ is the tangential gradient.
# The vector transport equation is
# ```math
#  \partial_t u + \nabla_\Gamma (c \otimes u) = 0
# ```
#
# Here is an example of the transport of a scalar bump on a torus meshed with P2 triangles:
#
# ![](../assets/transport-torus-mesh2-degree1.gif)
using Plots
using Bcube
using StaticArrays
using LinearAlgebra
using Printf
using WriteVTK
using DelimitedFiles
using Random
using ProgressMeter
using Profile
using BenchmarkTools

const out_dir = joinpath(@__DIR__, "../../../myout/transport_hypersurface")
rm(out_dir; force = true, recursive = true)
mkpath(out_dir)

""" Norm alias for AbstractLazy """
mynorm(a) = sqrt(a ⋅ a)

""" Tangential divergence operator """
divₛ(a) = tr(∇ₛ(a))

mutable struct VtkHandler
    basename::Any
    ite::Any
    mesh::Any
    dΩ::Any
    U::Any
    θ::Any
    θ_centers::Any
    θ_vertices::Any
    c::Any
    c_centers::Any
    c_vertices::Any
    ν::Any
    ν_centers::Any
    ν_vertices::Any
    function VtkHandler(basename, dΩ, U, c)
        @info "Writing to $basename.vtu"

        mesh = get_mesh(get_domain(dΩ))
        θ = PhysicalFunction(x -> atan(x[2], x[1]))
        θ_centers = var_on_centers(θ, mesh)
        θ_vertices = var_on_vertices(θ, mesh)

        ν = Bcube.CellNormal(mesh)
        ν_centers = transpose(var_on_centers(ν, mesh))
        ν_vertices = transpose(var_on_vertices(ν, mesh))

        c_centers = transpose(var_on_centers(c, mesh))
        c_vertices = transpose(var_on_vertices(c, mesh))

        new(
            basename,
            0,
            mesh,
            dΩ,
            U,
            θ,
            θ_centers,
            θ_vertices,
            c,
            c_centers,
            c_vertices,
            ν,
            ν_centers,
            ν_vertices,
        )
        ## new(basename, 0, mesh, [atan(n.x[2], n.x[1]) for n in Bcube.get_nodes(mesh)])
    end
end

function plot_solution(i, t, u, mesh, xcenters, ycenters, xnodes, ynodes)
    ## Build animation
    uplot = var_on_centers(u, mesh)
    lmax = 1.5

    plt = plot(
        [xnodes..., xnodes[1]],
        [ynodes..., ynodes[1]];
        aspect_ratio = :equal,
        primary = false,
        xlim = (-lmax, lmax),
        ylim = (-lmax, lmax),
    )
    annotate!(0, 0.25, "i  = $i")
    annotate!(0, 0, @sprintf "t = %.2e" t)

    if ndims(uplot) == 1
        scatter!(xcenters, ycenters; marker_z = uplot, label = "u", clims = (-1, 1))
        annotate!(0, -0.25, @sprintf "|u|_max = %.2e" maximum(abs.(uplot)))
    elseif ndims(uplot) == 2
        L = maximum(x -> norm(x), eachrow(uplot))
        scale = 0.75
        uplot .*= scale
        quiver!(
            xcenters,
            ycenters;
            quiver = (uplot[:, 1], uplot[:, 2]),
            label = "u",
            xlabel = "x",
            ylabel = "y",
        )
        annotate!(0, -0.25, @sprintf "|u|_max = %.2e" L)
    end

    return plt
end

"""
    rotMat(θx, θy, θz)

Build the 3D rotation matrix from the angles given for each axis.
"""
function rotMat(θx, θy, θz)
    Rx = @SMatrix[
        1.0 0.0 0.0
        0.0 cos(θx) sin(θx)
        0.0 (-sin(θx)) cos(θx)
    ]
    Ry = @SMatrix[
        cos(θy) 0.0 (-sin(θy))
        0.0 1.0 0.0
        sin(θy) 0.0 cos(θy)
    ]
    Rz = @SMatrix[
        cos(θz) sin(θz) 0.0
        -sin(θz) cos(θz) 0.0
        0.0 0.0 1.0
    ]
    return Rx * Ry * Rz
end

"""
Linear transport of a scalar quantity on a circle

nrot = number of "rotations" to run (in time)
volumic_bilinear = true means the volumic term is assembled once with a bilinear assembly while false means
    a linear assembly (and hence the limiter is applied)
"""
function scalar_circle(;
    degree,
    CFL,
    nθ,
    nrot = 2,
    nout = 100,
    nitemax = Int(1e9),
    volumic_bilinear = false,
    isLimiterActive = true,
)
    function append_vtk(vtk, u::Bcube.AbstractFEFunction, lim_u, u_mean, t)
        ## Build animation
        values_vertices = var_on_vertices(u, mesh)
        values_centers = var_on_centers(u, mesh)
        ## values_nodes = var_on_nodes_discontinuous(u, mesh, degree)
        θ_centers = var_on_centers(vtk.θ, mesh)
        θ_vertices = var_on_vertices(vtk.θ, mesh)

        ## Write
        Bcube.write_vtk(
            vtk.basename,
            vtk.ite,
            t,
            vtk.mesh,
            Dict(
                "u_centers" => (values_centers, VTKCellData()),
                "u_vertices" => (values_vertices, VTKPointData()),
                "lim_u" => (get_values(lim_u), VTKCellData()),
                "u_mean" => (get_values(u_mean), VTKCellData()),
                "θ_centers" => (θ_centers, VTKCellData()),
                "θ_vertices" => (θ_vertices, VTKPointData()),
            ),
            ;
            append = vtk.ite > 0,
        )

        ## Update counter
        vtk.ite += 1
    end

    ## Settings
    radius = 1.0
    C = 1.0 ## velocity norm
    tmax = nrot * 2π / C

    ## Mesh
    qOrder = 2 * degree + 1 ## we shall use Gauss-Lobatto for degree > 0, but in 1D this is ok
    mesh = circle_mesh(nθ; radius = radius, order = 1)
    dΩ = Measure(CellDomain(mesh), qOrder)
    Γ = InteriorFaceDomain(mesh)
    dΓ = Measure(Γ, qOrder)
    nΓ = get_face_normals(Γ)

    ## FESpace
    fs = FunctionSpace(:Lagrange, degree)
    U = TrialFESpace(fs, mesh, :discontinuous)
    V = TestFESpace(U)

    ## Transport velocity
    _c = PhysicalFunction(x -> C * SA[-x[2], x[1]] / radius)
    P = Bcube.tangential_projector(mesh)
    c = (x -> C * normalize(x)) ∘ (P * _c) ## useless in theory since velocity is already tangent

    ## Find quadrature weight (mesh is composed of a unique "shape" so first element is enough)
    quad = Bcube.get_quadrature(dΩ)
    s = Bcube.shape(Bcube.cells(mesh)[1])
    qrule = Bcube.QuadratureRule(s, quad)
    ω_quad = degree > 0 ? Bcube.get_weights(qrule)[1] : 1.0

    ## Time step and else
    dl = 2π * radius / nθ ## analytic length
    dl = 2 * radius * sin(2π / nθ / 2) ## discretized length
    Δt = CFL * dl * ω_quad / C / (2 * degree + 1)
    t = 0.0
    nite = min(floor(Int, tmax / Δt), nitemax)
    @show nite
    _nout = min(nite, nout)

    @show dl
    @show Δt

    ## Limitation
    DMPrelax = 0.0 * dl

    ## Output
    tail = isLimiterActive ? "lim" : "nolim"
    filename = "scalar-on-circle-d$(degree)-$(tail)"
    vtk = VtkHandler(joinpath(out_dir, filename), dΩ, U, c)
    dofOverTime = zeros(nite + 1, 2) ## t, u
    i_dof_out = 1

    ## FEFunction and "boundary / source" condition
    u = FEFunction(U)
    if false
        u.dofValues[1] = 1.0
    else
        projection_l2!(u, PhysicalFunction(x -> cos(atan(x[2], x[1]))), mesh)
    end

    ## Forms
    m(u, v) = ∫(u ⋅ v)dΩ ## Mass matrix
    a_Ω(u, v) = ∫(u * (c ⋅ ∇ₛ(v)))dΩ ## bilinear volumic convective term

    function upwind(ui, uj, ci, nij)
        cij = ci ⋅ nij
        if cij > zero(cij)
            fij = cij * ui
        else
            fij = cij * uj
        end
        fij
    end

    ## Mass
    M = assemble_bilinear(m, U, V)
    invM = inv(Matrix(M)) ##WARNING : really expensive !!!
    if volumic_bilinear
        K = assemble_bilinear(a_Ω, U, V)
        invMK = invM * K
    end

    ## Anim
    anim = Animation()
    xnodes = [get_coords(node, 1) for node in get_nodes(mesh)]
    ynodes = [get_coords(node, 2) for node in get_nodes(mesh)]
    xcenters = [center[1] for center in Bcube.get_cell_centers(mesh)]
    ycenters = [center[2] for center in Bcube.get_cell_centers(mesh)]

    ## Initial solution
    lim_u, _u = linear_scaling_limiter(u, dΩ; DMPrelax, mass = M)
    isLimiterActive && (u.dofValues .= _u.dofValues)

    u_mean = cell_mean(u, dΩ)
    t = 0.0
    plt = plot_solution(0, t, u, mesh, xcenters, ycenters, xnodes, ynodes)
    append_vtk(vtk, u, lim_u, u_mean, t)
    frame(anim, plt)
    dofOverTime[1, :] .= t, u.dofValues[i_dof_out]

    b = Bcube.allocate_dofs(U)
    for ite in 1:nite
        b .= 0.0

        ## Apply limitation
        if isLimiterActive
            lim_u, _u = linear_scaling_limiter(u, dΩ; DMPrelax, mass = M)
            set_dof_values!(u, get_dof_values(_u))
        end

        ## Define linear forms
        flux = upwind ∘ (side⁻(u), side⁺(u), side⁻(c), side⁻(nΓ))
        l_Γ(v) = ∫(-flux * jump(v))dΓ
        l_Ω(v) = ∫(u * (c ⋅ ∇ₛ(v)))dΩ ## linear Volumic convective term
        l(v) = l_Ω(v) + l_Γ(v)

        if volumic_bilinear
            ## Version bilinear volumic term
            assemble_linear!(b, l_Γ, V)
            u.dofValues .= (I + Δt .* invMK) * u.dofValues + Δt .* invM * b
        else
            ## Version linear volumic term
            assemble_linear!(b, l, V)
            u.dofValues .+= Δt .* invM * b
        end

        t += Δt

        ## Output results
        if ite % (nite ÷ _nout) == 0
            u_mean = cell_mean(u, dΩ)
            append_vtk(vtk, u, lim_u, u_mean, t)
            plt = plot_solution(vtk.ite, t, u, mesh, xcenters, ycenters, xnodes, ynodes)
            frame(anim, plt)
        end
        dofOverTime[ite + 1, :] .= t, u.dofValues[i_dof_out]
    end

    ## Output final result and anim
    path = joinpath(out_dir, filename * ".csv")
    @info "Writing to $path"
    open(path, "w") do io
        println(io, "t,u")
        writedlm(io, dofOverTime, ",")
    end
    println("Computation is done, building gif...")
    g = gif(anim, joinpath(out_dir, "$filename.gif"); fps = 4)
    display(g)
end

"""
Linear transport of a vector quantity on a circle
"""
function vector_circle(; degree, nite, CFL, nθ)
    ## Settings
    radius = 1.0
    C = 1.0 ## velocity norm

    ## Mesh
    mesh = circle_mesh(nθ; radius = radius, order = 1)
    dΩ = Measure(CellDomain(mesh), 2 * degree + 1)
    Γ = InteriorFaceDomain(mesh)
    dΓ = Measure(Γ, 2 * degree + 1)
    nΓ = get_face_normals(Γ)

    ## Operators
    P = Bcube.tangential_projector(mesh)
    R = Bcube.CoplanarRotation()

    ## Transport velocity : it must be coplanar to each element, so we use the
    ## tangential projector operator and force the projected velocity to have
    ## the same norm as the "analytical" one
    _c = PhysicalFunction(x -> C * SA[-x[2], x[1]] / radius)
    c = (x -> C * normalize(x)) ∘ (P * _c) ## useless in theory since velocity is already tangent

    ## FESpace
    fs = FunctionSpace(:Lagrange, degree)
    U = TrialFESpace(fs, mesh, :discontinuous; size = 2)
    V = TestFESpace(U)

    ## FEFunction and "boundary / source" condition
    u = FEFunction(U)

    ## Initial condition
    _u0 = PhysicalFunction(x -> begin
        θ = atan(x[2], x[1])
        if 0 <= θ <= 2π / nθ
            return radius * SA[-x[2], x[1]]
        else
            return SA[0, 0]
        end
    end)
    projection_l2!(u, _u0, mesh)

    ## Forms
    m(u, v) = ∫(u ⋅ v)dΩ ## Mass matrix
    l_Ω(v) = ∫((u ⊗ c) ⊡ ∇ₛ(v))dΩ ## Volumic convective term

    function upwind(ui, uj, Ri, Rj, vi, vj, ci, nij)
        _uj = Ri * uj

        cij = ci ⋅ nij
        if cij > zero(cij)
            fi = cij * ui
        else
            fi = cij * _uj
        end

        return fi ⋅ (vi - transpose(Rj) * vj)
    end

    function flux(v)
        upwind ∘
        (side⁻(u), side⁺(u), side⁻(R), side⁺(R), side⁻(v), side⁺(v), side⁻(c), side⁻(nΓ))
    end
    l_Γ(v) = ∫(-flux(v))dΓ

    l(v) = l_Ω(v) + l_Γ(v)

    ## Time step
    dl = 2π * radius / nθ ## analytic length
    dl = 2 * radius * sin(2π / nθ / 2) ## discretized length
    Δt = CFL * dl / C
    @show dl
    @show Δt

    ## Mass
    M = assemble_bilinear(m, U, V)
    invM = inv(Matrix(M)) ##WARNING : really expensive !!!
    ## display(invM)

    ## Anim
    anim = Animation()
    xnodes = [get_coords(node, 1) for node in get_nodes(mesh)]
    ynodes = [get_coords(node, 2) for node in get_nodes(mesh)]
    xcenters = [center[1] for center in Bcube.get_cell_centers(mesh)]
    ycenters = [center[2] for center in Bcube.get_cell_centers(mesh)]

    ## Initial solution
    plt = plot_solution(0, 0.0, u, mesh, xcenters, ycenters, xnodes, ynodes)
    frame(anim, plt)

    t = 0.0
    b = Bcube.allocate_dofs(U)
    for i in 1:nite
        b .= 0.0
        assemble_linear!(b, l, V)

        u.dofValues .+= Δt .* invM * b

        t += Δt

        ## Build animation
        plt = plot_solution(i, t, u, mesh, xcenters, ycenters, xnodes, ynodes)
        frame(anim, plt)
    end

    g = gif(anim, joinpath(out_dir, "vector_on_circle_d$degree.gif"); fps = 2)
    display(g)
end

"""
Linear transport of a scalar quantity on a cylinder
"""
function scalar_cylinder(;
    degree,
    CFL,
    lz,
    nz,
    nθ,
    tmax,
    ϕ, ## velocity angle with respect to z axis
    C, ## velocity norm
    radius = 1,
    nout = 100,
    nitemax = Int(1e9),
    isLimiterActive = true,
    progressBar = true,
    meshOrder = 1,
    profile = false,
)
    function append_vtk(vtk, u::Bcube.AbstractFEFunction, lim_u, t)
        vars = Dict(
            "u" => u,
            "u_mean" => cell_mean(u, vtk.dΩ),
            "lim_u" => lim_u,
            "c" => vtk.c,
            "cellnormal" => Bcube.CellNormal(vtk.mesh),
            "u_warp" => u * Bcube.CellNormal(vtk.mesh),
        )
        Bcube.write_vtk_lagrange(
            vtk.basename * "_lag",
            vars,
            vtk.mesh,
            vtk.U,
            vtk.ite,
            t;
            collection_append = vtk.ite > 0,
        )

        ## Update counter
        vtk.ite += 1
    end

    ## Mesh
    mesh_path = joinpath(out_dir, "mesh.msh")
    Bcube.gen_cylinder_shell_mesh(
        mesh_path,
        nθ,
        nz;
        lz,
        radius,
        lc = 1e-1,
        recombine = true,
        transfinite = true,
        order = meshOrder,
    )
    mesh = read_mesh(mesh_path)
    rng = Random.MersenneTwister(33)
    θ = rand(rng, 3) .* 2π
    println("θx, θy, θz = $(rad2deg.(θ))")
    Rmat = rotMat(θ...)
    RmatInv = inv(Rmat)
    transform!(mesh, x -> Rmat * x)

    ## Domains
    ## quad = Quadrature(QuadratureLobatto(), 2 * degree + 1)
    quad = Quadrature(QuadratureLegendre(), 2 * degree + 1)
    dΩ = Measure(CellDomain(mesh), quad)
    Γ = InteriorFaceDomain(mesh)
    dΓ = Measure(Γ, quad)
    nΓ = get_face_normals(Γ)
    Γ_bnd = BoundaryFaceDomain(mesh, ("zmin", "zmax"))
    dΓ_bnd = Measure(Γ_bnd, quad)
    nΓ_bnd = get_face_normals(Γ_bnd)

    ## FESpace
    fs = FunctionSpace(:Lagrange, degree)
    U = TrialFESpace(fs, mesh, :discontinuous)
    V = TestFESpace(U)

    ## Transport velocity
    Cz = C * cos(ϕ)
    Cθ = C * sin(ϕ)
    _c = PhysicalFunction(x -> begin
        _x = RmatInv * x
        Rmat * SA[-Cθ * _x[2] / radius, Cθ * _x[1] / radius, Cz]
    end)
    P = Bcube.tangential_projector(mesh) ##Bcube.TangentialProjector()
    c = (x -> C * normalize(x)) ∘ (P * _c)

    ## Find quadrature weight (mesh is composed of a unique "shape" so first element is enough)
    quad = Bcube.get_quadrature(dΩ)
    s = Bcube.shape(Bcube.cells(mesh)[1])
    qrule = Bcube.QuadratureRule(s, quad)
    ω_quad = degree > 0 ? Bcube.get_weights(qrule)[1] : 1.0

    ## Time step and else
    dlθ = 2π * radius / nθ ## analytic length
    dlθ = 2 * radius * sin(2π / nθ / 2) ## discretized length
    dlz = lz / (nz - 1)
    println("Timestep constrained by $(dlθ < dlz ? "θ" : "z") discretization")
    dl = min(dlθ, dlz)
    Δt = CFL * dl * ω_quad / C / (2 * degree + 1)
    t = 0.0
    nite = min(floor(Int, tmax / Δt), nitemax)
    _nout = min(nite, nout)

    @show nite
    @show dl
    @show Δt
    @show get_ndofs(U)

    ## Limitation
    DMPrelax = 0.0 * dl

    ## Output
    tail = isLimiterActive ? "lim" : "nolim"
    filename = "scalar-on-cylinder-d$(degree)-$(tail)"
    U_export =
        TrialFESpace(FunctionSpace(:Lagrange, max(degree, meshOrder)), mesh, :discontinuous)
    vtk = VtkHandler(joinpath(out_dir, filename), dΩ, U_export, c)

    ## FEFunction and initial solution (P3 Gaussian bump)
    u = FEFunction(U)
    _θ0 = 0
    x0 = Rmat * SA[radius * cos(_θ0), radius * sin(_θ0), 0.2 * lz] ## bump center (in rotated frame)
    _r = 1 ## bump radius
    _umax = 1 ## bump amplitude
    _a, _b = SA[
        _r^3 _r^2
        3*_r^2 2*_r
    ] \ SA[-_umax; 0]
    f = PhysicalFunction(x -> begin
        dx = norm(x - x0)
        dx < _r ? _a * dx^3 + _b * dx^2 + _umax : 0.0
    end)

    projection_l2!(u, f, mesh)

    ## Forms
    m(u, v) = ∫(u ⋅ v)dΩ ## Mass matrix

    function upwind(ui, uj, ci, nij)
        cij = ci ⋅ nij
        if cij > zero(cij)
            fij = cij * ui
        else
            fij = cij * uj
        end
        fij
    end

    ## Mass
    M = factorize(assemble_bilinear(m, U, V))

    ## Initial solution
    lim_u, _u = linear_scaling_limiter(u, dΩ; DMPrelax, mass = M)
    isLimiterActive && (u.dofValues .= _u.dofValues)

    t = 0.0
    append_vtk(vtk, u, lim_u, t)

    b = Bcube.allocate_dofs(U)
    du = similar(b)
    progressBar && (progress = Progress(nitemax))
    for ite in 1:nitemax
        b .= 0.0

        ## Apply limitation
        if isLimiterActive
            lim_u, _u = linear_scaling_limiter(u, dΩ; DMPrelax, mass = M)
            set_dof_values!(u, get_dof_values(_u))
        end

        ## Define linear forms
        flux = upwind ∘ (side⁻(u), side⁺(u), side⁻(c), side⁻(nΓ))
        l_Γ(v) = ∫(-flux * jump(v))dΓ
        flux_bnd = upwind ∘ (side⁻(u), side⁺(u), side⁻(c), side⁻(nΓ_bnd))
        l_Γ_bnd(v) = ∫(-flux_bnd * jump(v))dΓ_bnd
        l_Ω(v) = ∫(u * (c ⋅ ∇ₛ(v)))dΩ ## linear Volumic convective term
        l(v) = l_Ω(v) + l_Γ(v) + l_Γ_bnd(v)

        ## Version linear volumic term
        assemble_linear!(b, l, V)
        du .= M \ b
        @. u.dofValues += Δt * du

        t += Δt
        progressBar && next!(progress)

        ## Output results
        if ite % (nitemax ÷ _nout) == 0
            append_vtk(vtk, u, lim_u, t)
        end

        if ite == nitemax && profile
            println("ndofs total = ", Bcube.get_ndofs(U))
            Profile.init(; n = 10^7) ## returns the current settings
            Profile.clear()
            Profile.clear_malloc_data()
            @profile begin
                for i in 1:1000
                    assemble_linear!(b, l, V)
                end
            end
            @btime assemble_linear!($b, $l, $V)
        end
    end
end

"""
Linear transport of a vector quantity on a cylinder
"""
function vector_cylinder(;
    degree,
    CFL,
    lz,
    nz,
    nθ,
    tmax,
    ϕ_vel, ## velocity angle with respect to z axis
    C_vel, ## velocity norm
    ϕ_u, ## transported vector angle with respect to z axis
    C_u, ## transported vector norm
    radius = 1,
    nout = 100,
    nitemax = Int(1e9),
    isLimiterActive = true,
    progressBar = true,
)
    function append_vtk(vtk, u::Bcube.AbstractFEFunction, lim_u, u_mean, t)
        ## Build animation
        values_vertices = transpose(var_on_vertices(u, mesh))
        values_centers = transpose(var_on_centers(u, mesh))
        ## values_nodes = var_on_nodes_discontinuous(u, mesh, degree)

        ## Write
        Bcube.write_vtk(
            vtk.basename,
            vtk.ite,
            t,
            vtk.mesh,
            Dict(
                "u_centers" => (values_centers, VTKCellData()),
                "u_vertices" => (values_vertices, VTKPointData()),
                "lim_u" => (get_values(lim_u), VTKCellData()),
                "u_mean" => (get_values(u_mean), VTKCellData()),
                "θ_centers" => (vtk.θ_centers, VTKCellData()),
                "θ_vertices" => (vtk.θ_vertices, VTKPointData()),
                "c_centers" => (vtk.c_centers, VTKCellData()),
                "c_vertices" => (vtk.c_vertices, VTKPointData()),
                "ν_centers" => (vtk.ν_centers, VTKCellData()),
                "ν_vertices" => (vtk.ν_vertices, VTKPointData()),
            ),
            ;
            append = vtk.ite > 0,
        )

        ## Update counter
        vtk.ite += 1
    end

    ## Mesh
    mesh_path = joinpath(out_dir, "mesh.msh")
    Bcube.gen_cylinder_shell_mesh(
        mesh_path,
        nθ,
        nz;
        lz,
        radius,
        lc = 1e-1,
        recombine = true,
        transfinite = true,
    )
    mesh = read_mesh(mesh_path)
    rng = Random.MersenneTwister(33)
    θ = rand(rng, 3) .* 2π
    println("θx, θy, θz = $(rad2deg.(θ))")
    Rmat = rotMat(θ...)
    RmatInv = inv(Rmat)
    transform!(mesh, x -> Rmat * x)

    ## Domains
    ## quad = Quadrature(QuadratureLobatto(), 2 * degree + 1)
    quad = Quadrature(QuadratureLegendre(), 2 * degree + 1)
    dΩ = Measure(CellDomain(mesh), quad)
    Γ = InteriorFaceDomain(mesh)
    dΓ = Measure(Γ, quad)
    nΓ = get_face_normals(Γ)
    Γ_bnd = BoundaryFaceDomain(mesh, ("zmin", "zmax"))
    dΓ_bnd = Measure(Γ_bnd, quad)
    nΓ_bnd = get_face_normals(Γ_bnd)

    ## FESpace
    fs = FunctionSpace(:Lagrange, degree)
    U = TrialFESpace(fs, mesh, :discontinuous; size = Bcube.spacedim(mesh))
    V = TestFESpace(U)

    ## Operators
    P = Bcube.tangential_projector(mesh)
    R = Bcube.CoplanarRotation()

    ## Transport velocity
    Cz = C_vel * cos(ϕ_vel)
    Cθ = C_vel * sin(ϕ_vel)
    _c = PhysicalFunction(x -> begin
        _x = RmatInv * x
        Rmat * SA[-Cθ * _x[2] / radius, Cθ * _x[1] / radius, Cz]
    end)
    c = (x -> C_vel * normalize(x)) ∘ (P * _c) ## useless in theory because velocity is already tangent

    ## Find quadrature weight (mesh is composed of a unique "shape" so first element is enough)
    quad = Bcube.get_quadrature(dΩ)
    s = Bcube.shape(Bcube.cells(mesh)[1])
    qrule = Bcube.QuadratureRule(s, quad)
    ω_quad = degree > 0 ? Bcube.get_weights(qrule)[1] : 1.0

    ## Time step and else
    dlθ = 2π * radius / nθ ## analytic length
    dlθ = 2 * radius * sin(2π / nθ / 2) ## discretized length
    dlz = lz / (nz - 1)
    println("Timestep constrained by $(dlθ < dlz ? "θ" : "z") discretization")
    dl = min(dlθ, dlz)
    Δt = CFL * dl * ω_quad / C_vel / (2 * degree + 1)
    t = 0.0
    nite = min(floor(Int, tmax / Δt), nitemax)
    _nout = min(nite, nout)

    @show nite
    @show dl
    @show Δt
    @show get_ndofs(U)

    ## Limitation
    DMPrelax = 0.0 * dl

    ## Output
    tail = isLimiterActive ? "lim" : "nolim"
    filename = "vector-on-cylinder-d$(degree)-$(tail)"
    vtk = VtkHandler(joinpath(out_dir, filename), dΩ, U, c)

    ## FEFunction and initial solution
    u = FEFunction(U)
    _θ0 = 0
    x0 = Rmat * [radius * cos(_θ0), radius * sin(_θ0), 0.2 * lz] ## bump center (in rotated frame)
    _r = 1 ## bump radius
    _umax = C_u ## bump amplitude
    _a, _b = [
        _r^3 _r^2
        3*_r^2 2*_r
    ] \ [-_umax; 0]
    function norm_bump_p3(x)
        dx = norm(x - x0)
        y = dx < _r ? _a * dx^3 + _b * dx^2 + _umax : 0.0
        return y * cos(ϕ_u), y * sin(ϕ_u)
    end
    _f = PhysicalFunction(
        x -> begin
            _Cz, _Cθ = norm_bump_p3(x)
            _x = RmatInv * x
            return Rmat * SA[-_Cθ * _x[2] / radius, _Cθ * _x[1] / radius, _Cz]
        end,
    )
    f = C_u * (P * _f) / (mynorm(P * _f) + eps())
    projection_l2!(u, f, mesh)

    ## Forms
    m(u, v) = ∫(u ⋅ v)dΩ ## Mass matrix

    function upwind(ui, uj, Ri, Rj, vi, vj, ci, nij)
        _uj = Ri * uj

        cij = ci ⋅ nij
        if cij > zero(cij)
            fi = cij * ui
        else
            fi = cij * _uj
        end

        return fi ⋅ (vi - transpose(Rj) * vj)
    end

    function flux(v, n)
        upwind ∘
        (side⁻(u), side⁺(u), side⁻(R), side⁺(R), side⁻(v), side⁺(v), side⁻(c), side⁻(n))
    end

    ## Mass
    M = assemble_bilinear(m, U, V)

    ## Initial solution
    if isLimiterActive
        lim_u, _u = linear_scaling_limiter(u, dΩ; DMPrelax, mass = M)
        u.dofValues .= _u.dofValues
    else
        lim_u = MeshCellData(zero(get_dof_values(u))) ## dummy, just for the output
    end

    u_mean = cell_mean(u, dΩ)
    t = 0.0
    append_vtk(vtk, u, lim_u, u_mean, t)

    b = Bcube.allocate_dofs(U)
    du = similar(b)
    progressBar && (progress = Progress(nitemax))
    for ite in 1:nitemax
        b .= 0.0

        ## Apply limitation
        if isLimiterActive
            lim_u, _u = linear_scaling_limiter(u, dΩ; DMPrelax, mass = M)
            set_dof_values!(u, get_dof_values(_u))
        end

        ## Define linear forms
        l_Γ(v) = ∫(-flux(v, nΓ))dΓ
        l_Γ_bnd(v) = ∫(-flux(v, nΓ_bnd))dΓ_bnd
        l_Ω(v) = ∫((u ⊗ c) ⊡ ∇ₛ(v))dΩ ## linear Volumic convective term
        l(v) = l_Ω(v) + l_Γ(v) + l_Γ_bnd(v)

        ## Version linear volumic term
        assemble_linear!(b, l, V)
        du .= M \ b
        @. u.dofValues += Δt * du

        t += Δt
        progressBar && next!(progress)

        ## Output results
        if ite % (nitemax ÷ _nout) == 0
            u_mean = cell_mean(u, dΩ)
            append_vtk(vtk, u, lim_u, u_mean, t)
        end
    end
end

"""
Linear transport of a scalar quantity on a torus
"""
function scalar_torus(;
    degree,
    CFL,
    rint,
    rext,
    lc,
    tmax,
    ϕ, ## velocity angle with respect to z axis
    C, ## velocity norm
    nout = 100,
    nitemax = Int(1e9),
    isLimiterActive = true,
    progressBar = true,
    meshOrder = 1,
)
    function append_vtk(vtk, u::Bcube.AbstractFEFunction, lim_u, t)
        vars = Dict(
            "u" => u,
            "u_mean" => cell_mean(u, vtk.dΩ),
            "lim_u" => lim_u,
            "c" => vtk.c,
            "cellnormal" => Bcube.CellNormal(vtk.mesh),
            "u_warp" => u * Bcube.CellNormal(vtk.mesh),
        )
        Bcube.write_vtk_lagrange(
            vtk.basename * "_lag",
            vars,
            vtk.mesh,
            vtk.U,
            vtk.ite,
            t;
            collection_append = vtk.ite > 0,
        )

        ## Update counter
        vtk.ite += 1
    end

    ## Mesh
    mesh_path = joinpath(out_dir, "mesh.msh")
    Bcube.gen_torus_shell_mesh(
        mesh_path,
        rint,
        rext;
        lc,
        order = meshOrder,
        verbose = false,
    )
    mesh = read_mesh(mesh_path)
    rng = Random.MersenneTwister(33)
    θ = zeros(3)
    ## θ = rand(rng, 3) .* 2π
    println("θx, θy, θz = $(rad2deg.(θ))")
    Rmat = rotMat(θ...)
    RmatInv = inv(Rmat)
    transform!(mesh, x -> Rmat * x)

    ## Domains
    ## quad = Quadrature(QuadratureLobatto(), 2 * degree + 1)
    quad = Quadrature(QuadratureLegendre(), 2 * degree + 1)
    dΩ = Measure(CellDomain(mesh), quad)
    Γ = InteriorFaceDomain(mesh)
    dΓ = Measure(Γ, quad)
    nΓ = get_face_normals(Γ)

    ## FESpace
    fs = FunctionSpace(:Lagrange, degree)
    U = TrialFESpace(fs, mesh, :discontinuous)
    V = TestFESpace(U)

    ## Transport velocity
    Cθ = C * cos(ϕ)
    Cφ = C * sin(ϕ)
    rc = (rint + rext) / 2
    r = (rext - rint) / 2
    ez = SA[0, 0, 1]
    _c = PhysicalFunction(coords -> begin
        _x = RmatInv * coords
        x, y, z = _x

        ## In the (ex, ey) plane
        r_xy = √(x * x + y * y)
        cosθ = x / r_xy
        sinθ = y / r_xy
        er = SA[cosθ, sinθ, 0]
        eθ = SA[-sinθ, cosθ, 0]

        ## In the (er, ez) plane
        l = _x ⋅ er - rc
        r_rz = √(z * z + l * l)
        cosφ = l / r_rz
        sinφ = z / r_rz
        eφ = -sinφ * er + cosφ * ez

        ## direction vector
        v = Cθ * eθ + Cφ * eφ

        ## Rotate back
        Rmat * v
    end)
    ##P = Bcube.tangential_projector(mesh)
    ## c = (x -> C * normalize(x)) ∘ (P * _c) ## use this if `_c` is not necessarily tangent
    c = _c ## `_c` is anatically tangent, so no need to project

    ## Find quadrature weight (mesh is composed of a unique "shape" so first element is enough)
    quad = Bcube.get_quadrature(dΩ)
    s = Bcube.shape(Bcube.cells(mesh)[1])
    qrule = Bcube.QuadratureRule(s, quad)
    ω_quad = degree > 0 ? Bcube.get_weights(qrule)[1] : 1.0

    ## Time step and else
    dl = lc
    Δt = CFL * dl * ω_quad / C / (2 * degree + 1)
    t = 0.0
    nite = min(floor(Int, tmax / Δt), nitemax)
    _nout = min(nite, nout)

    @show nite
    @show dl
    @show Δt
    @show get_ndofs(U)

    ## Limitation
    DMPrelax = 0.0 * dl

    ## Output
    tail = isLimiterActive ? "lim" : "nolim"
    filename = "scalar-on-torus-d$(degree)-$(tail)"
    U_export =
        TrialFESpace(FunctionSpace(:Lagrange, max(degree, meshOrder)), mesh, :discontinuous)
    vtk = VtkHandler(joinpath(out_dir, filename), dΩ, U_export, c)

    ## FEFunction and initial solution (P3 Gaussian bump)
    u = FEFunction(U)
    _θ0 = π / 2
    x0 = Rmat * SA[rc + r * cos(_θ0), 0.0, r * sin(_θ0)] ## bump center (in rotated frame)
    _r = 0.5 ## bump radius
    _umax = 1 ## bump amplitude
    _a, _b = SA[
        _r^3 _r^2
        3*_r^2 2*_r
    ] \ SA[-_umax; 0]
    f = PhysicalFunction(x -> begin
        dx = norm(x - x0)
        dx < _r ? _a * dx^3 + _b * dx^2 + _umax : 0.0
    end)

    projection_l2!(u, f, mesh)

    ## Forms
    m(u, v) = ∫(u ⋅ v)dΩ ## Mass matrix

    function upwind(ui, uj, ci, nij)
        cij = ci ⋅ nij
        if cij > zero(cij)
            fij = cij * ui
        else
            fij = cij * uj
        end
        fij
    end

    ## Mass
    M = factorize(assemble_bilinear(m, U, V))

    ## Initial solution
    lim_u, _u = linear_scaling_limiter(u, dΩ; DMPrelax, mass = M)
    isLimiterActive && (u.dofValues .= _u.dofValues)

    t = 0.0
    append_vtk(vtk, u, lim_u, t)

    b = Bcube.allocate_dofs(U)
    du = similar(b)
    progressBar && (progress = Progress(nite))
    for ite in 1:nite
        b .= 0.0

        ## Apply limitation
        if isLimiterActive
            lim_u, _u = linear_scaling_limiter(u, dΩ; DMPrelax, mass = M)
            set_dof_values!(u, get_dof_values(_u))
        end

        ## Define linear forms
        flux = upwind ∘ (side⁻(u), side⁺(u), side⁻(c), side⁻(nΓ))
        l_Γ(v) = ∫(-flux * jump(v))dΓ
        l_Ω(v) = ∫(u * (c ⋅ ∇ₛ(v)))dΩ ## linear Volumic convective term
        l(v) = l_Ω(v) + l_Γ(v)

        ## Version linear volumic term
        assemble_linear!(b, l, V)
        du .= M \ b
        @. u.dofValues += Δt * du

        t += Δt
        progressBar && next!(progress)

        ## Output results
        if ite % (nite ÷ _nout) == 0
            append_vtk(vtk, u, lim_u, t)
        end
    end
end

## Run
scalar_circle(; degree = 1, nrot = 5, CFL = 0.1, nθ = 25, isLimiterActive = false)
vector_circle(; degree = 0, nite = 100, CFL = 1, nθ = 20)
@time scalar_cylinder(;
    degree = 1,
    CFL = 0.1,
    lz = 10,
    nθ = 50,
    nz = 70,
    ϕ = 0.5 * π / 2,
    C = 1.0,
    tmax = 10.0,
    nout = 100,
    nitemax = 2000,#Int(1e9),
    isLimiterActive = false,
    progressBar = true,
    meshOrder = 2,
)
@time vector_cylinder(;
    degree = 0,
    CFL = 0.1,
    lz = 10,
    nθ = 50,
    nz = 70,
    ϕ_vel = 0.5 * π / 2,
    C_vel = 1.0,
    ϕ_u = -0.5 * π / 2,
    C_u = 1.0,
    tmax = 10.0,
    nout = 100,
    nitemax = 50,#Int(1e9),
    isLimiterActive = false,
    progressBar = true,
)
@time scalar_torus(;
    degree = 1,
    CFL = 0.2, # d=0, CFL=0.4 OK
    rint = 1.0,
    rext = 1.5,
    lc = 0.08,
    ϕ = 0.5 * π / 2,
    C = 1.0,
    tmax = 11.0, # 11s  = 1 turn for rint=1, rext=1.5
    nout = 100,
    nitemax = 100000, # d = 4, n = 400 OK,
    isLimiterActive = false,
    progressBar = true,
    meshOrder = 2,
)

end #hide
