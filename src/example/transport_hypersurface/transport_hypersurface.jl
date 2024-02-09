module TransportHypersurface

using Plots
using Bcube
using StaticArrays
using LinearAlgebra

# Common settings
const degree = 0
const nite = 25
const CFL = 1

const out_dir = joinpath(@__DIR__, "../../../myout/transport_hypersurface")
mkpath(out_dir)

""" Hack waiting for Ghislain to finish his branch """
mynorm(a) = sqrt(a ⋅ a)

function scalar_circle()
    function plot_solution(u, mesh, degree, xplot, yplot, xnodes, ynodes)
        # Build animation
        if degree > 0
            uplot = var_on_vertices(u, mesh)
        else
            uplot = var_on_centers(u, mesh)
        end
        plt =
            plt = plot(
                [xnodes..., xnodes[1]],
                [ynodes..., ynodes[1]];
                aspect_ratio = :equal,
                primary = false,
            )
        scatter!(xplot, yplot; marker_z = uplot, label = "u")

        return plt
    end

    # Settings
    nθ = 10
    radius = 1.0
    C = 1.0 # velocity norm

    # Mesh
    mesh = circle_mesh(nθ; radius = radius, order = 1)
    dΩ = Measure(CellDomain(mesh), 2 * degree + 1)
    Γ = InteriorFaceDomain(mesh)
    dΓ = Measure(Γ, 2 * degree + 1)
    nΓ = get_face_normals(Γ)

    # Transport velocity
    _c = PhysicalFunction(x -> C * SA[-x[2], x[1]] / radius)
    P = Bcube.TangentialProjector()
    c = C * (P * _c) / mynorm(P * _c)

    # FESpace
    fs = FunctionSpace(:Lagrange, degree)
    U = TrialFESpace(fs, mesh, :discontinuous)
    V = TestFESpace(U)

    # FEFunction and "boundary / source" condition
    u = FEFunction(U)
    u.dofValues[1] = 1.0

    # Forms
    m(u, v) = ∫(u ⋅ v)dΩ # Mass matrix
    l_Ω(v) = ∫(u * (c ⋅ ∇(v)))dΩ # Volumic convective term

    function upwind(ui, uj, ci, nij)
        cij = ci ⋅ nij
        if cij > zero(cij)
            fij = cij * ui
        else
            fij = cij * uj
        end
        fij
    end
    flux = upwind ∘ (side⁻(u), side⁺(u), side⁻(c), side⁻(nΓ))
    l_Γ(v) = ∫(-flux * jump(v))dΓ

    l(v) = l_Ω(v) + l_Γ(v)

    # Time step
    dl = 2π * radius / nθ # analytic length
    dl = 2 * radius * sin(2π / nθ / 2) # discretized length
    Δt = CFL * dl / C
    @show dl
    @show Δt

    # Mass
    M = assemble_bilinear(m, U, V)
    invM = inv(Matrix(M)) #WARNING : really expensive !!!
    # display(invM)
    # display(Δt .* invM)

    # Anim
    anim = Animation()
    xnodes = [Bcube.coords(node, 1) for node in get_nodes(mesh)]
    ynodes = [Bcube.coords(node, 2) for node in get_nodes(mesh)]
    if degree > 0
        xplot = xnodes
        yplot = ynodes
    else
        xplot = [center[1] for center in Bcube.get_cell_centers(mesh)]
        yplot = [center[2] for center in Bcube.get_cell_centers(mesh)]
    end

    # Initial solution
    plt = plot_solution(u, mesh, degree, xplot, yplot, xnodes, ynodes)
    frame(anim, plt)

    b = Bcube.allocate_dofs(U)
    for i in 1:nite
        b .= 0.0
        assemble_linear!(b, l, V)
        u.dofValues .+= Δt .* invM * b

        # Build animation
        plt = plot_solution(u, mesh, degree, xplot, yplot, xnodes, ynodes)
        frame(anim, plt)
    end

    g = gif(anim, joinpath(out_dir, "scalar_on_circle.gif"); fps = 2)
    display(g)
end

function vector_circle()
    function plot_solution(u, mesh, degree, xplot, yplot, xnodes, ynodes)
        # Build animation
        if degree > 0
            uplot = var_on_vertices(u, mesh)
        else
            uplot = var_on_centers(u, mesh)
        end
        L = maximum(x -> norm(x), eachrow(uplot))
        println("Max norm of u : $L")
        plt =
            plt = plot(
                [xnodes..., xnodes[1]],
                [ynodes..., ynodes[1]];
                aspect_ratio = :equal,
                primary = false,
            )
        quiver!(xplot, yplot; quiver = (uplot[:, 1], uplot[:, 2]), label = "u")
        xlims!(-2, 2)
        ylims!(-2, 2)

        return plt
    end

    # Settings
    nθ = 10
    radius = 1.0
    C = 1.0 # velocity norm

    # Mesh
    mesh = circle_mesh(nθ; radius = radius, order = 1)
    dΩ = Measure(CellDomain(mesh), 2 * degree + 1)
    Γ = InteriorFaceDomain(mesh)
    dΓ = Measure(Γ, 2 * degree + 1)
    nΓ = get_face_normals(Γ)

    # Transport velocity : it must be coplanar to each element, so we use the
    # tangential projector operator and force the projected velocity to have
    # the same norm as the "analytical" one
    _c = PhysicalFunction(x -> C * SA[-x[2], x[1]] / radius)
    P = Bcube.TangentialProjector()
    c = -C * (P * _c) / mynorm(P * _c)

    # Operators
    P = Bcube.TangentialProjector()
    R = Bcube.CoplanarRotation()

    # FESpace
    fs = FunctionSpace(:Lagrange, degree)
    U = TrialFESpace(fs, mesh, :discontinuous; size = 2)
    V = TestFESpace(U)

    # FEFunction and "boundary / source" condition
    u = FEFunction(U)

    # Initial condition
    _u0 = PhysicalFunction(x -> begin
        θ = atan(x[2], x[1])
        if 0 <= θ <= 2π / nθ
            return radius * SA[-x[2], x[1]]
        else
            return SA[0, 0]
        end
    end)
    projection_l2!(u, _u0, mesh)

    # Forms
    m(u, v) = ∫(u ⋅ v)dΩ # Mass matrix
    l_Ω(v) = ∫((u ⊗ c) ⊡ ∇(v))dΩ # Volumic convective term

    function upwind(ui, uj, Ri, Rj, vi, vj, ci, nij)
        _uj = Ri * uj

        cij = ci ⋅ nij
        if cij > zero(cij)
            fi = cij * ui
        else
            fi = cij * _uj
        end

        fj = Rj * fi

        return fi ⋅ vi - fj ⋅ vj
    end

    function flux(v)
        upwind ∘
        (side⁻(u), side⁺(u), side⁻(R), side⁺(R), side⁻(v), side⁺(v), side⁻(c), side⁻(nΓ))
    end
    l_Γ(v) = ∫(-flux(v))dΓ

    l(v) = l_Ω(v) + l_Γ(v)

    # Time step
    dl = 2π * radius / nθ # analytic length
    dl = 2 * radius * sin(2π / nθ / 2) # discretized length
    Δt = CFL * dl / C
    @show dl
    @show Δt

    # Mass
    M = assemble_bilinear(m, U, V)
    invM = inv(Matrix(M)) #WARNING : really expensive !!!
    # display(invM)

    # Anim
    anim = Animation()
    xnodes = [Bcube.coords(node, 1) for node in get_nodes(mesh)]
    ynodes = [Bcube.coords(node, 2) for node in get_nodes(mesh)]
    if degree > 0
        xplot = xnodes
        yplot = ynodes
    else
        xplot = [center[1] for center in Bcube.get_cell_centers(mesh)]
        yplot = [center[2] for center in Bcube.get_cell_centers(mesh)]
    end

    # Initial solution
    plt = plot_solution(u, mesh, degree, xplot, yplot, xnodes, ynodes)
    frame(anim, plt)

    b = Bcube.allocate_dofs(U)
    for i in 1:nite
        b .= 0.0
        assemble_linear!(b, l, V)

        u.dofValues .+= Δt .* invM * b

        # Build animation
        plt = plot_solution(u, mesh, degree, xplot, yplot, xnodes, ynodes)
        frame(anim, plt)
    end

    g = gif(anim, joinpath(out_dir, "vector_on_circle.gif"); fps = 2)
    display(g)
end

scalar_circle()
vector_circle()

end