module TransportHypersurface

using Plots
using Bcube
using StaticArrays
using LinearAlgebra
using Printf

# Common settings
const degree = 1
const nite = 50
const CFL = 0.05

const out_dir = joinpath(@__DIR__, "../../../myout/transport_hypersurface")
mkpath(out_dir)

""" Hack waiting for Ghislain to finish his branch """
mynorm(a) = sqrt(a ⋅ a)

function scalar_circle()
    function plot_solution(i, t, u, mesh, degree, xplot, yplot, xnodes, ynodes)
        # Build animation
        if degree > 0
            uplot = var_on_vertices(u, mesh)
        else
            uplot = var_on_centers(u, mesh)
        end
        lmax = 1.2
        plt =
            plt = plot(
                [xnodes..., xnodes[1]],
                [ynodes..., ynodes[1]];
                aspect_ratio = :equal,
                primary = false,
                xlim = (-lmax, lmax),
                ylim = (-lmax, lmax),
            )
        scatter!(xplot, yplot; marker_z = uplot, label = "u", clims = (-1, 1))
        annotate!(0, 0.25, "i  = $i")
        annotate!(0, 0, @sprintf "t = %.2e" t)
        annotate!(0, -0.25, @sprintf "|u|_max = %.2e" maximum(abs.(uplot)))

        return plt
    end

    # Settings
    nθ = 10
    radius = 1.0
    C = 1.0 # velocity norm

    # Mesh
    qOrder = 2 * degree + 1
    mesh = circle_mesh(nθ; radius = radius, order = 1)
    dΩ = Measure(CellDomain(mesh), qOrder)
    Γ = InteriorFaceDomain(mesh)
    dΓ = Measure(Γ, qOrder)
    nΓ = get_face_normals(Γ)

    # Transport velocity
    _c = PhysicalFunction(x -> C * SA[-x[2], x[1]] / radius)
    P = Bcube.TangentialProjector()
    c = C * (P * _c) / mynorm(P * _c)
    # c = PhysicalFunction(x -> C * SA[-x[2], x[1]] / radius)
    # println("DEBUG !!!!")

    # FESpace
    fs = FunctionSpace(:Lagrange, degree)
    U = TrialFESpace(fs, mesh, :discontinuous)
    V = TestFESpace(U)

    # FEFunction and "boundary / source" condition
    u = FEFunction(U)
    if false
        u.dofValues[1] = 1.0
    else
        projection_l2!(u, PhysicalFunction(x -> cos(atan(x[2], x[1]))), mesh)
    end

    # Forms
    m(u, v) = ∫(u ⋅ v)dΩ # Mass matrix
    a_Ω(u, v) = ∫(u * (c ⋅ ∇ₛ(v)))dΩ # bilinear volumic convective term
    l_Ω(v) = ∫(u * (c ⋅ ∇ₛ(v)))dΩ # linear Volumic convective term

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
    Δt = CFL * dl / C / (2 * degree + 1)
    @show dl
    @show Δt

    # Mass
    M = assemble_bilinear(m, U, V)
    K = assemble_bilinear(a_Ω, U, V)
    invM = inv(Matrix(M)) #WARNING : really expensive !!!
    invMK = invM * K
    # display(K)
    # display(invMK)
    # display(invM)
    # display(Δt .* invM)

    # Anim
    anim = Animation()
    xnodes = [Bcube.coords(node, 1) for node in Bcube.get_nodes(mesh)]
    ynodes = [Bcube.coords(node, 2) for node in Bcube.get_nodes(mesh)]
    if degree > 0
        xplot = xnodes
        yplot = ynodes
    else
        xplot = [center[1] for center in Bcube.get_cell_centers(mesh)]
        yplot = [center[2] for center in Bcube.get_cell_centers(mesh)]
    end

    # Initial solution
    t = 0.0
    plt = plot_solution(0, t, u, mesh, degree, xplot, yplot, xnodes, ynodes)
    frame(anim, plt)

    b = Bcube.allocate_dofs(U)
    for i in 1:nite
        b .= 0.0

        # Version linear volumic term
        # assemble_linear!(b, l, V)
        # u.dofValues .+= Δt .* invM * b

        # Version bilinear volumic term
        assemble_linear!(b, l_Γ, V)
        u.dofValues .= (I + Δt .* invMK) * u.dofValues + Δt .* invM * b

        t += Δt

        # Build animation
        plt = plot_solution(i, t, u, mesh, degree, xplot, yplot, xnodes, ynodes)
        frame(anim, plt)
    end

    g = gif(anim, joinpath(out_dir, "scalar_on_circle.gif"); fps = 4)
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
    l_Ω(v) = ∫((u ⊗ c) ⊡ ∇ₛ(v))dΩ # Volumic convective term

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
# vector_circle()

end