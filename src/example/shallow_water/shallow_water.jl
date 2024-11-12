module ShallowWater #hide
println("Running shallow_water example...") #hide
# # Shallow water equations
#
# WARNING : the below explanations are not up to date, the text has been copied from the old-api example.
#
# Following "A conservative Saint-Venant type model to describe the dynamics of thien partially wetting films with
# regularized forces at the contact line".
# The gravity is noted ``g = g_n \vec{e_n} + \vec{g_t}`` (note that ``g_n`` is a scalar while ``g_t`` is a vector). The goal is to solve:
# ```math
# \begin{cases}
#   \partial_t \rho h + \nabla \cdot \rho h u = 0 \\
#   \partial_t \rho h u + \nabla \cdot \mathcal{F}_{\rho h u} = h \left( \rho g_t - \nabla P_{gaz} \right) - \tilde{\tau}_{wall} + \tilde{\tau}_{air}
# \end{cases}
# ```
#
# To simplify a little bit, we assume a constant density. The systems becomes:
# ```math
# \begin{cases}
#   \partial_t h + \nabla \cdot h u = 0 \\
#   \partial_t h u + \nabla \cdot \mathcal{F}_{hu} = h \left( g_t - \nabla P_{gaz} \right) - \tau_{wall} + \tau_{air}
# \end{cases}
# ```
# where
# ```math
# \tau_{wall} = \frac{3 \nu u}{h + b} - \frac{\tau_{air}h}{2(h+b)}
# ```
# ``b`` being a slip length and
# ```math
# \mathcal{F}_{h u} = h u \otimes u + g_n \frac{h^2}{2} \mathcal{I} + \frac{1}{\rho}\left[h \partial_h e_d(h) - e_d(h) \right]
# - \frac{\gamma_{lg}}{\rho \sqrt{1 + ||\nabla h||^2}} + \frac{1}{\rho} \gamma_{lg} h \kappa
# ```
#
# ## Explicit time integration
# ```math
# \begin{cases}
#   h^{n+1} = h^n - \Delta t \nabla \cdot h u^n \\
#   h u^{n+1} =  hu^n - \Delta t \left[
#     \nabla \cdot \mathcal{F}_{hu}(h^n,hu^n)
#     - h^n \left( g_t - \nabla P_{gaz} \right) + \tau_{wall}(h^n, hu^n) - \tau_{air}
#   \right]
# \end{cases}
# ```
#
# ## Implicit time integration
# ```math
# \begin{cases}
#   h^{n+1} = h^n - \Delta t \nabla \cdot h u^{n+1} \\
#   h u^{n+1} =  hu^n - \Delta t \left[
#     \nabla \cdot \mathcal{F}_{hu}(h^{n+1},hu^{n+1})
#     - h^{n+1} \left( g_t - \nabla P_{gaz} \right) + \tau_{wall}(h^{n+1}, hu^{n+1}) - \tau_{air}
#   \right]
# \end{cases}
# ```
#
# ## IMEX time integration (not finished / not relevant)
# The wall friction term, ``\tau_w`` is singular when ``h \rightarrow 0``. To overcome this difficulty, an implicit-explicit (IMEX)
# scheme is used : all terms are integrated explicitely except the wall friction. More precisely, the system is first written:
# ```math
# \begin{cases}
#   h^{n+1} = h^n - \Delta t \nabla \cdot h u^n \\
#   h u^{n+1} =  hu^n - \Delta t \left[
#     \nabla \cdot \mathcal{F}_{hu}(h^n,hu^n)
#     - h^n \left( g_t - \nabla P_{gaz} \right) + \tau_{wall}(h^{n+1}, hu^{n+1}) - \tau_{air}
#   \right]
# \end{cases}
# ```
# At each time step, the mass equation can be solved explicitely independantly from the momentum equation. Besides, the wall
# friction can be expressed as:
# ```math
# \tau_{wall} = \frac{3 \nu hu^{n+1}}{{h^{n+1}}^2} - \frac{\tau_{air}}{2}
# ```
# where the slipping length, ``b``, is neglected (which is fine when working with an implicit formulation). The momentum
# equation can then be rearranged to obtain:
# ```math
#   h u^{n+1}\left( 1 + \frac{3  \nu \Delta t }{{h^{n+1}}^2} \right) =  hu^n - \Delta t \left[
#     \nabla \cdot \mathcal{F}_{hu}(h^n,hu^n)
#     - h^n \left( g_t - \nabla P_{gaz} \right) - \frac{3}{2}\tau_{air}
#   \right]
# ```
# Moving the multiplying factor to the right-hand-side, we finally obtain:
# ```math
#   h u^{n+1} = \frac{{h^{n+1}}^2}{{h^{n+1}}^2 + 3 \nu \Delta t} \left[
#     hu^n - \Delta t \left[
#       \nabla \cdot \mathcal{F}_{hu}(h^n,hu^n)
#       - h^n \left( g_t - \nabla P_{gaz} \right) - \frac{3}{2}\tau_{air}
#     \right]
#  \right]
# ```
#
#
# ## Weak form
# First we write the different equation with a full explicit scheme to improve clarity.
#
# ### Mass conservation equation
# We multiply the equation by a test function ``\phi_{h}`` and integrate on a control volume ``\Omega``. After an integration by parts,
# we obtain (for an explicit integration time scheme):
# ```math
# \int_{\Omega} h^{n+1} \phi_{h} \mathrm{\,d}\Omega = \int_{\Omega} h^n \phi_{h} \mathrm{\,d}\Omega
# + \Delta t \left[ \int_{\Omega} h u^n \cdot \nabla \phi_{h} \mathrm{\,d}\Omega
# - \oint_{\Gamma} F_{h}^*(h^n, h u^n) \phi_{h} \mathrm{\,d} \Gamma \right]
# ```
# where ``F^*_{h}`` is the numerical flux corresponding to ``hu``.
#
# ### Momentum conservation equation
# We first consider the case without contact line forces nor curvature. Multiplying by a test function ``\phi_{h u}`` and integrating
# by parts leads to:
# ```math
#   \int_{\Omega} h u^{n+1} \phi_{h u} \mathrm{\,d}\Omega = \int_{\Omega} h u^n \phi_{h u} \mathrm{\,d}\Omega
#   + \Delta t \left[
#     \int_{\Omega} \left[
#         \mathcal{F}^n \cdot \nabla \phi_{h u}
#         + \left( h^n(g_t - \nabla P_g) - \tau_w + \tau_a \right) \phi_{h u}
#     \right] \mathrm{\,d}\Omega
#     - \oint_{\Gamma} F_{h u}^*(h^n, h u^n) \phi_{h u} \mathrm{\,d} \Gamma
#   \right]
# ```
# where ``F^*_{h u}`` is the numerical flux corresponding to ``h u \otimes u + g_n h^2 /2 \mathcal{I}``.
#

const dir = string(@__DIR__, "/")
using Bcube
using BcubeVTK
using LinearAlgebra
using StaticArrays
using BenchmarkTools
using Roots
using SparseArrays
using SparseDiffTools
using Profile
#using Symbolics
using InteractiveUtils

const eps_h = 1.0e-10

function compute_residual(_q, V, params, cache)
    # alias on measures
    dΓ = params.dΓ
    dΩ = params.dΩ
    dΓ_perio_x = params.dΓ_perio_x
    dΓ_perio_y = params.dΓ_perio_y

    q = get_fe_functions(_q)

    # face normals for each face domain (lazy, no computation at this step)
    n_Γ = get_face_normals(dΓ)
    nΓ_perio_x = get_face_normals(dΓ_perio_x)
    nΓ_perio_y = get_face_normals(dΓ_perio_y)

    function l(v)
        ∫(flux_Ω(q, v))dΩ +
        -∫(flux_Γ(q, v, n_Γ))dΓ +
        -∫(flux_Γ(q, v, nΓ_perio_x))dΓ_perio_x +
        -∫(flux_Γ(q, v, nΓ_perio_y))dΓ_perio_y
    end

    rhs = assemble_linear(l, V)
    return rhs
end

#velocity(h,hu) = (hu/max(h,eps_h))*(h>eps_h)
velocity(h, hu) = (hu * 2 * h) / (h * h + max(h * h, 1.0e-6))  #desingularization

"""
    flux_Ω(q, v)

Compute volume residual using the lazy-operators approach
"""
flux_Ω(q, v) = _flux_Ω ∘ (q, map(∇, v))

function _flux_Ω(q, ∇v)
    ∇λ_h, ∇λ_hu = ∇v
    f_h, f_hu = flux_sw(q)
    return ∇λ_h ⋅ f_h + ∇λ_hu ⊡ f_hu
end

function flux_sw(q)
    h, hu = q
    u = velocity(h, hu)
    huu = hu * transpose(u)
    g = stateInit.gravity
    p_grav = 0.5 * g * h * h
    return h .* u, huu + p_grav * I
end

"""
    flux_Γ(q, v, n)

Flux at the interface is defined by a composition of two functions:
* the input states which are needed for the flux using operator notations
* flux_rusanov(q, v, n) defines the face flux for values returned by facevar (as usual)
"""
flux_Γ(q, v, n) = flux_HLL ∘ (side⁻(q), side⁺(q), jump(v), side⁻(n))

function flux_HLL(q1, q2, δv, n12)
    g = stateInit.gravity
    δv_h, δv_hu = δv

    f_λ = x -> shallow_water_eigval(x, n12, g)
    flux = _flux_HLL(q1, q2, n12, flux_sw, f_λ)

    flux_h, flux_hu = flux
    return flux_h ⋅ δv_h + flux_hu ⋅ δv_hu
end

function _flux_rusanov(a, b, n, flux, f_λ)
    λ = max(f_λ(a), f_λ(b))
    f_rusanov(a, b, fa, fb) = 0.5 * (dotn(fa + fb, n) - λ * (b - a))
    map(f_rusanov, a, b, flux(a), flux(b))
end

function _flux_HLL(qL, qR, n, flux, f_λ)
    λL, λR = f_λ(qL), f_λ(qR)
    λ⁻ = min(minimum(λL), minimum(λR), zero(λL[1]))
    λ⁺ = max(maximum(λL), maximum(λR), zero(λL[1]))
    function f_HLL(qL, qR, fL, fR)
        if abs(λ⁺ - λ⁻) > 1.0e-12
            fLn, fRn = dotn(fL, n), dotn(fR, n)
            f = (λ⁺ * fLn - λ⁻ * fRn + λ⁻ * λ⁺ * (qR - qL)) / (λ⁺ - λ⁻)
        else
            f = 0.5 * (fL(qL) + fR(qR))
        end
        return f
    end
    map(f_HLL, qL, qR, flux(qL), flux(qR))
end

dotn(f::AbstractVector, n::AbstractVector) = f ⋅ n
dotn(f::AbstractMatrix, n::AbstractVector) = f * n

"""
    flux_Γ_wall(q, v, n)
"""
flux_Γ_wall(q, v, n) = flux_HLL ∘ (side⁻(q), side⁻(q), side⁻(v), side⁻(n))

# function _flux_Γ_wall(q1, v1, n12)
#     g = stateInit.gravity
#     h1, hu1 = q1
#     λ_h1, λ_hu1 = v1

#     flux_h  = zero(h1)
#     flux_hu = 0.5 * g * h1^2 * n12

#     return λ_h1 * flux_h + λ_hu1 * flux_hu
# end

function rhs(u, U, V, params, cache)
    rhs = compute_residual(u, V, params, cache)
    return cache.mass \ rhs
end

""" Inversion of mass matrix (expensive version!!) """
function compute_mass_matrix(U, V, dΩ)
    m(u, v) = ∫(u ⋅ v)dΩ
    M = assemble_bilinear(m, U, V)
    return factorize(M)
end

"""
    rk3_ssp(q, f::Function, t, Δt)

`f(q, t)` is the function to integrate.
"""
function rk3_ssp(q, f::Function, t, Δt)
    stepper(q, t) = forward_euler(q, f, t, Δt)
    _q0 = get_dof_values(q)

    _q1 = stepper(q, Δt)

    set_dof_values!(q, _q1)
    _q2 = (3 / 4) .* _q0 .+ (1 / 4) .* stepper(q, t + Δt)

    set_dof_values!(q, _q2)
    _q1 .= (1 / 3) * _q0 .+ (2 / 3) .* stepper(q, t + Δt / 2)

    return _q1
end

"""
Time integration of `f(q, t)` over a timestep `Δt`.
"""
forward_euler(q, f::Function, t, Δt) = get_dof_values(q) .+ Δt .* f(q, t)

mutable struct VtkHandler
    basename::String
    basename_residual::String
    ite::Int
    VtkHandler(basename) = new(basename, basename * "_residual", 0)
end

"""
    Write solution (at cell centers) to vtk
    Wrapper for `write_vtk`
"""
function append_vtk(vtk, mesh, vars, t, params)
    h, hu = vars

    vtk_degree = maximum(x -> get_degree(Bcube.get_function_space(get_fespace(x))), vars)
    vtk_degree = max(1, mesh_degree, vtk_degree)

    u = velocity ∘ (h, hu)

    # Write
    dict_vars_dg = Dict(
        "h" => h,
        "hu" => hu,
        "u" => u,
        "h_mean" => cell_mean(h, params.dΩ),
        "hu_mean" => cell_mean(hu, params.dΩ),
        "lim_h" => params.limh,
        "centers" => MeshCellData(params.xc),
    )

    write_file(
        vtk.basename * "_DG.pvd",
        mesh,
        dict_vars_dg,
        vtk.ite,
        t;
        mesh_degree = vtk_degree,
        collection_append = vtk.ite > 0,
    )

    # Update counter
    vtk.ite += 1
    return nothing
end

function init!(q, dΩ, initstate)
    x0 = SA[initstate.x0, initstate.y0]
    Lstep = initstate.Lstep
    hstep = initstate.hstep

    f_h(x) = norm(x - x0) < Lstep / 2 ? hstep : 0.0
    # f_h(x)  = hstep * exp.(-(abs(x[1] - x0) ./ Lstep)^2 ./ 2)
    f_hu(x) = SA[0.0, 0.0]
    f = map(PhysicalFunction, (f_h, f_hu))
    projection_l2!(q, f, dΩ)
    return nothing
end

function run_simulation(stateInit)
    # Then generate the mesh of a rectangle using Gmsh and read it
    tmp_path = "tmp.msh"
    nx, ny, lx, ly = stateInit.nx, stateInit.ny, stateInit.lx, stateInit.ly
    gen_rectangle_mesh(
        tmp_path,
        :quad;
        nx = nx,
        ny = ny,
        lx = lx,
        ly = ly,
        xc = 0.0,
        yc = 0.0,
    )
    mesh = read_msh(tmp_path, 2) # '2' indicates the space dimension (3 by default)
    rm(tmp_path)

    dimcar = compute_dimcar(mesh)

    # Create a `CellVariable`
    fs = FunctionSpace(fspace, degree)
    Q_sca = TrialFESpace(fs, mesh, :discontinuous; size = 1) # DG, scalar
    Q_vec = TrialFESpace(fs, mesh, :discontinuous; size = 2) # DG, vectoriel
    V_sca = TestFESpace(Q_sca)
    V_vec = TestFESpace(Q_vec)

    Q = MultiFESpace(Q_sca, Q_vec)
    V = MultiFESpace(V_sca, V_vec)
    q = FEFunction(Q)

    # select an initial configurations:
    init!(q, mesh, stateInit)

    DMPrelax = DMPcurv₀ .* dimcar .^ 2

    # Then we create a `NamedTuple` to hold the simulation parameters.
    params = (degree = degree, stateInit = stateInit, DMPrelax = DMPrelax)

    # Define measures for cell and interior face integrations
    dΩ = Measure(CellDomain(mesh), QUAD)
    dΓ = Measure(InteriorFaceDomain(mesh), QUAD)

    # Declare boundary conditions and
    # create associated domains and measures
    periodicBCType_x = PeriodicBCType(Translation(SA[-lx, 0.0]), ("East",), ("West",))
    periodicBCType_y = PeriodicBCType(Translation(SA[0.0, ly]), ("South",), ("North",))
    Γ_perio_x = BoundaryFaceDomain(mesh, periodicBCType_x)
    Γ_perio_y = BoundaryFaceDomain(mesh, periodicBCType_y)
    dΓ_perio_x = Measure(Γ_perio_x, QUAD)
    dΓ_perio_y = Measure(Γ_perio_y, QUAD)

    xc = get_cell_centers(mesh) # used for VTK outputs

    params = (
        params...,
        dΓ = dΓ,
        dΩ = dΩ,
        dΓ_perio_x = dΓ_perio_x,
        dΓ_perio_y = dΓ_perio_y,
        xc = xc,
    )

    # create CellData to store limiter values
    limh = MeshCellData(ones(ncells(mesh)))
    limAll = MeshCellData(ones(ncells(mesh)))
    params = (params..., limh = limh, limAll = limAll)

    # Init vtk handler
    mkpath(outputpath)
    vtk = VtkHandler(output)

    # Init time
    time = 0.0

    # cache mass matrices
    cache = (
        mass = compute_mass_matrix(Q, V, dΩ),
        mass_sca = compute_mass_matrix(Q_sca, V_sca, dΩ),
        mass_vec = compute_mass_matrix(Q_vec, V_vec, dΩ),
        cacheCellMean = Bcube.build_cell_mean_cache(q, dΩ),
    )
    Δt = Δt₀

    limiter_projection && apply_limitation!(q, params, cache)

    # Save initial solution
    append_vtk(vtk, mesh, q, time, params)

    # Let's loop to solve the equation.
    for i in 1:nite
        Δt = compute_timestep!(q, dΩ, dimcar, CFL)

        ## Infos
        if (i % Int(max(floor(nite / (nout * 10)), 1)) == 0)
            println("---")
            println("Iteration ", i)
            @show Δt, CFL
        end

        ## Step forward in time
        _rhs(q, t) = rhs(q, Q, V, params, cache)
        if timeScheme == :ForwardEuler
            qnew = forward_euler(q, _rhs, time, Δt)
        elseif timeScheme == :RK3_SPP
            qnew = rk3_ssp(q, _rhs, time, Δt)
        else
            error("Unknown time scheme")
        end
        set_dof_values!(q, qnew)

        limiter_projection && apply_limitation!(q, params, cache)

        time += Δt

        ## Write solution to file (respecting max. number of output)
        if (i % Int(max(floor(nite / nout), 1)) == 0)
            append_vtk(vtk, mesh, q, time, params)
        end
    end

    append_vtk(vtk, mesh, q, time, params)

    println("Benchmarking 'forward_euler':")
    _rhs1(q, t) = rhs(q, Q, V, params, cache)
    @btime forward_euler($q, $_rhs1, $time, $Δt)
    @btime apply_limitation!($q, $params, $cache)
    println("ndofs total = ", Bcube.get_ndofs(Q))

    Profile.init(; n = 10^7) # returns the current settings
    Profile.clear()
    Profile.clear_malloc_data()
    @profile begin
        for i in 1:100
            forward_euler(q, _rhs1, time, Δt)
            limiter_projection && apply_limitation!(q, params, cache)
        end
    end
end

function compute_timestep!(q, dΩ, dimcar, CFL)
    q_mean = map(Base.Fix2(cell_mean, dΩ), get_fe_functions(q))
    λ(x...) = shallow_water_maxeigval(x, stateInit.gravity)
    λmax = map(λ, get_values.(q_mean)...)
    Δt = CFL * minimum(dimcar ./ λmax)
    return Δt
end

function shallow_water_maxeigval(q, gravity)
    h, hu = q
    u = velocity(h, hu)
    return norm(u) + √(abs(gravity) * max(h, eps_h))
end

function shallow_water_maxeigval(q, n, gravity)
    h, hu = q
    un = velocity(h, hu) ⋅ n
    return abs(un) + √(abs(gravity) * max(h, eps_h))
end

function shallow_water_eigval(q, n, gravity)
    h, hu = q
    un = velocity(h, hu) ⋅ n
    c = √(abs(gravity) * max(h, eps_h))
    return un - c, un + c
end

function compute_dimcar(mesh)
    fs = FunctionSpace(:Lagrange, 0)
    V = TestFESpace(fs, mesh; size = 1, isContinuous = false)

    # Define measures for cell and interior face integrations
    dΩ = Measure(CellDomain(mesh), QUAD)
    dΓ = Measure(InteriorFaceDomain(mesh), QUAD)
    dΓ_bc = Measure(BoundaryFaceDomain(mesh), QUAD)

    f1 = PhysicalFunction(x -> 1.0)
    l(v) = ∫(f1 ⋅ v)dΩ
    l_face(v, dω) = ∫(side⁻(f1) ⋅ side⁻(v) + side⁺(f1) ⋅ side⁺(v))dω

    vol = assemble_linear(l, V)
    surf = assemble_linear(Base.Fix2(l_face, dΓ), V)
    assemble_linear!(surf, Base.Fix2(l_face, dΓ_bc), V)
    return vol ./ surf
end

function apply_limitation!(q, params, cache)
    h, hu = q
    dΩ = params.dΩ

    q_mean = cell_mean(q, cache.cacheCellMean)

    _limh, _h_proj = linear_scaling_limiter(
        h,
        params.dΩ;
        bounds = (hmin₀, hmax₀),
        DMPrelax = params.DMPrelax,
        mass = cache.mass_sca,
    )
    set_dof_values!(h, get_dof_values(_h_proj))

    h_mean, hu_mean, = q_mean
    limited_var(a, a̅, lim_a) = a̅ + lim_a * (a - a̅)
    projection_l2!(hu, limited_var(hu, hu_mean, _limh), params.dΩ; mass = cache.mass_vec)

    if eltype(_limh) == eltype(params.limh) # skip Dual number case
        set_values!(params.limh, get_values(_limh))
    end
    return nothing
end

# Function-space get_degree (Taylor(0) = first order Finite Volume)
const degree = 1
const mesh_degree = 1
const fspace = Bcube.Lagrange(:Lobatto) #:Lagrange
# The degree of quadrature is chosen such as mass matrix are integrated exactly.
const QUAD = Quadrature(QuadratureLobatto(), max(1, 2 * degree))
const limiter_projection = true
const hmin₀ = 1.e-8
const hmax₀ = 1.0e10
const DMPcurv₀ = 10.0

const stateInit = (
    gravity = 9.81,
    nx = 65,
    ny = 65,
    lx = 2.0,
    ly = 2.0,
    x0 = 0.0,
    y0 = 0.0,
    Lstep = 1.0,
    hstep = 2.5,
)
const nite = 5000 #300000 # Number of time iteration(s)
const timeScheme = :ForwardEuler # :ForwardEuler, :RK3_SPP
const CFL = 0.4 / (2 * degree + 1)
const nout = 100 # Number of time steps to save
const outputpath = joinpath(@__DIR__, "../../../myout/shallow_water/")
const output = outputpath * "sw_deg$degree"
const Δt₀ = 1.e-7

@show get_degree(QUAD)

mkpath(outputpath)
run_simulation(stateInit)

end #hide
