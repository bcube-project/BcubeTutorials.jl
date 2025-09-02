module FlowAroundCylinderCBS
using Bcube, BcubeGmsh, BcubeVTK
using LinearAlgebra
using StaticArrays

function run_inviscid()
    degree_u = 1
    degree_p = 1
    θ1 = 1.0
    θ2 = 1.0
    f_inlet(x) = (0.05 ≤ x[2] ≤ 0.36) ? -1.0 : 0.0
    un_inlet(t) = PhysicalFunction(x -> f_inlet(x)) # u ⋅ n on inlet
    CFL = 0.1

    degmax = max(degree_u, degree_p)
    outdir = joinpath(@__DIR__, "..", "..", "..", "myout", "flow_around_cyl_cbs")
    outpath = joinpath(outdir, "flow_around_cyl_cbs.pvd")
    rm(outdir; force = true, recursive = true)
    mkpath(outdir)

    meshpath = joinpath(
        @__DIR__,
        "..",
        "..",
        "..",
        "input",
        "mesh",
        "cylinder_navier_stokes_tri.msh",
    )
    mesh = read_mesh(meshpath)
    dΩ = Measure(CellDomain(mesh), degmax)
    dΓ_inlet = Measure(BoundaryFaceDomain(mesh, "left"), degmax)
    dΓ_outlet = Measure(BoundaryFaceDomain(mesh, "right"), degmax)
    # dΓ_farfield = Measure(BoundaryFaceDomain(mesh, ("top", "bottom")), degmax)
    # dΓ_wall = Measure(BoundaryFaceDomain(mesh, "wall"), degmax)
    dΓ_wall = Measure(BoundaryFaceDomain(mesh, ("top", "bottom", "cylinder")), degmax)
    nΓ_inlet = get_face_normals(dΓ_inlet)
    nΓ_wall = get_face_normals(dΓ_wall)
    h = minimum(compute_dimcar(mesh))
    @show h

    Δt = CFL * h / 2
    @show Δt

    diri_u = Dict{String, Any}(
        "left" => (x, t) -> f_inlet(x),
        "top" => 0.0,
        "bottom" => 0.0,
        "cylinder" => 0.0,
    )
    diri_p = Dict{String, Any}("right" => 0.0)

    U_u = TrialFESpace(FunctionSpace(:Lagrange, degree_u), mesh, diri_u; size = 2)
    U_p = TrialFESpace(FunctionSpace(:Lagrange, degree_p), mesh, diri_p)
    V_u, V_p = TestFESpace.((U_u, U_p))

    # U = MultiFESpace(U_u, U_p)
    # V = MultiFESpace(TestFESpace.(U_u, U_p)...)
    # q = FEFunction(U)

    a_Mu(u, φ) = ∫(u ⋅ φ)dΩ
    Mu = assemble_bilinear(a_Mu, U_u, V_u)
    Mu_fac = factorize(Mu)

    a_H(p, φ) = ∫(∇(p) ⋅ ∇(φ))dΩ
    H = assemble_bilinear(a_H, U_p, V_p)
    H_fac = copy(H)
    apply_dirichlet_to_matrix!(H_fac, U_p, V_p, mesh)
    H_fac = factorize(H_fac)

    a_G(φ_u, φ_p) = ∫(∇(φ_p) ⋅ φ_u)dΩ
    G = assemble_bilinear(a_G, U_u, V_p) # Note that TrialFESpace=U_u and TestFESpace=V_p

    sn = side_n
    α = 0.0 / h
    function a_D(u, φ)
        ∫(α * (sn(u) ⋅ sn(nΓ_inlet)) * (sn(φ) ⋅ sn(nΓ_inlet)))dΓ_inlet +
        ∫(α * (sn(u) ⋅ sn(nΓ_wall)) * (sn(φ) ⋅ sn(nΓ_wall)))dΓ_wall
    end
    Du = assemble_bilinear(a_D, U_u, V_u)
    l_d(φ) = ∫(-α * sn(un_inlet(0.0)) * (sn(φ) ⋅ sn(nΓ_inlet)))dΓ_inlet
    du = assemble_linear(l_d, V_u)

    A_step1 = Mu + Δt * Du
    apply_un!(A_step1, U_u, V_u, mesh)
    A_step1 = factorize(A_step1)

    u = FEFunction(U_u)
    p = FEFunction(U_p)

    # projection_l2!(u, PhysicalFunction(x -> SA[1.0, 0.0]), mesh)

    t = 0.0
    write_file(outpath, mesh, Dict("u" => u, "p" => p))

    nitemax = 2
    nout = max(round(Int, nitemax / 500), 1)
    for i in 1:nitemax
        udofs = get_dof_values(u)
        pdofs = get_dof_values(p)

        # Step 1
        Cu = assemble_linear(φ -> ∫(f_Cu ∘ (u, ∇(u), φ))dΩ, V_u)
        Ku = assemble_linear(φ -> ∫(f_Ku ∘ (u, ∇(u), φ, ∇(φ)))dΩ, V_u)
        b = -Δt * (Cu + du - Δt * Ku)
        apply_un!(b, U_u, V_u, mesh)
        Δu = A_step1 \ b
        @show maximum(abs.(Δu))

        # Step 2
        # ustar = FEFunction(U_u, get_dof_values(u) + θ1 * Δu_star)
        # G = assemble_linear(φ -> ∫(f_G ∘ (ustar, ∇(φ))))
        fp =
            Δt *
            assemble_linear(φ -> ∫(f_fp ∘ (side_n(φ), side_n(un_inlet(t))))dΓ_inlet, V_p) # Note the multiplying Δt
        b = 1 / (Δt * θ1 * θ2) * (G * (udofs + θ1 * Δu) - Δt * θ1 * H * pdofs - fp)
        apply_dirichlet_to_vector!(b, U_p, V_p, mesh)
        Δp = H_fac \ b
        @show maximum(abs.(Δp))

        # Step 3
        P = assemble_linear(φ -> ∫(f_P ∘ (u, ∇(u), φ, ∇(φ), ∇(p)))dΩ, V_u)
        Δu_tilde = Δu - Mu_fac \ (Δt * (transpose(G) * (pdofs + θ2 * Δp) + Δt / 2 * P))
        @show maximum(abs.(Δu_tilde))

        # Update solution
        set_dof_values!(u, udofs + Δu_tilde)
        set_dof_values!(p, pdofs + Δp)
        t += Δt

        # Write
        if i % nout == 0
            write_file(
                outpath,
                mesh,
                Dict("u" => u, "p" => p),
                i,
                t;
                collection_append = true,
            )
        end
    end
    write_file(
        outpath,
        mesh,
        Dict("u" => u, "p" => p),
        nitemax,
        t;
        collection_append = true,
    )
end

"""
Corresponds to φ ⋅ [∇⋅(u⊗u)]

Rq : div(a⊗b)= ∇a*b + a div(b)
"""
f_Cu(u, ∇u, φ) = φ ⋅ (∇u * u + u * tr(∇u))

"""
Corresponds to [∇⋅(φ⊗u)] ⋅ [∇⋅(u⊗u)]

Rq : div(a⊗b)= ∇a*b + a div(b)
"""
f_Ku(u, ∇u, φ, ∇φ) = -0.5 * (∇φ * u + φ * tr(∇u)) ⋅ (∇u * u + u * tr(∇u))

f_G(u, ∇φ_p) = ∇φ_p ⋅ u

f_fp(φ_p, un) = φ_p * un

"""
Corresponds to [∇⋅(φ⊗u)] ⋅ ∇p

Rq : div(a⊗b)= ∇a*b + a div(b)
"""
f_P(u, ∇u, φ, ∇φ, ∇p) = (∇φ * u + φ * tr(∇u)) ⋅ ∇p

div(∇u) = tr(∇u)

function compute_dimcar(mesh)
    fs = FunctionSpace(:Lagrange, 0)
    V = TestFESpace(fs, mesh; size = 1, isContinuous = false)

    # Define measures for cell and interior face integrations
    degquad = 1
    dΩ = Measure(CellDomain(mesh), degquad)
    dΓ = Measure(InteriorFaceDomain(mesh), degquad)
    dΓ_bc = Measure(BoundaryFaceDomain(mesh), degquad)

    f1 = PhysicalFunction(x -> 1.0)
    l(v) = ∫(f1 ⋅ v)dΩ
    l_face(v, dω) = ∫(side⁻(f1) ⋅ side⁻(v) + side⁺(f1) ⋅ side⁺(v))dω

    vol = assemble_linear(l, V)
    surf = assemble_linear(Base.Fix2(l_face, dΓ), V)
    surf += assemble_linear(Base.Fix2(l_face, dΓ_bc), V)
    return vol ./ surf
end

function apply_un!(
    array,
    U::TrialFESpace{S, FE},
    V::TestFESpace{S, FE},
    mesh,
    t::Number = 0.0,
) where {S, FE}
    # Alias
    _mesh = parent(mesh)
    fs_V = Bcube.get_function_space(V)
    dhl_V = Bcube._get_dhl(V)

    # Loop over the boundaries
    for bndTag in Bcube.get_dirichlet_boundary_tags(U)
        # Function to apply, giving the Dirichlet value(s) on each node
        f = Bcube.get_dirichlet_values(U, bndTag)
        f_t = Base.Fix2(f, t)

        # Loop over the face of the boundary
        for kface in Bcube.boundary_faces(_mesh, bndTag)
            apply_un_on_face!(array, kface, _mesh, fs_V, dhl_V, f_t)
        end
    end
end

function apply_un_on_face!(array, kface::Int, mesh, fs_V, dhl_V, f_t::Function)
    # Alias
    sizeV = Bcube.get_ncomponents(dhl_V)
    c2n = Bcube.connectivities_indices(mesh, :c2n)
    f2n = Bcube.connectivities_indices(mesh, :f2n)
    f2c = Bcube.connectivities_indices(mesh, :f2c)
    cellTypes = Bcube.cells(mesh)

    # Interior cell
    icell = f2c[kface][1]
    ctype = cellTypes[icell]
    _c2n = c2n[icell]
    cnodes = get_nodes(mesh, _c2n)
    side = Bcube.cell_side(ctype, c2n[icell], f2n[kface])
    cshape = Bcube.shape(ctype)
    ξcell = get_coords(fs_V, cshape) # ref coordinates of the FunctionSpace in the cell
    F = Bcube.mapping(ctype, cnodes)

    # CHEATING : FOR LINEAR 2D-ELEMENTS THE FACE-NORMALS ARE CONSTANT
    n = Bcube.normal(ctype, cnodes, side, SA[0.0])

    # local indices of dofs lying on the face (assuming scalar FE)
    idofs_loc = Bcube.idof_by_face_with_bounds(fs_V, Bcube.shape(ctype))[side]

    # Loop over the dofs concerned by the Dirichlet condition
    for idof_loc in idofs_loc
        ξ = ξcell[idof_loc]
        values = f_t(F(ξ)) # dirichlet value(s)

        # Loop over components
        for icomp in 1:sizeV
            # Absolute number of the dof for this component
            idof_glo = Bcube.get_dof(dhl_V, icell, icomp, idof_loc)

            if array isa AbstractMatrix
                if values isa Number
                    array[Bcube.get_dof(dhl_V, icell, 2, idof_loc), idof_glo] = n[icomp]
                else
                    array[:, idof_glo] .= 0.0
                    array[idof_glo, idof_glo] = 1.0
                end
            else
                array[Bcube.get_dof(dhl_V, icell, 2, idof_loc)] = values
            end
        end
    end
end

run_inviscid()
end