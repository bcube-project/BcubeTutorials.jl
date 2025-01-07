module HeatEquationTwoLayers #hide
println("Running heat equation two layers example...") #hide
# # Heat equation
# # Theory
# This example shows how to solve the heat equation with eventually variable physical properties in steady and unsteady formulations:
# ```math
#   \rho C_p \partial_t u - \nabla . ( \lambda \nabla u) = f
# ```
# We shall assume that $$f, \, \rho, \, C_p, \, \lambda \, \in L^2(\Omega)$$. The weak form of the problem is given by: find $$ u \in \tilde{H}^1_0(\Omega)$$
# (there will be at least one Dirichlet boundary condition) such that:
# ```math
#   \forall v \in  \tilde{H}^1_0(\Omega), \, \, \, \underbrace{\int_\Omega \partial_t u . v dx}_{m(\partial_t u,v)} + \underbrace{\int_\Omega \nabla u . \nabla v dx}_{a(u,v)} = \underbrace{\int_\Omega f v dx}_{l(v)}
# ```
# To numerically solve this problem we seek an approximate solution using Lagrange $$P^1$$ or $$P^2$$ elements.
# Here we assume that the domain can be split into two domains having different material properties.

const dir = joinpath(@__DIR__, "..", "..", "..") # BcubeTutorials dir
using Bcube
using BcubeVTK
using LinearAlgebra
using Test #src

function T_analytical_two_layers(x, λ1, λ2, T0, T1, L)
    if x < 0.5 * L
        T = 2.0 * (λ2 / (λ1 + λ2)) * ((T1 - T0) / L) * x + T0
    else
        T = 2.0 * (λ1 / (λ1 + λ2)) * ((T1 - T0) / L) * (x - L) + T1
    end
    return T
end

const outputpath = joinpath(dir, "myout", "heat_equation_two_layers")
mkpath(outputpath)

"""
Two layers case using MeshCellData that depends on the subdomain. A unique assembly is then performed over the whole mesh
(by opposition to method2)
"""
function run_steady_two_layers_method1(; degree)
    println("Running steady two layers case - method 1")
    T0 = 260.0
    T1 = 300.0
    λ1 = 150.0
    λ2 = 10.0

    # Read mesh
    mesh = read_mesh(joinpath(dir, "input", "mesh", "domainTwoLayer_tri.msh"))

    fs = FunctionSpace(:Lagrange, degree)
    U = TrialFESpace(fs, mesh, Dict("West" => T0, "East" => T1))
    V = TestFESpace(U)

    # Define measures for cell integration
    dΩ = Measure(CellDomain(mesh), 2 * degree + 1)

    material_1 = Bcube.get_zone_element_indices(mesh, "Domain_1")
    material_2 = Bcube.get_zone_element_indices(mesh, "Domain_2")

    λ = zeros(ncells(mesh))
    λ[material_1] .= λ1
    λ[material_2] .= λ2

    qtmp = zeros(ncells(mesh))

    q = MeshCellData(qtmp)
    η = MeshCellData(λ)

    # Compute matrices associated to bilinear and linear forms
    a(u, v) = ∫(η * ∇(u) ⋅ ∇(v))dΩ
    l(v) = ∫(q * v)dΩ

    system = AffineFESystem(a, l, U, V)
    ϕ = Bcube.solve(system)

    T_analytical = PhysicalFunction(x -> T_analytical_two_layers(x[1], λ1, λ2, T0, T1, 0.2))

    Tcn = var_on_centers(ϕ, mesh)
    Tca = var_on_centers(T_analytical, mesh)
    dict_vars =
        Dict("Temperature" => MeshCellData(Tcn), "Temperature_a" => MeshCellData(Tca))
    write_file(
        joinpath(outputpath, "result_steady_heat_equation_two_layers_method1.pvd"),
        mesh,
        dict_vars,
    )

    err = norm(Tca .- Tcn, Inf) / norm(Tca, Inf)
    @show err

    if get(ENV, "TestMode", "false") == "true" #src
        @test err < 2.1e-5                     #src
    end                                        #src
end

"""
Two layers case using subdomain integration.
"""
function run_steady_two_layers_method2(; degree)
    println("Running steady two layers case - method 2")
    T0 = 260.0
    T1 = 300.0
    λ1 = 150.0
    λ2 = 10.0

    # Read mesh
    mesh = read_mesh(joinpath(dir, "input", "mesh", "domainTwoLayer_tri.msh"))

    fs = FunctionSpace(:Lagrange, degree)
    U = TrialFESpace(fs, mesh, Dict("West" => T0, "East" => T1))
    V = TestFESpace(U)

    material_1 = Bcube.get_zone_element_indices(mesh, "Domain_1")
    material_2 = Bcube.get_zone_element_indices(mesh, "Domain_2")

    # Define measures for cell and interior face integrations

    dΩ = Measure(CellDomain(mesh), 2 * degree + 1)
    dΩ1 = Measure(CellDomain(mesh, material_1), 2 * degree + 1)
    dΩ2 = Measure(CellDomain(mesh, material_2), 2 * degree + 1)

    qtmp = zeros(ncells(mesh))

    q = MeshCellData(qtmp)

    # compute matrices associated to bilinear and linear forms
    a(u, v) = ∫(λ1 * ∇(u) ⋅ ∇(v))dΩ1 + ∫(λ2 * ∇(u) ⋅ ∇(v))dΩ2
    l(v) = ∫(q * v)dΩ

    system = AffineFESystem(a, l, U, V)
    ϕ = Bcube.solve(system)

    T_analytical = PhysicalFunction(x -> T_analytical_two_layers(x[1], λ1, λ2, T0, T1, 0.2))

    Tcn = var_on_centers(ϕ, mesh)
    Tca = var_on_centers(T_analytical, mesh)
    dict_vars =
        Dict("Temperature" => MeshCellData(Tcn), "Temperature_a" => MeshCellData(Tca))
    write_file(
        joinpath(outputpath, "result_steady_heat_equation_two_layers_method2.pvd"),
        mesh,
        dict_vars,
    )

    err = norm(Tca .- Tcn, Inf) / norm(Tca, Inf)
    @show err

    if get(ENV, "TestMode", "false") == "true" #src
        @test err < 2.1e-5                     #src
    end                                        #src
end

run_steady_two_layers_method1(; degree = 2)
run_steady_two_layers_method2(; degree = 2)

end #hide
