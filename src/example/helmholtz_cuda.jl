# Solving Helmholtz eigenvalue problem on GPU
module Helmholtz

using Bcube
using LinearAlgebra
using LinearMaps
using Krylov
using KrylovKit
using CUDA
using CUDA.CUSPARSE
CUDA.allowscalar(false)

struct ILU0gpu{T}
    ILU0::T
    ILU0gpu(A) = new{typeof(A)}(A)
end

# import LinearAlgebra
function LinearAlgebra.ldiv!(
    y::CuArray{T},
    P::ILU0gpu{CuSparseMatrixCSR{T, M}},
    x::CuArray{T},
) where {T, M}
    A = P.ILU0
    copyto!(y, x)
    sv2!('N', 'L', 'U', 1.0, A, y, 'O')
    sv2!('N', 'U', 'N', 1.0, A, y, 'O')
    return y
end

const σ = π
const nev = 10
# mesh = rectangle_mesh(21, 21)
mesh = line_mesh(1000)

degree = 1

U = TrialFESpace(FunctionSpace(:Lagrange, degree), mesh)
V = TestFESpace(U)

dΩ = Measure(CellDomain(mesh), 2 * degree + 1)

a(u, v) = ∫(∇(u) ⋅ ∇(v))dΩ
b(u, v) = ∫(u ⋅ v)dΩ

A = assemble_bilinear(a, U, V)
B = assemble_bilinear(b, U, V)

A = CuSparseMatrixCSR(A)
B = CuSparseMatrixCSR(B)
C = A - σ * B

x0 = CUDA.rand(Float64, get_ndofs(U))
P = ILU0gpu(ilu02(C, 'O'))
solver = GmresSolver(size(C)..., 20, typeof(x0))

linsolve =
    (_x, _A, _b) -> begin
        copyto!(_x, gmres!(solver, _A, _b; M = P, ldiv = true, restart = true).x)
    end

op = InverseMap(C; solver = linsolve) * B
vp, vecs, info = eigsolve(op, x0, nev, :LM; krylovdim = 20)
show(info)
vp = σ .+ 1 ./ vp
@show sqrt.(sort(abs.(vp[1:10])))

end