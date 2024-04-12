"""
Custom way to "include" a file to print infos.
"""
function custom_include(path)
    filename = split(path, "/")[end]
    basename = first(splitext(filename))
    @testset "$basename" begin
        include(path)
    end
    println("done.")
end

#------------------------------------------------#
#  Functions for tests based on checksum         #
#------------------------------------------------#

"""
    scientific_format(x, digits = 10)

Return the value of `x` as a `string` in a scientific format
(decimal exponent notation). Optional argument `digits` specifies
the exact number of digits to appear after the decimal point character.
The default precision is `10`.

# Examples
```julia-repl
julia> scientific_format(-0.00001)
"-1.0000000000e-05"
julia> scientific_format(pi)
"3.1415926536e+00"
```
"""
scientific_format(x, digits = 10) = format(Format("%10.$(digits)e"), x)

"""
    get_ref_checksum(key)

Read the checksum value corresponding to the `key` in
file "checksums.sha256" in the current direction.
If `key` is not found, it returns `nothing`.
"""
function get_ref_checksum(key)
    f = readdlm(joinpath(@__DIR__, "checksums.sha256"), String)
    fname2sum = Dict(r[2] => r[1] for r in eachrow(f))
    return get(fname2sum, key, nothing)
end

"""
    compute_checksum(array::AbstractArray{<:Number}; digits = 10)

Compute the sha256 checksum of the given `array`, which be must be
a subtype of `AbstractArray{<:Number}`. The array values are rounded
and converted to string in decimal exponent notation before
computing the checksum. Optional argument `digits` specifies the exact
number of digits that are kept after the decimal point character.
The default precision is `10`.
"""
function compute_checksum(array::AbstractArray{<:Number}; digits = 10)
    path, _ = mktemp()
    str_array = map(Base.Fix2(scientific_format, digits), round.(array, digits = digits))
    writedlm(path, str_array)
    checksum = bytes2hex(open(sha256, path))
    return checksum
end

"""
    compare_checksum(key, value; digits = 10, verbose = true)

Compare the checksum of `value` with the reference checksum
corresponding to the `key` in file "checksums.sha256".
Keyword argument `digits` can be use to control rounding of `value`
before computing the checksum. This option allows to control test
failure due to floating-point precision.
Note that the reference checksum must be generated with the same number of `digits`.
"""
function compare_checksum(key, value; digits = 10, verbose = true)
    _ref = get_ref_checksum(key)
    _cur = compute_checksum(value; digits = digits)
    result = _cur == _ref
    if (!result) && verbose
        printstyled("   ➡ TEST FAIL : key = ", key, " \n"; color = :light_black)
        printstyled("      ↪ checksum ref = ", _ref, " \n"; color = :light_black)
        printstyled("      ↪ checksum     = ", _cur, " \n"; color = :light_black)
    end
    return result
end

#------------------------------------------------#
#  Functions for tests based on ReferenceTest.jl #
#------------------------------------------------#

"""
    refpath(filename::String)

Return the full path of the reference file named `filename`.
It assumes that the file is located in the `test/references/`
directory of the current project.
"""
refpath(filename::String) = joinpath(@__DIR__, "references/", filename)

"""
    compare(
        d1::Dict{String, <:AbstractArray},
        d2::Dict{String, <:AbstractArray},
        rtol,
        atol,
    )

Return `true` if dictionaries `d1` and `d2` have the
same keys and if these keys are associated with arrays
of equal values, with a relative tolerance `rtol` and
an absolute tolerance `atol`.
"""
function compare(
    d1::Dict{String, <:AbstractArray},
    d2::Dict{String, <:AbstractArray},
    rtol,
    atol,
)
    # From the official doc for `Dict`:
    # -
    # "When the values are stored internally in a hash table,
    # as is the case for Dict, the order in which they are returned may vary".
    # -
    # Then `keys(d1) == keys(d2)` is not recommended
    # and we do the following two tests instead :
    length(keys(d1)) != length(keys(d2)) && return false
    any(x -> !haskey(d1, x), keys(d2)) && return false

    for key in keys(d1)
        any((!isapprox).(d1[key], d2[key]; atol = atol, rtol = rtol)) && return false
    end
    return true
end
compare(; atol::Real = 1.0e-12, rtol::Real = 1.0e-12) = (a, b) -> compare(a, b, atol, rtol)

"""
    test_ref(filename_ref::String, data, comp::Function = compare())

Test that the values `data` with reference `filename_ref` (stored in `./test/references`)
using equality test strategy given by `comp`.

By default, `comp=compare()` assumes `data` is a subtype of
`Dict{String, <:AbstractArray}`.

If `data` is not a valid subtype for the defaut test strategy, one could :
- defined the method `_as_dict(data)` to convert `data` to a valid subtype.
- provide another test function with the keyword argument `comp`.
The function must have the same signature as `Base.isequal` function.
"""
function test_ref(filename_ref::String, data, comp::Function = compare())
    @test_reference refpath(filename_ref) _as_dict(data) by = comp
end
_as_dict(a::Dict) = a
_as_dict(a::AbstractArray) = Dict("array" => a)
_as_dict(a::AbstractSparseMatrix) = Dict(zip("sparse_" .* ("I", "J", "V"), findnz(a)))