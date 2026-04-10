# BasicFunctions.jl
# This module contains general functions that are required by all other modules.
# v2.0

module BasicFunctions

#! format: off
export
	ft,
	ift,
	struct_fcn2_ft,
	meshgrid,
	zeropad,
	inner_circ,
	rect,
	cart2pol,
	pol2cart,
	print_paras,
	msum,
	norm_mode,
	wopt,
	rayleigh_waist,
	norm_coeff
#! format: on

using FFTW
using Printf
using Statistics


# I. Fourier transforms and some applications

"""
    ft(g, delta=1)

Multi dimensional Fourier transfrom with optional gridspacing `d = 1``.
"""
function ft(
    g::Array,
    delta::Number
    )
    return fftshift(fft(fftshift(g))) .* delta^ndims(g) 
end

"""
    ift(g, delta_f=1)

Inverse multi dimensional Fourier transfrom. See also [`ft`](@ref).
"""
function ift(
    g::Array,
    delta_f::Number
    )
    return ifftshift(ifft(ifftshift(g))) .* (size(g, 1) * delta_f)^ndims(g)
end

"""
    struct_fcn2_ft(ph, mask, delta)

Structure function based on FT from [Schmidt (2010)]. Very sensitive to
boundaries, needs zeropadding. Not tested. Recommended not to use.
"""
function struct_fcn2_ft(
    ph::Matrix,
    mask::Matrix,
    delta::Number
    )
    n = size(ph, 1)
    ph = ph .* mask
    p = ft(ph, delta)
    s = ft(ph.^2, delta)
    w = ft(mask, delta)
    delta_f = 1 / (n * delta)
    w2 = ift(w .* conj(w), delta_f)
    out = 2 * ift(real(s .* conj(w)) - abs.(p).^2, delta_f) ./ w2 .* mask
    return out
end


# II. Help functions

"""
    meshgrid(x, y)

Rectangular meshgrid array from the vectors `x` and `y`.
"""
function meshgrid(
    x::Vector,
    y::Vector
    )
    Nx = size(x, 1)
    Ny = size(y, 1)
    xx = x' .* ones(Ny)
    yy = ones(Nx)' .* y
    return xx, yy
end

"""
    zeropad(a, s)

Zeropad a matrix to the desired side length `s`.
"""
function zeropad(
    a::Matrix,
    s::Integer
    )
    n = size(a, 1)
    if s > n
        out = zeros(eltype(a), s, s)
        i1 = floor(Int, (s - n)/2) - 1
        i2 = i1 + n
        out[i1+1:i2, i1+1:i2] .= a
        return out
    else
        return a
    end
end

"""
    inner_circ(a; nans=false)

Set all elements outside of an inscribed circle to zero or NaN (`nans=true`).
"""
function inner_circ(
    a::Matrix;
    nans::Bool=false
    )
    out = zeros(eltype(a), size(a))
    if nans
        fill!(out, NaN)
    end
    n = size(a, 1)
    for i in 1:n
        for j in 1:n
            r = hypot(i - n/2 -.5, j - n/2- .5)
            if r <= n/2
                out[i, j] = a[i, j]
            end
        end
    end
    return out
end

"""
    rect(x; w=1, nans=false)

Rectangle function returns one if x is within
(-w/2, w/2) else zero or NaN (`nans=true`).
"""
function rect(
    x::Number;
    w::Number=1,
    nans::Bool=false
    )
    ax = abs(x)
    w2 = w/2
    if ax < w2
        return 1
    elseif ax == w2
        return 1/2
    elseif ax > w2
        if nans
            return NaN
        else
            return 0
        end
    end
end


"""
# III. Print and help functions
"""

"""
    cart2pol(x, y)

Convert Cartesian coordinates to polar coordinates.
"""
function cart2pol(
    x::Number,
    y::Number
    )
    th = atan(y, x)
    r = hypot(x, y)
    return th, r
end


"""
    pol2cart(th, r)

Convert polar to Cartesian coordinates.
"""
function pol2cart(
    th::Number,
    r::Number
    )
    x = r * cos(th)
    y = r * sin(th)
    return x, y
end


"""
    print_paras(datdir,
        Rphys, outin, zzR,
        N, Nin,
        Rin, RoutC, Rout, din, dout, wvl, z,
        w0, nmax, jlmax, zR, nhigh, ratio,
        spec, r0, l0, L0, sh, jzmax,
        rs, ns, sr2,
        vmean, vstd, tv, deltat, Ts)

Print parameters from `parameters.jl`.
"""
function print_paras(
    outin::Real,
    zzR::Number,
    N::Integer,
    Nin::Integer,
    Rin::Real,
    Rout::Real,
    RoutD::Real,
    din::Real,
    dout::Real,
    wvl::Real,
    z::Real,
    w0::Real,
    nmax::Integer,
    jlmax::Integer,
    zR::Real,
    nhigh::Real,
    ratio::Real,
    spect::String,
    r0::Real,
    l0::Real,
    L0::Real,
    sh::Integer,
    jzmax::Integer,
    rs::Real,
    ns::Integer,
    sr2::Real,
    vmean::Real,
    vstd::Real,
    tc::Real,
    deltat::Real,
    Ts::Integer
    )
    Rphys = 1
    @printf("########## PARAMETER INFO #############\n")

    @printf("\n-- Geometry --\n")
    @printf("outin = %.2f\n", outin)
    @printf("Rin   = %.2fcm\n", Rin * Rphys * 1e2)
    @printf("Rout = %.2fcm\n", Rout * Rphys * 1e2)
    @printf("RoutD  = %.2fcm\n", RoutD * Rphys * 1e2)
    @printf("wvl   = %.0fnm\n", wvl * Rphys * 1e9)
    @printf("zzR   = %.2f\n", zzR)
    @printf("zR    = %.2fm\n", zR * Rphys)
    @printf("z     = %.0fm\n", z * Rphys)
    
    @printf("\n-- Grid --\n")
    @printf("N     = %i\n", N)
    @printf("Nin   = %i\n", Nin)
    @printf("din   = %.2fmm\n", din * Rphys * 1e3)
    @printf("dout  = %.2fmm\n", dout * Rphys * 1e3)

    @printf("\n-- LG basis --\n")
    @printf("w0    = %.2fmm\n", w0 * Rphys * 1e3)
    @printf("nmax  = %d\n", nmax)
    @printf("jlmax = %d\n", jlmax)
    @printf("nhigh = %.2f\n", nhigh)
    @printf("ratio = %.2f\n", ratio)

    @printf("\n-- Turbulence --\n")
    @printf("spec  = %s\n", spect)
    if r0 == Inf
        @printf("r0    = infinity\n")
    else
        @printf("r0    = %.2fmm\n", r0 * Rphys * 1e3)
    end
    @printf("l0    = %.2fmm\n", l0 * Rphys * 1e3)
    if L0 == Inf
        @printf("L0    = infinity\n")
    else
        @printf("L0    = %.2fm\n", L0 * Rphys)
    end
    @printf("sh    = %d\n", sh)
    @printf("jzmax = %d\n", jzmax)
    @printf("rs    = %.d\n", rs)
    @printf("ns    = %i\n", ns)
    @printf("sr2   = %.3f\n", sr2)

    @printf("\n-- Wind --\n")
    @printf("vmean  = %.2fm/s\n", vmean * Rphys)
    @printf("vstd   = %.2fm/s\n", vstd * Rphys)
    @printf("tc     = %.2fms\n", tc * 1e3)
    @printf("deltat = %.2fms\n", deltat * 1e3)
    @printf("Ts     = %d\n", Ts)

    @printf("########################################\n")
end 

"""
    msum(mode)

Total power of a mode.
"""
function msum(mode::Matrix)
    return sum(abs2.(mode))
end

"""
    norm_mode(mode)

Normalize the intensity of a mode such that the total power yields one.
"""
function norm_mode(mode::Matrix)
    return mode / sqrt(sum(abs2.(mode)))
end


"""
    wopt(D, Df)

Theoretical beam width for optimal vacuum propagation.
"""
function wopt(
    D::Number,
    Df::Number
    )
    return sqrt(2) * D / 2 / (1 + 4 * Df)^(1/4)
end

"""
    rayleigh_waist(z, wvl)

Determine beam waist such that `z = zR` (Rayleigh length).
This leads to the smallest possible output aperture.
"""
function rayleigh_waist(
    z::Number,
    wvl::Number
    )
    return sqrt(z * wvl / π)
end

"""
    norm_coeff(vec)

Normalizes the vector vec s.t. `sum(abs.2(vec))=1`.
"""
function norm_coeff(
    vec::Vector
    )
    c = sum(abs2.(vec))
    return vec ./ sqrt(c)
end

end