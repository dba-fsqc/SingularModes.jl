# PhaseScreens.jl
# This module contains all functions for the generation of phase screens.
# Depends on `BasicFunctions.jl` and `ModeGeneration.jl`.
# v5.1

module PhaseScreens

export zk_covmat, ph_ft, ph_sh, ph_zk, ph_hd, ph_kl, ph_x
export phz_x, fourier_shift, wind_shift, get_phzT

using Random
using LinearAlgebra
using Statistics
using HypergeometricFunctions
using Printf

include("BasicFunctions.jl")
include("ModeGeneration.jl")


# I. Zernike covariances

"""
    pFq(alist, blist, z)

Regularized hypergeometric function. Relies on `HypergeometricFunctions.jl`.
"""
function pFq_reg(
    alist::Vector,
    blist::Vector,
    z::Number)
    term = pFq(alist, blist, z)
    return term / prod(gamma.(blist))
end

"""
    zk_cov_atm(j, jp, R; L0=Inf, alpha=5/3)

Analytical Zernike covariances for the atmospheric phase power spectra.
"""
function zk_cov_atm(
    j::Integer,
    jp::Integer,
    R::Number;
    L0::Number=Inf,
    alpha::Number=5/3,
    kwargs...
    )
    n, m = j2nm(j)
    np, mp = j2nm(jp)
    if mp == m && ((j - jp) % 2 == 0 || m == 0)
        prefac1 = (2^(2/3) * gamma(11/6)^2) / (π^2)
        prefac2 = (24/5 * gamma(6/5))^(5/6)
        if L0 == Inf
            term1a = (-1)^((n + np - 2*m)/2)
            term1b = sqrt((n + 1) * (np + 1)) * π
            term1 = term1a * term1b
            term2 = gamma(alpha + 3) * gamma((n + np - alpha)/2)
            term3a = gamma((n - np + alpha + 4)/2)
            term3b = gamma((np - n + alpha + 4)/2)
            term3c = gamma((n + np + alpha + 6)/2)
            term3 = term3a * term3b * term3c
            term4 = (1 / 2)^alpha
            result = prefac1 * prefac2 * term1 * term2 / term3 * term4
            return result * R^alpha 
        else
            kappa0 = 2 * π / L0
            prefac3 =  8π / R^2 * sqrt((n + 1) * (np + 1))
            prefac4 = (-1)^(0.5 * (n + np - m - mp))
            prefactor = prefac1 * prefac2 * prefac3 * prefac4
            term1a = (-1 / (3 * kappa0^(5/3) * gamma(11/6)))
            term1b = 2^(-5. - n - np) * π * R^2 
            term1c = (1 / sin(1 / 6 * (1 + 3n + 3np) * π))
            term1 = term1a * term1b * term1c
            term2a = 11 * 2^(1/3 + n + np) * (R * kappa0)^(5/3)
            term2b = gamma(11/6) * gamma(11/3)
            term2c = pFq_reg(
                [11/6, 7/3, 17/6],
                [
                    1/6 * (11 - 3n - 3np),
                    1/6 * (17 + 3n - 3np),
                    1/6 * (17 - 3n + 3np),
                    1/6 * (23 + 3n + 3np)
                ],
                R^2 * kappa0^2
                )
            term2 = term2a * term2b * term2c
            term3a = 12 * (R * kappa0)^(n + np)
            term3b = gamma(0.5 * (2 + n + np))
            term3c = gamma(3 + n + np)
            term3d = pFq_reg(
                [
                    0.5 * (2 + n + np),
                    0.5 * (3 + n + np),
                    0.5 * (4 + n + np)
                ],
                [
                    2 + n,
                    1/6 * (1 + 3n + 3np),
                    2 + np,
                    3 + n + np
                ],
                R^2 * kappa0^2
                )
            term3 = term3a * term3b * term3c * term3d
            return prefactor * term1 * (term2 - term3)
        end
    else
        return 0 
    end
end

"""
    zk_cov_num(power_spec, j, jp, R; kwargs...)

Numerical Zernike covariances for arbitrary `power_spec`. All `kwargs...`
are passed to `power_spec`.
"""
function zk_cov_num(
    pow_spec::Function,
    j::Integer,
    jp::Integer,
    R::Number;
    kwargs...
    )
    n, m = j2nm(j)
    np, mp = j2nm(jp)
    if mp == m && (j - jp) % 2 == 0 || m == 0
        term1a = 8 * π * sqrt((n + 1) * (np + 1))
        term1b = (-1)^((n + np - m - mp) / 2)
        term1 = term1a * term1b
        term2, _ = quadgk(
            k -> 1 / (R^2 * k) * pow_spec(k; kwargs...) 
            * besselj(n + 1, k * R) 
            * besselj(np + 1, k * R),
            0,
            Inf
            )
        return term1 * term2
    else
        return 0
    end
end

"""
    zk_covmat(
        pow_spec, nmax, R;
        analytics=true, jmax=0, l0=0, alpha=5/3, kwargs...
        )

Zernike covariance matrix obtained analytically for Kolmogorov
and von Karman turbulence , else (or if `analytics=false`) numerically.
"""
function zk_covmat(
    pow_spec::Function,
    nmax::Integer,
    R::Number;
    analytics::Bool=true,
    jmax::Integer=0,
    l0::Number=0,
    alpha::Number=5/3,
    kwargs...
    )
    if jmax != 0
        nmax = j2nmax(jmax)
    end
    jmax = nmax2j(nmax)
    check = string(pow_spec) == "atm_spec" && l0 == 0 && alpha == 5/3
    cov = zeros(jmax - 1, jmax - 1)
        for i in 2:jmax, j in 2:jmax
            if i <= j
                if analytics && check
                    cov[i - 1, j - 1] = zk_cov_atm(
                        i,
                        j,
                        R;
                        alpha=alpha,
                        kwargs...
                        )
                else
                    cov[i - 1, j - 1] = zk_cov_num(
                        pow_spec,
                        i,
                        j,
                        R;
                        alpha=alpha,
                        l0=l0,
                        kwargs...
                        )
                end
            else
                cov[i - 1, j - 1] = cov[j - 1, i - 1]
            end

        end
    return cov
end


# II. Phase screens and arrays

"""
    ph_ft(
        pow_spec, N, R, rs;
        silent=false, kwargs...
        )

Fourier phase screen for a given `pow_spec` on an `N x N` array
with physical radius `R` and radom seed `rs`.
"""
function ph_ft(
    pow_spec::Function,
    N::Integer,
    R::Number,
    rs::Number;
    silent::Bool=false,
    kwargs...
    )
    Random.seed!(rs)
    nmax = 0
    jzmax = 0
    if !silent
        fname = nameof(var"#self#")
        @printf "%s: N=%d, R=%.3f, nmax=%d (jzmax=%d), rs=%d\n" fname N R nmax jzmax rs
    end
    df = 1 / (2*R)
    n2 = floor(Integer, N/2)
    fx = collect(-n2:n2 - 1) * df
    fx, fy = meshgrid(fx, fx)
    f = abs.(fx + 1im * fy)
    ps = pow_spec.(2π * f; kwargs...) * (2π)^2
    ps[floor(Integer, n2 + 1), floor(Integer, n2 + 1)] = 0
    cn = randn(ComplexF64, (N, N)) .* sqrt.(ps) * df
    return inner_circ(real(ift(cn, 1)))
end

"""
    ph_sh(
        pow_spec, N, R, rs;
        sh=5, kwargs...
        )

Subharmonic Fourier phase screen for a given `pow_spec` on an `N x N` array
with physical radius `R` and radom seed `rs`.
Defaults to `sh=5` subharmonic levels.
"""
function ph_sh(
    pow_spec::Function,
    N::Integer,
    R::Number,
    rs::Number;
    sh::Integer=5,
    kwargs...
    )
    Random.seed!(rs)
    fname = nameof(var"#self#")
    @printf "%s: N=%d, R=%.3f, sh=%d, rs=%d\n" fname N R sh rs
    ph1 = ph_ft(pow_spec, N, R, rs, silent=true; kwargs...)
    x0 = collect(-N/2:N/2 - 1) * 2 * R / N
    xx, yy = meshgrid(x0, x0)
    ph2 = zeros(ComplexF64, size(ph1))
    for p in 1:sh
        df = 1 / (3^p * 2*R)
        fx = vec([-1 0 1]) * df
        fx, fy = meshgrid(fx, fx)
        f = abs.(fx + 1im * fy)
        ps = pow_spec.(2*π * f; kwargs...) * (2*π)^2
        ps[2, 2] = 0
        cn = randn(ComplexF64, (3, 3)) .* sqrt.(ps) * df
        SH = zeros(ComplexF64, N, N)
        for j in 1:3, i in 1:3
            SH .+= cn[i, j] * exp.(1im*2π * (fx[i, j] * xx + fy[i, j] * yy))
        end
        ph2 .+= SH
    end
    ph2 = inner_circ(real(ph2) .- mean(real(ph2)))
    return ph1 + ph2
end

"""
    ph_zk(
        pow_spec, N, R, rs;
        nmax=9, jzmax=0, zks=zks, ascorr=ascorr, silent=false, kwargs...
        )

Zernike phase screen for a given `pow_spec` on an `N x N` array
with physical radius `R` and radom seed `rs`.
Defaults to `nmax=9` pyramidal levels.
"""
function ph_zk(
    pow_spec::Function,
    N::Integer,
    R::Number,
    rs::Number;
    nmax::Integer=9,
    jzmax::Integer=0,
    zks::Array=zeros(1,1,1),
    ascorr::Vector=zeros(1),
    silent::Bool=false,
    kwargs...
    )
    if jzmax !== 0
        nmax = j2nmax(jzmax)
    end
    jzmax = nmax2j(nmax)
    if abs(maximum(ascorr)) > 0
        jzmax = size(ascorr, 1) + 1
    else
        ascorr = zeros(jzmax - 1)
    end
    nmax = j2nmax(jzmax)
    
    Random.seed!(rs)
    if !silent
        fname = nameof(var"#self#")
        @printf "%s: N=%d, R=%.3f, nmax=%d (jzmax=%d), rs=%d\n" fname N R nmax jzmax rs
    end
    if zks == zeros(1, 1, 1)
        zks = mode_pyramid_load(zk, nmax, N, R)
    end
    covmat = zk_covmat(pow_spec, nmax, R; kwargs...)
    u, s, _ = svd(covmat, full=true)
    bs = randn(jzmax - 1) .* sqrt.(s)
    as = u * bs
    ph = zeros(N, N)
    for i in 1:jzmax - 1
        ph += (as[i] - ascorr[i]) * zks[i + 1, :, :]
    end
    return ph
end

"""
    zk_overlap(
        ph, zks
        )

Computes the overlap of a phase screen `ph` with the Zernike polynomials `zks`.
"""
function zk_overlap(
    ph::Matrix,
    zks::Array,
    )
    N = size(zks, 2)
    Nz = size(zks, 1)
    as = zeros(Nz - 1)
    for i in 1:Nz-1
        as[i] = sum(ph .* zks[i + 1, :, :]) * 4 / N^2 / π
    end
    return as
end

"""
    ph_hd(
        pow_spec, N, R, rs;
        nmax=9, jzmax=0, zks=zks, phft=phft, nans=false, kwargs...
        )

Hybrid phase screen for a given `pow_spec` on an `N x N` array
with physical radius `R` and radom seed `rs`.
Defaults to `nmax=9` pyramidal levels for the constituent Zernike screen.
"""
function ph_hd(
    pow_spec::Function,
    N::Integer,
    R::Number,
    rs::Number;
    nmax::Integer=9,
    jzmax::Integer=0,
    zks::Array=zeros(1,1,1),
    phft::Matrix=zeros(1,1),
    nans::Bool=false,
    kwargs...)
    Random.seed!(rs)
    rsft, rszk = rand(1:10^6, 2)
    if jzmax !== 0
        nmax = j2nmax(jzmax)
    end
    jzmax = nmax2j(nmax)
    jlimit = nmax2j(1)
    if jzmax < jlimit
        jzmax = jlimit
    end
    nmax = j2nmax(jzmax)
    fname = nameof(var"#self#")
    @printf "%s: N=%d, R=%.3f, nmax=%d (jzmax=%d), rs=%d\n" fname N R nmax jzmax rs
    if zks == zeros(1, 1, 1)
        zks = mode_pyramid_load(zk, nmax, N, R)
    end
    if phft == zeros(1, 1)
        phft = ph_ft(pow_spec, N, R, rsft, silent=true; kwargs...)
    end
    ascorr = zk_overlap(phft, zks)
    phzkmod =  ph_zk(pow_spec, N, R, rszk, zks=zks, ascorr=ascorr, silent=true; kwargs...)
    return inner_circ(phft + phzkmod, nans=nans)
end

"""
    ph_kl(
        pow_spec, N, R, rs;
        nmax=9, jzmax=0, zks=zks, kls=kls, phft=phft, klsReturn=false,
        kwargs...
        )

Hybrid Karhunen-Loève phase screen for a given `pow_spec` on an `N x N` array
with physical radius `R` and radom seed `rs`.
Defaults to `nmax=9` pyramidal levels for the constituent Zernike screen.
"""
function ph_kl(
    pow_spec::Function,
    N::Integer,
    R::Number,
    rs::Number;
    nmax::Integer=9,
    jzmax::Integer=0,
    zks::Array=zeros(1,1,1),
    kls::Array=zeros(1,1,1),
    phft::Matrix=zeros(1,1),
    klsReturn::Bool=false,
    kwargs...)
    Random.seed!(rs)
    rsft, rszk = rand(1:10^6, 2)
    if jzmax !== 0
        nmax = j2nmax(jzmax)
    end
    jzmax = nmax2j(nmax)
    jlimit = nmax2j(1)
    if jzmax < jlimit
        jzmax = jlimit
    end
    nmax = j2nmax(jzmax)
    fname = nameof(var"#self#")
    @printf "%s: N=%d, R=%.3f, nmax=%d (jzmax=%d), rs=%d\n" fname N R nmax jzmax rs
    if zks == zeros(1, 1, 1)
        zks = mode_pyramid_load(zk, nmax, N, R)
    end
    covmat = zk_covmat(pow_spec, nmax, R; kwargs...)
    u, s, _ = svd(covmat, full=true)
    bs = randn(jzmax - 1) .* sqrt.(s)
    if kls == zeros(1, 1, 1)
        kls = zeros(jzmax, N, N)
        for i in 1:jzmax-1
            for j in 1:jzmax-1
                kls[i + 1, :, :] += u[j, i] .* zks[j + 1, :, :]
            end
        end
    end
    if phft == zeros(1, 1)
        phft = ph_ft(pow_spec, N, R, rsft, silent=true; kwargs...)
    end
    bscorr = zk_overlap(phft, kls)
    phmod = zeros(N, N)
    for i in 1:jzmax - 1
        phmod += (bs[i] - bscorr[i]) * kls[i + 1, :, :]
    end
    if klsReturn
        return inner_circ(phft + phmod), kls
    else
        return inner_circ(phft + phmod)
    end
end

"""
    ph_x(
        pow_spec;
        N=1024, R=1, rs=0, kind="hd", kwargs...
        )

Wrapper function for single phase screens. Supported kinds are "ft", "sh",
"zk", "hd", and "kl", see keyword arguments of the corresponding functions.
Defaults to a hybrid phase screen with `N=1024`, `R=1`, and `rs=0`.
"""
function ph_x(
    pow_spec::Function;
    N::Integer=1024,
    R::Number=1,
    rs::Integer=0,
    kind::String="hd",
    kwargs...
    )
    if kind == "hd"
        return ph_hd(pow_spec, N, R, rs; kwargs...)
    elseif kind == "ft"
        return ph_ft(pow_spec, N, R, rs; kwargs...)
    elseif kind == "sh"
        return ph_sh(pow_spec, N, R, rs; kwargs...)
    elseif kind == "zk"
        return ph_zk(pow_spec, N, R, rs; kwargs...)
    elseif kind == "kl"
        return ph_kl(pow_spec, N, R, rs; kwargs...)
    else
        println("Possible kinds are 'zk', 'ft', 'sh', 'hd' (default), or 'kl'.")
    end
end

"""
    ph_x(
        pow_spec, Rarray;
        N=1024, R=1, rs=0, kind="hd", r0tot=0, kwargs...
        )

Wrapper function for phase screen arrays along the z-axis specied via `Rarray`.
Supported kinds are "ft", "sh", "zk", "hd", and "kl",
see keyword arguments of the corresponding functions.
Defaults to a hybrid phase screens with `N=1024`, `R=1`, and `rs=0`.
The total Fried parameter has either to be specified via `r0tot` or multiplied
manually for non-Kolmogorov turbulence.
"""
function phz_x(
    pow_spec::Function,
    Rarray::Vector;
    N::Integer=1024,
    kind::String="hd",
    rs::Integer=0,
    r0tot::Number=0,
    kwargs...
    )
    ns = floor(Integer, (size(Rarray, 1) - 1)/2)
    Random.seed!(rs)
    rss = rand(1:10^6, ns)
    phz = zeros(ns, N, N)
    Threads.@threads for i in 1:ns
        phz[i, :, :] = 1/ns * ph_x( # fix this!!!!!!!!!!
            pow_spec,
            kind=kind,
            N=N,
            R=Rarray[i],
            rs=rss[i];
            kwargs...
            )
    end
    if r0tot == 0
        return phz
    else 
        return phz * r0tot^(-5/3)
    end
end


# III. Phase screen arrays with shifting and cropping

"""
    screen_factor(
        vmean, vstd, tmax, R;
        offset=0.1
        )

Computes the required size of a large, time-shifted phase screen based on
the mean transverse wind velocity `vmean`, its standard deivation `vstd`, the
maximally considered time `tmax`, and the smaller subscreen's radius `R`.
"""
function screen_factor(
    vmean::Number,
    vstd::Number,
    tmax::Number,
    R::Number;
    offset::Number=0.1,
    )
    L = 2*R
    vmax = vmean + 4 * vstd 
    dmax = L + vmax * tmax + 2 * offset * L
    i = 0
    while dmax >= 2^i * L
        i += 1
    end
    return 2^i
end

"""
    fourier_shift(
        a, dshift;
        d=1
        )

Subpixel Fourier shift of a two-dimensional array `a` by distance `dshift`. 
"""
function fourier_shift(
    a::Matrix,
    dshift::Number;
    d::Number=1,
    )
    if dshift == 0
        return a
    else
        N = size(a, 1)
        fa = ft(a, d)
        xi0 = collect(-N/2:N/2-1) * (1 / (N * d))
        fs = fa .* exp.(1im * 2π * dshift * xi0')
        return real(ift(fs, 1 / (N * d)))
    end
end

"""
    wind_shift(
        phzB, Norg, Rarray, vmean, vstd, time;
        rsw=0, offset=0.1, threads=true, kwargs...
        )

Returns a z-array of wind-shifted phase screens with resolution `Norg x Norg`
based on a large underlaying screen array `phzB` with mean transverse wind
velocity `vmean` and its standard deivation `vstd` at time `time`.
The random seed for the wind velocity distribution defaults to `rsw=0`, the
`offset=0.1` defaults to a 10% margin from the large screens boundaries.
"""
function wind_shift(
    phzB::Array,
    Norg::Integer,
    Rarray::Vector,
    vmean::Number,
    vstd::Number,
    time::Number;
    rsw::Integer=0,
    offset::Number=0.1,
    threads::Bool=true,
    kwargs...
    )
    ns = size(phzB, 1)
    N0 = size(phzB, 2)
    Random.seed!(rsw)
    deltax = (randn(ns) * vstd .+ vmean) * time
    phzout = zeros(ns, Norg, Norg)
    offset_pixel = floor(Integer, Norg * offset)
    function to_loop(i::Integer)
        if deltax[i] !== 0
            di = N / (2 * Rarray[2*i + 1])
            phshift = fourier_shift(
                phzB[i, :, :],
                deltax[i] * di
                )
        else
            phshift = phzB[i, :, :]
        end
        idx1 = floor(Integer, N0 / 2) - floor(Integer, Norg / 2) + 1
        phzout[i, :, :] = inner_circ(
            phshift[
                idx1:idx1 + Norg - 1,
                offset_pixel:offset_pixel + Norg - 1
                ]
            )
    end
    if threads
        Threads.@threads for i in 1:ns
            to_loop(i)
        end
    else
        for i in 1:ns
            to_loop(i)
        end
    end
    return phzout
end

"""
    get_phzT(
        phzB, Rarray, Tarray, vmean, vstd;
        Norg=1024, rsw=0, kwargs...
        )

Returns a time array of z-arrays of wind-shifted phase screens
based on a large underlaying screen z-array `phzB` with mean transverse wind
velocity `vmean` and its standard deivation `vstd` at time `time`.
The random seed for the wind velocity distribution defaults to `rsw=0`, and
the subscreen side length to `Norg=1024` pixels.
"""
function get_phzT(
    phzB::Array,
    Rarray::Vector,
    Tarray::Vector,
    vmean::Real,
    vstd::Real;
    Norg::Integer=1024,
    rsw::Integer=0,
    kwargs...
    )
    ns = size(phzB, 1)
    Ts = size(Tarray, 1)
    phzT = zeros(Ts, ns, Norg, Norg)
    Threads.@threads for i in 1:Ts
        phzT[i, :, : ,:] = wind_shift(
            phzB,
            Norg,
            Rarray,
            vmean,
            vstd,
            Tarray[i],
            rsw=rsw,
            threads=false
            )
    end
    return phzT
end

end