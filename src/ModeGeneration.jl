# ModeGeneration.jl
# This module contains all functions related to optical mode generation.
# Depends on `BasicFunctions.jl`.
# v2.2

module ModeGeneration


export pl2j, nm2j, mn2j, j2pl, j2nm, j2mn, j2nmax, nmax2j, noll_pyramid
export zk, lg, hg, besselbeam, besselgaussian, points2grid
export lg_mode, hg_mode, bg_mode, zk_mode, mode_pyramid, zka_load
export mode_pyramid_load
export pcl, pcl_mod, pcl_arc, pcl_geo, pcl_lim, bcl
export rho_mean_bg, pcl_bg

using SpecialFunctions
using Polynomials
using SpecialPolynomials
using HypergeometricFunctions
using JLD2
using Printf

using ..BasicFunctions


# I. Index conventions

"""
    pl2j(p, l; zk=false)

Double to single index conversion on basis of Noll's convention. Use `p, l`
for LG modes, or `n, m` for ZK modes (with `zk=true`) or `m, n` for HG modes
(with `hg=true`).
"""
function pl2j(
    p::Integer,
    l::Integer;
    zk::Bool=false,
    hg::Bool=false
    )
    if hg
        zk = true
        (p, l) = (p + l, p - l)
    end
    if ! zk
        n = 2 * p + abs(l)
    else
        n = p
    end
    out = n * (n + 1) / 2 + abs(l)
    if l >= 0
        if (n % 4 == 2) || (n % 4 == 3)
            out += 1
        end
    end
    if l <=0
        if (n % 4 == 0) || (n % 4 == 1)
            out += 1
        end
    end
    return trunc(Integer, out)
end

"""
    nm2j(n, m) = pl2j(n, m, zk=true)

Shorthand for ZK Noll indices. See also [`pl2j`](@ref).
"""
function nm2j(
    n::Integer,
    m::Integer
    )
    pl2j(n, m, zk=true)
end

"""
    mn2j(m, n) = pl2j(m, n, hg=true)

Shorthand for HG Noll indices. See also [`pl2j`](@ref).
"""
function mn2j(
    m::Integer,
    n::Integer
    )
    pl2j(m, n, hg=true)
end

"""
    j2pl(j, zk=false)

Reconversion of Noll superindex `j` to individual indices `p, l` for LG modes
and `n, m` for ZK modes (with `zk=true`).
"""
function j2pl(
    j::Integer;
    zk::Bool=false,
    hg::Bool=false
    )
    if hg
        zk = true
    end
    n = trunc(Integer, (-1 + sqrt(1 + 8 * (j - 1))) / 2)
    if n < 2
        l = n
    elseif n > 1 && n % 2 == 0
        rj = j - trunc(Integer, n * (n + 1) / 2)
        l = trunc(Integer, rj / 2) * 2
    elseif n > 1 && n % 2 != 0
        rj = j - trunc(Integer, n * (n + 1) / 2)
        l = trunc(Integer, (rj + 1) / 2) * 2 - 1
    end
    if j % 2 != 0
        l = -l
    end
    if ! zk
        p = trunc(Integer, (n - abs(l)) / 2)
        return p, l
    elseif hg
        return (trunc(Integer, (n + l) / 2), trunc(Integer, (n - l) / 2))
    else
        return n, l
    end
end

"""
    j2nm(j) = j2pl(j, zk=true)

Shorthand for reconversion of ZK Noll indices. See also [`j2pl`](@ref).
"""
function j2nm(
    j::Integer
    )
    j2pl(j, zk=true)
end

"""
    j2mn(j) = j2mn(j, hg=true)

Shorthand for reconversion of HG Noll indices. See also [`j2pl`](@ref).
"""
function j2mn(
    j::Integer
    )
    j2pl(j, hg=true)
end

"""
    j2nmax(j)

Maximal radial order `nmax` for a given superindex `j`.
See also [`j2pl`](@ref).
"""
function j2nmax(
    j::Integer
    )
    n, _ = j2pl(j, zk=true)
    return n
end

"""
    nmax2j(j)

Maximal superindex `j` for a given radial order `nmax`.
See also [`j2pl`](@ref).
"""
function nmax2j(
    nmax::Integer
    )
    return maximum([pl2j(0, -nmax), pl2j(0, nmax)])
end


# II. Modes

"""
    noll_pyramid(nmax)

Print a pyramid of Noll indices up to the radial order `nmax`.
See also [`j2pl`](@ref), [`pl2j`](@ref).
"""
function noll_pyramid(
    nmax::Integer;
    tuples::Bool=false,
    tuple2j::Function=nm2j,
    j2tuple::Function=j2nm
    )
    if ! tuples
        maxi = tuple2j(nmax, nmax)
        maxi2 = tuple2j(nmax, -nmax)
        if maxi2 > maxi
            maxi = maxi2
        end
        ms = length(string(maxi))
    else
        ms1 = length(string(tuple2j(nmax, nmax))*"=($nmax,$nmax)")
        ms2 = length(string(tuple2j(nmax, -nmax))*"=($nmax,-$nmax)")
        ms = maximum([ms1, ms2])
    end
    mss = ""
    for i in 1:ms
        mss *= " "
    end
    for n in 0:nmax
        line = ""
        for m in -nmax:nmax
            if abs(m) > n
                line *= mss
            elseif (n % 2 == 0 && m % 2 == 0) || (n % 2 != 0 && m % 2 != 0)
                (an, am) = j2tuple(nm2j(n,m))
                if tuples
                    line *= lpad(string(tuple2j(an, am))*"=($an,$am)", ms, " ")
                else
                    line *= lpad(tuple2j(an, am), ms, " ")
                end
            else
                line *= mss
            end
        end
        println(line)
    end
end

"""
    r_z(r, n, m)

Radial function of the Zernike polynomials based on [Noll (1976)], Eq. (1).
"""
function r_z(
    r::Number,
    n::Integer,
    m::Integer
    )
    val = 0
    max_sum = trunc(Integer, (n - abs(m)) / 2) + 1
    for s in 0:max_sum-1
        # Int128 works up to factorial 33 (= number of radial levels)
        term1 = (-1)^s * factorial(Int128(n - s)) * r^(n - 2s)
        term2 = factorial(Int128(s))
        term3 = factorial(Int128(trunc(Integer, (n + m) / 2) - s))
        term4 = factorial(Int128(trunc(Integer, (n - m) / 2) - s))
        val += term1 / (term2 * term3 * term4)
    end
    return val
end

"""
    zk(x, y; R=1, j=1, angle=0)

Zernike polynomial with Noll index `j`, radius `R` and rotation angle `angle`.
See also [`r_z`](@ref).
"""
function zk(
    x::Number,
    y::Number;
    R::Number=1,
    j::Integer=1,
    angle::Number=0,
    d::Number=1,
    kwargs...
    )
    n, m = j2nm(j)
    th, r = cart2pol(x, y)
    if r <= R
        if m == 0
            out = sqrt(n + 1) * r_z(r/R, n, 0)
        elseif m != 0 && j % 2 == 0
            if m >= 0
                out = sqrt(n + 1) * r_z(r/R, n, m) * cos(m*th - angle) * sqrt(2)
            else
                out = sqrt(n + 1) * r_z(r/R, n, m) * cos(m*th + angle) * sqrt(2)
            end
        elseif m != 0 && j % 2 != 0
            if m >= 0
                out = sqrt(n + 1) * r_z(r/R, n, m) * sin(m*th - angle) * sqrt(2)
            else
                out = sqrt(n + 1) * r_z(r/R, n, m) * sin(m*th + angle) * sqrt(2)
            end
        end
        return out * d
    else
        return 0
    end
end

"""
    laguerre(n, alpha, x)

Associated/generalized Laguerre polynomial with radial order `n` and
azimuthal index `alpha`. Definition from wikipedia.
"""
function laguerre(
    n::Integer,
    alpha::Number,
    x::Number
    )
    p0, p1 = 1, 1 + alpha - x
    n == 0 && return p0
    for k = 1:n-1
        p1, p0 = ((2*k + 1 + alpha -x)*p1 - (k + alpha)*p0)/(k + 1), p1
    end
    p1
end

"""
    lg(x, y; w0=1, z=0, wvl=1550e-9, p=0, l=0, d=1, j=0)

Pointwise Lagurre-Gaussian mode function. See also [`laguerre`](@ref).
"""
function lg(
    x::Number,
    y::Number;
    w0::Number=1,
    z::Number=0,
    wvl::Number=1550e-9,
    p::Integer=0,
    l::Integer=0,
    d::Number=1,
    j::Integer=0,
    noabs::Bool=false,
    fullbeam::Bool=false,
    kwargs...
    )
    if j !== 0
        p, l = j2pl(j)
    end
    ph, r = cart2pol(x, y)  # from BasicFunctions
    al = abs(l)
    if noabs
        al = l
    end
    if z == 0
        L = laguerre(p, al, 2 * r^2 / w0^2)
        u = (sqrt(2) * r / w0)^al * L * exp(-r^2 / w0^2) * exp(1im * l * ph)
        c = convert(
            Float64,
            sqrt(abs(2 * factorial(BigInt(p)) / (π * factorial(BigInt(p + al)))))
            )
        u = u * c / w0 * d
        return u
    else
        z_R = π * w0^2 / wvl
        w0 = w0 * sqrt(1 + (z / z_R)^2)
        c = convert(
            Float64,
            sqrt(abs(2 * factorial(BigInt(p)) / (π * factorial(BigInt(p + al)))))
            )
        term1 = (sqrt(2) * r / w0)^al
        L = laguerre(p, al, 2 * r.^2 / w0^2)
        term2 = exp(-r.^2 / w0^2)
        term3 = exp(1im * 2 * π / wvl * r^2 * z / 2 / (z^2 + z_R^2))
        term4 = exp(1im * l * ph)
        term5 = exp(-1im * (2 * p + al + 1) * atan(z / z_R))
        u = c / w0 * term1 * L * term2 * term3 * term4 * term5 * d
        if fullbeam
            u = u * exp(1im * 2 * π / wvl * z)
        end
        return u
    end
end

"""
    hg_1d(x, y; w0=1, z=0, wvl=1550e-9, m=0, d=1)

Pointwise one-dimensional Hermite-Gaussian mode function.
"""
function hg_1d(
    x::Number;
    w0::Number=1,
    z::Number=0,
    wvl::Number=1550e-9,
    m::Integer=0,
    d::Number=1
    )
    if z == 0
        term1 = 1 / sqrt(2^(m - 1/2) * factorial(big(m)) * sqrt(π)) / sqrt(w0)
        term2 = evalpoly(sqrt(2) * x / w0, basis(Hermite, m)) * exp(-x^2 / w0^2) * d
        return Complex{Float64}(term1 * term2)
    else
        z_R = π * w0^2 / wvl
        k = 2 * π / wvl
        w0 = w0 * sqrt(1 + (z / z_R)^2)
        Rz = z * (1 + (zR / z)^2)
        term1 = 1 / sqrt(2^(m - 1/2) * factorial(big(m)) * sqrt(π)) / sqrt(w0)
        term2 = evalpoly(sqrt(2) * x / w0, basis(Hermite, m)) * exp(-x^2 / w0^2) * d
        term3 = exp(1im * (k*z - (m + 1/2) * atan(z / z_R) + k * x^2 / (2 * Rz)))
        return Complex{Float64}(term1 * term2 * term3)
    end
end

"""
    hg(x, y; w0x=1, w0y=1, w0=1, z=0, wvl=1550e-9, m=0, n=0, d=1)

Pointwise Hermite-Gaussian mode function.
"""
function hg(
    x::Number,
    y::Number;
    w0x::Number=1,
    w0y::Number=1,
    w0::Number=1,
    z::Number=0,
    wvl::Number=1550e-9,
    m::Integer=0,
    n::Integer=0,
    d::Number=1,
    j::Integer=0,
    kwargs...
    )
    if j != 0
        n, m = j2mn(j)
    end
    if w0x == w0y
        w0x = w0
        w0y = w0
    end
    return hg_1d(x, z=z, w0=w0x, m=m, d=sqrt(d)) * hg_1d(y, z=z, w0=w0y, m=n, d=sqrt(d))
end

"""
    besselbeam(x, y; z=0, wvl=1550e-9, beta=0, l=0, paraxial=true)

Pointwise Bessel mode function. For non-paraxial wave equation set `paraxial=false`.
"""
function besselbeam(
    x::Number,
    y::Number;
    z::Number=0,
    wvl::Number=1550e-9,
    beta::Number=0,
    l::Number=0,
    paraxial::Bool=true
    )
    k = 2*π/wvl
    ph, r = cart2pol(x, y)
    if paraxial
        kz = k - beta^2/(2k)
    else
        kz = sqrt(k^2 - beta^2)
    end
    return besseljx.(l, beta*r) * exp(1im*l*ph) * exp(1im * kz * z) 
end

"""
    besselgaussian(x, y; w0=1, z=0, wvl=1550e-9, beta=0, l=0)

Pointwise Bessel-Gaussian mode function. See also [`besselbeam`](@ref).
"""
function besselgaussian(
    x::Number,
    y::Number;
    w0::Number=1,
    z::Number=0,
    wvl::Number=1550e-9,
    beta::Number=0,
    l::Number=0
    )
    k = 2*π/wvl
    zR = k * w0^2 / 2
    ph, r = cart2pol(x, y)
    al = abs(l)
    term0a = 1 / sqrt(2 * π) * 2 * exp(beta^2 * w0^2 / 8)
    term0b = w0 * sqrt(besselix(al, beta^2 * w0^2 / 4))
    if z != 0
        wz = w0 * sqrt(1 + z^2 / zR^2)
        zeta = atan(z / zR)
        Rz = z * (1 + zR^2 / z^2)
        term1 = (-1)^al * w0 / wz * exp(1im*l*ph)
        term2 = exp(-1im * (beta^2 * z / 2 / k + zeta))
        term3 = besseljx.(al, beta * r * zR / (zR + 1im * z))
        term4 = exp((-1 / wz^2 + 1im * k / 2 / Rz) * (r^2 + beta^2 * z^2 / k^2))
        return term0a / term0b * term1 * term2 * term3 * term4
    else
        term1 = (-1)^al * exp(1im*l*ph)
        term2 = besseljx.(al, beta * r)
        term3 = exp(-r^2 / w0^2 )
        return term0a / term0b * term1 * term2 * term3
    end
end


# III. Mode arrays

"""
    make_grid(N, R)

Meshgrid of size `N x N` with sidelength `2R`.
"""
function make_grid(
    N::Integer=1024,
    R::Number=1 
    )
    x = -N/2:N/2 - 1
    d = 2*R/N
    xx = x' .* ones(N) .* d
    yy = ones(N)' .* x .* d
    return xx, yy
end


"""
    _default_N(N, fpointsname)

Default gridsize `N` for LG and ZK modes.
See also [`lg`](@ref), [`zk`](@ref), [`points2grid`](@ref).
"""
function _default_N(
    N::Integer=0,
    fpointsname::String="lg"
    )
    if N == 0
        if fpointsname == "lg"
            return 256
        elseif fpointsname == "hg"
            return 256
        elseif fpointsname == "besselgaussian"
            return 256
        else
            return 1024
        end
    else
        return N
    end
end

"""
    points2grid(fpoints; N=1024, R=1, kwargs...)

Evaluate the function `fpoints` on an `N x N` grid of sidelength `2R`
(default `N=1024`). Keywords arguments `kwargs...` are passed to `fpoints`.
See also [`_default_N`](@ref).
"""
function points2grid(
    fpoints::Function;
    N::Integer=0,
    R::Number=1,
    kwargs...
    )
    N = _default_N(N, string(fpoints))
    xx, yy = make_grid(N, R)
    fpoints.(xx, yy; kwargs...) 
end

"""
    lg_mode(; N=256, R=1, kwargs...)

Shorthand for a single LG mode on a grid.
See also [`lg`](@ref), [`points2grid`](@ref).
"""
function lg_mode(
    ;
    N::Integer=0,
    R::Number=1,
    kwargs...
    )
    points2grid(lg; N, R, kwargs...)
end

"""
    hg_mode(; N=256, R=1, kwargs...)

Shorthand for a single HG mode on a grid.
See also [`hg`](@ref), [`points2grid`](@ref).
"""
function hg_mode(
    ;
    N::Integer=0,
    R::Number=1,
    kwargs...
    )
    points2grid(hg; N, R, kwargs...)
end

"""
    bg_mode(; N=256, R=1, kwargs...)

Shorthand for a single BG mode on a grid.
See also [`besselgaussian`](@ref), [`points2grid`](@ref).
"""
function bg_mode(
    ;
    N::Integer=0,
    R::Number=1,
    kwargs...
    )
    points2grid(besselgaussian; N, R, kwargs...)
end

"""
    zk_mode(; N=1024, R=1, kwargs...)

Shorthand for a single ZK mode on a grid.
See also [`zk`](@ref), [`points2grid`](@ref).
"""
function zk_mode(
    ;
    N::Integer=0,
    R::Number=1,
    kwargs...
    )
    points2grid(zk; N, R, kwargs...)
end

"""
    mode_pyramid(fpoints, nmax, N, R; jmax=0, threads=true, kwargs...)

Generate a pyramid array of a given mode via `fpoints`, i.e.,
`kwargs...` -> `fpoints`. See also [`lg`](@ref), [`zk`](@ref). 
"""
function mode_pyramid(
    fpoints::Function,
    nmax::Integer,
    N::Integer,
    R::Number;
    jmax::Integer=0,
    threads::Bool=true,
    kwargs...
    )
    N = _default_N(N, string(fpoints))  # only acts if N=0
    if jmax !== 0
        nmax = j2nmax(jmax)
    end
    jmax = nmax2j(nmax)
    xx, yy = make_grid(N, R)
    modes = zeros(typeof((fpoints(0, 0))), jmax, N, N)
    if threads
        Threads.@threads for k in 1:jmax
            modes[k, :, :] = fpoints.(xx, yy; R=R, j=k, d=2*R/N, kwargs...)
        end
    else
        for k in 1:jmax
            modes[k, :, :] = fpoints.(xx, yy; R=R, j=k, d=2*R/N, kwargs...)
        end
    end
    return modes
end

"""
    mode_pyramid_load(fpoints, nmax, N, R; kwargs...)

Load the mode pyramid from file (currently disabled).
See also [`mode_pyramid`](@ref).
"""
function mode_pyramid_load(
    fpoints::Function,
    nmax::Integer,
    N::Integer,
    R::Number;
    w0::Number=1,
    kwargs...
    )
    # currently disabled to avoid double naming
    modes = mode_pyramid(fpoints, nmax, N, R, w0=w0; kwargs...)
    # println("w0=$w0")
    # N = _default_N(N, string(fpoints))
    # savename = @sprintf "%ss-%d-%.3f-%d-%.3f.jld" string(fpoints) N R nmax w0
    # savepath = "static/$savename"
    # if isfile(savepath)
    #     @load savepath modes
    #     println("$savename loaded")
    # else
    #     mkpath("static")
    #     modes = mode_pyramid(fpoints, nmax, N, R, w0=w0; kwargs...)
    #     @save savepath modes
    #     println("$savename generated and saved")
    # end
    return modes
end

"""
    zka_load(nmax, darray; N=1024)

Generate an array of ZK pyramids for the multi-step propagation.
Load from file if it exsists.
"""
function zka_load(
    nmax::Integer,
    darray::Vector;
    N::Integer=1024
    )
    ns = floor(Integer, (size(darray, 1) - 1)/2)
    savename = @sprintf "zka-%d-%d-%.2f-%d.jld" N ns R nmax
    savepath = "static/$savename"
    if isfile(savepath)
        @load savepath zka
        println("$savename loaded")
    else
        mkpath("static")
        jmax = nmax2j(nmax)
        zka = zeros(Float64, ns, jmax, N, N)
        Threads.@threads for i in 1:ns
            R = N*darray[2*i]
            println(R)
            zka[i, :, :, :] = mode_pyramid(zk, nmax, N, R; threads=false)
        end
        @save savepath zka
        println("$savename generated and saved")
    end
    return zka
end


# IV. Phase correlation length

"""
    pcl(l; w0=1)

Phase correlation length based on Eq.~(11) in [Leonhard (2015)].
"""
function pcl(
    l::Integer;
    w0::Number=1
    )
    term1 = sin(π / (2 * abs(l))) * w0 / sqrt(2) 
    term2 = gamma(abs(l) + 3 / 2) / gamma(abs(l) + 1)
    return term1 * term2
end

"""
    pcl_mod(l; w0=1)

Modified phase correlation length based on an isosceles triangle.
"""
function pcl_mod(
    l::Integer;
    w0::Number=1
    )
    term1 = 2 * sin(π / (4 * abs(l))) * w0 / sqrt(2) 
    term2 = gamma(abs(l) + 3 / 2) / gamma(abs(l) + 1)
    return term1 * term2
end

"""
    pcl_arc(l; w0=1)

Phase correlation arc.
"""
function pcl_arc(
    l::Integer;
    w0::Number=1
    )
    term1 = π / (2 * abs(l)) * w0 / sqrt(2) 
    term2 = gamma(abs(l) + 3 / 2) / gamma(abs(l) + 1)
    return term1 * term2
end


"""
    pcl_geo(l; w0=1)

Geometric phase correlation length.
"""
function pcl_geo(
    l::Integer;
    w0::Number=1
    )
    return 2*w0 / sqrt(2) * sqrt(abs(l) + 1) * sin(π/(4*abs(l)))
end


"""
    pcl_lim(l; w0=1)

Large l limit of the phase correlation length.
"""
function pcl_lim(
    l::Integer;
    w0::Number=1
    )
    return π * w0/(2 * sqrt(2* abs(l)))
end

"""
    bcl(uin)

Fourier-based beam correlation length of an input beam profile `uin``.
"""
function bcl(
    uin::Matrix
    )
    keff2 = 0
    N = size(uin)[1]
    uft = ft(uin, 1/N)
    for nx in 1:N
        for ny in 1:N
            keff2 += ((nx - (N+1)/2)^2 + (ny - (N+1)/2)^2) * abs2(uft[nx, ny])
        end
    end
    return 1/sqrt(keff2)/4
end

"""
   rho_mean_bg(l, beta; w0=1)

Mean radius of a Bessel-Gaussian beam.
"""
function rho_mean_bg(
    l::Integer,
    beta::Number;
    w0::Number=1
    )
    z1 = -0.5 * beta^2 * w0^2
    z2 = (beta^4 * w0^4) / 64
    l = abs(l)
    num = exp((beta^2 * w0^2) / 4) * w0 *
          gamma(1.5 + l) *
          pFq([0.5 + l, 1.5 + l], [1.0 + l, 1.0 + 2l], z1)
    denom = sqrt(2) * gamma(1.0 + l) *
            pFq([], [1.0 + l], z2)  # _0F_1 is pFq([], [b], z)
    return num / denom
end

"""
   pcl_bg(l, beta; w0=1, mod=false)

Phase correlation length for Bessel-Gaussian beams, as used in
[Bachmann (2024)].
"""
function pcl_bg(
    l::Integer,
    beta::Number;
    w0::Number=1,
    mod::Bool=false
    )
    rmean = rho_mean_bg(l, beta, w0=w0)
    if !mod
        return sin(π / (2 * abs(l))) * rmean
    else
        return 2 * sin(π / (4 * abs(l))) * rmean
    end
end

end