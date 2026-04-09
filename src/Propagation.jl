# Propagation.jl
# This module contains all functions related to the numerical propagation
# of light.
# Depends on `BasicFunctions.jl` and `ModeGeneration.jl`.
# Propagation through random media requires `PhaseScreens.jl`.
# v3.1

module Propagation

export fraunhofer_prop, lens_against_ft, ang_spec_prop
export rayleigh_length, sigma_R2, min_ns, multi_step_variable
export propagate, modeset_propagation, transmission_square
export svd_lgs_square, svd_square_wrapper, transmission_pixel
export svd_lgs_pixel

using LinearAlgebra
using ..BasicFunctions
using ..ModeGeneration


# I. Propagators

"""
    fraunhofer_prop(u1, wvl, d1, z)

Fraunhofer propagation based on Schmidt (2010).
"""
function fraunhofer_prop(
    u1::Matrix,
    wvl::Number,
    d1::Number,
    z::Number
    )
    N = size(u1, 1)
    k = 2 * π / wvl
    x = (-N/2:N/2 - 1) ./ (N * d1) .* wvl .* z
    xx = x' .* ones(N)
    yy = ones(N)' .* x
    term1 = exp.(1im * k / (2 * z) .* (xx.^2 + yy.^2))
    term2 = (1im * wvl * z) .* ft(u1, d1)
    return term1 ./ term2
end

"""
    lens_against_ft(u1, wvl, d1, f)

Propagation from the pupil plane to the focal plane of a lens
for an object placed against (or just before) the lens.
"""
function lens_against_ft(
    u1::Matrix,
    wvl::Number,
    d1::Number,
    f::Number
    )
    N = size(u1, 1)
    k = 2 * π / wvl
    x = (-N/2:N/2 - 1) ./ (N * d1) .* wvl .* f
    xx = x' .* ones(N)
    yy = ones(N)' .* x
    term1 = exp.(-1im * k / (2 * f) .* (xx.^2 + yy.^2))
    term2 = (1im * wvl * f) .* ft(u1, d1)
    return term1 ./ term2
end

"""
    ang_spec_prop(u1, wvl, d1, d2, z)

Fresnel angular spectrum propagation based on [Schmidt (2010)].
"""
function ang_spec_prop(
    u1::Matrix,
    wvl::Number,
    d1::Number,
    d2::Number,
    z::Number
    )
    N = size(u1, 1)
    k = 2 * π / wvl
    x = -N/2:N/2 - 1
    xx = x' .* ones(N)
    yy = ones(N)' .* x
    r1sq = (xx.^2 .+ yy.^2) .* d1^2
    df1 = 1/(d1 * N)
    fsq = (xx.^2 .+ yy.^2) .* df1^2
    m = d2/d1
    r2sq = (xx.^2 .+ yy.^2) .* d2^2
    q1 = exp.(1im * k / 2 * (1-m)/z .* r1sq)
    q2 = exp.(-1im * π^2 * 2 * z / m / k .* fsq)
    q3 = exp.(1im * k / 2 * (m-1)/(m * z) .* r2sq)
    uout = q3 .* ift(q2 .* ft(q1 .* u1 ./ m, d1), df1)
    return uout .* d2/d1
end


# II. Help functions

"""
    rayleigh_length(w0, wvl)

Rayleigh length of optical propagation.
"""
function rayleigh_length(
    w0::Number,
    wvl::Number
    )
    return π * w0^2 / wvl
end

"""
    sigma_R2(wvl, z, r0)

Rytov variance.
"""
function sigma_R2(
    wvl::Number,
    z::Number,
    r0::Number
    )
    k = 2 * π / wvl
    return 2.90 * k^(-5 / 6) * z^(5 / 6) * r0^(-5 / 3)
end

"""
    sigma_R2_1(wvl, z, r0, ns)

Fist (and largest) term of Rytov variance for multi step propagation
as given in Eq.~(5.72) in [Sorelli (2019)]. 
"""
function sigma_R2_1(
    wvl::Number,
    z::Number,
    r0::Number,
    ns::Number
    )
    k = 2 * π / wvl
    term1 = 5.32 * k^(-5 / 6) * z^(5 / 6) * r0^(-5 / 3)
    term2 = (1 - 1 / (2 * ns))^(5 / 6) / ns
    return term1 * term2
end

"""
    cn2(r0, wvl, z)

Refractive index structure constant C_n^2.
"""
function cn2(
    r0::Number,
    wvl::Number,
    z::Number
    )
    k = 2 * π / wvl
    return r0^(-5 / 3) / (0.423 * k^2 * z)
end

"""
    cn2sr(sr2, wvl, z)

Refractive index structure constant C_n^2 via Rytov variance `sr2`.
"""
function cn2sr(
    sr2::Number,
    wvl::Number,
    z::Number
    )
    k = 2 * π / wvl
    return sr2 / (1.23 * k^(7 / 6) * z^(11 / 6))
end

"""
    r0step(r0, ns)

Fried parameter of individual propagation steps (assumes equal step length),
based on Eq.~(5.68) in [Sorelli (2019)].
"""
function r0step(
    r0::Number,
    ns::Number
    )
    return r0 * ns^(3 / 5)
end

"""
    min_ns(z, wvl, r0, limit=0.5)

Determine the minimal number of propagation steps while complying with
the weak scintillaiton condition (default `sr2 < limit = 0.5`).
"""
function min_ns(
    z::Number,
    wvl::Number,
    r0::Number,
    limit::Number=0.5
    )
    ns = 0
    sr2 = 2
    while sr2 > limit
        ns += 1
        sr2 = sigma_R2_1(wvl, z, r0, ns)
    end
    return ns
end


# III. Multi step propagation

"""
    multi_step_variable(uin, wvl, darray, zarray, phz; track=false)

Multi step propagation along a channel described by `darray`, `zarray`
and `phz`.
"""
function multi_step_variable(
    uin::Matrix,
    wvl::Number,
    darray::Vector,
    zarray::Vector,
    phz::Array;
    track::Bool=false
    )
    ns = size(phz, 1)  # number of whole propagation steps
    if track
        N = size(uin, 1)
        usteps = zeros(eltype(uin), ns + 1, N, N)
        usteps[1, :, :] = uin
    end
    u0 = ang_spec_prop(
        uin,
        wvl,
        darray[1],
        darray[2],
        zarray[2] - zarray[1]
        ) 
    for i in 1:ns-1
        u1 = u0 .* exp.(1im .* phz[i, :, :])       
        u0 = ang_spec_prop(
            u1,
            wvl,
            darray[2*i],
            darray[2*i + 2],
            zarray[2*i + 2] - zarray[2*i]
            )
        if track
            usteps[i + 1, :, :] = u0
        end
    end 
    u1 = u0 .* exp.(1im .* phz[ns, :, :])  # last phase screen
    # last half step propagation
    uout = ang_spec_prop(
        u1,
        wvl,
        darray[2*ns],
        darray[2*ns + 1],
        zarray[2*ns + 1] - zarray[2*ns]
        )
    if track
        usteps[ns + 1, :, :] = inner_circ(uout)
        return usteps
    else
        return inner_circ(uout)
    end
end

"""
    propagate(mode, wvl, Rarray, zarray, phz)

Propagate through a channel (including numerical padding and cropping).
See also [`multi_step_variable`](@ref).
"""
function propagate(
    mode::Matrix,
    wvl::Number,
    Rarray::Vector,
    zarray::Vector,
    phz::Array
    )
    N = size(phz, 2)
    Nin = size(mode, 1)
    darray = 2 .* Rarray / Nin
    idxout1 = floor(Int, (N - Nin)/2)
    idxout2 = floor(Int, (N + Nin)/2) - 1
    uin = zeropad(inner_circ(mode), N)
    uout = multi_step_variable(uin, wvl, darray, zarray, phz)
    return inner_circ(uout[idxout1:idxout2, idxout1:idxout2])
end

"""
    modeset_propagation(modes, wvl, Rarray, zarray, phz)

Propagate a whole set of modes through a channel.
See also [`propagate`](@ref).
"""
function modeset_propagation(
    modes::Array,
    wvl::Number,
    Rarray::Vector,
    zarray::Vector,
    phz::Array;
    threads::Bool=true
    )
    jlmax = size(modes, 1)
    Nin = size(modes, 2)
    modes_prop = zeros(eltype(modes), jlmax, Nin, Nin)
    if threads
        Threads.@threads for i in 1:jlmax
            modes_prop[i, :, :] = propagate(
                modes[i, :, :],
                wvl,
                Rarray,
                zarray,
                phz
                )
            println("set_propagation: $i of $jlmax")
        end
    else
        for i in 1:jlmax
            modes_prop[i, :, :] = propagate(
                modes[i, :, :],
                wvl,
                Rarray,
                zarray,
                phz
                )
            println("set_propagation: $i of $jlmax")
        end
    end
    return modes_prop
end


# IV. Transmission operators and SVD

"""
    transmission_square(lgs_vac, lgs_turb)

Square transmission matrix from vacuum propagated (`lgs_vac`) and
turbulence propagated (`lgs_turb`) modes, i.e. matrix with entries
<lgs_vac | lgs_turb>.
"""
function transmission_square(
    lgs_vac::Array,
    lgs_turb::Array
    )
    jlmax = size(lgs_vac, 1)
    Nin = size(lgs_vac, 2)
    lgs_vac_vec = reshape(lgs_vac, jlmax, Nin^2)
    lgs_turb_vec = reshape(lgs_turb, jlmax, Nin^2)
    return conj(lgs_vac_vec * adjoint(lgs_turb_vec))
end

"""
    svd_lgs_square(lgs, lgs_turb, tmat)

Singular modes via SVD of a square transmission matrix `tmat`.
"""
function svd_lgs_square(
    lgs::Array,
    lgs_turb::Array,
    tmat::Matrix;
    return_coeffs::Bool=false
    )
    _, s, v = svd(tmat, full=false)
    jlmax = size(lgs, 1)
    Nin = size(lgs, 2)
    vins = zeros(eltype(lgs), jlmax, Nin, Nin)
    uouts = zeros(eltype(lgs), jlmax, Nin, Nin)
    for i in 1:jlmax
        for j in 1:jlmax
            vins[i, :, :] += v[j, i] .* lgs[j, :, :]
            uouts[i, :, :] += v[j, i] .* lgs_turb[j, :, :]
        end
    end
    if return_coeffs
        return uouts, s, vins, v
    else
        return uouts, s, vins
    end
end

"""
    svd_square_wrapper(Nin, R, w0, nmax, Rarray, zarray, phz)

Wrapper for all steps to obtain the singular modes of a square
transmission matrix.
"""
function svd_square_wrapper(
    Nin::Int,
    w0::Number,
    wvl::Number,
    nmax::Number,
    Rarray::Vector,
    zarray::Vector,
    phz::Array;
    lgs::Array=zeros(1,1,1),
    lgs_vac::Array=zeros(1,1,1),
    threads::Bool=true,
    return_coeffs::Bool=false
    )
    if lgs == zeros(1,1,1)
        lgs = mode_pyramid_load(lg, nmax, Nin, Rarray[1], w0=w0)
    end
    if lgs_vac == zeros(1,1,1)
        lgs_vac = mode_pyramid_load(
            lg,
            nmax,
            Nin,
            Rarray[end],
            z=zarray[end],
            wvl=wvl,
            w0=w0
            )
    end
    lgs_turb = modeset_propagation(
        lgs,
        wvl,
        Rarray,
        zarray,
        phz,
        threads=threads
        )
    tmat = transmission_square(lgs_vac, lgs_turb)
    if return_coeffs
        uouts, s, vins, v = svd_lgs_square(lgs, lgs_turb, tmat, return_coeffs=true)
        return uouts, s, vins, tmat, v
    else
        uouts, s, vins= svd_lgs_square(lgs, lgs_turb, tmat)
        return uouts, s, vins, tmat
    end  
end

"""
    transmission_pixel(lgs, phz, wvl, Rarray, zarray)

Transmission matrix of a channel with output pixel basis.
"""
function transmission_pixel(
    lgs::Array,
    phz::Array,
    wvl::Number,
    Rarray::Vector,
    zarray::Vector
    )
    N = size(phz, 2)
    jlmax = size(lgs, 1)
    Nin = size(lgs, 2)
    tmat = zeros(eltype(lgs), Nin^2, jlmax)
    Threads.@threads for i in 1:jlmax
        tmat[:, i] = vec(propagate(lgs[i, :, :], wvl, Rarray, zarray, phz))
        println("get_transmission: propagated mode $i")
    end
    return tmat
end

"""
    svd_lgs_pixel(lgs, tmat)

Singular modes via SVD of a transmission matrix `tmat` with output
pixel basis.
"""
function svd_lgs_pixel(
    lgs::Array,
    tmat::Matrix;
    return_coeffs::Bool=false
    )
    jlmax = size(lgs, 1)
    Nin = size(lgs, 2)
    u, s, v = svd(tmat, full=false)
    println(size(u))
    println(size(s))
    println(size(v))
    vins = zeros(eltype(lgs), jlmax, Nin, Nin)
    for i in 1:jlmax
        for j in 1:jlmax
            vins[i, :, :] += v[j, i] .* lgs[j, :, :]
        end
    end
    println(size(v[1, :]))
    println(size(lgs[1, :, :]))
    uouts = reshape(transpose(u), (jlmax, Nin, Nin))
    if return_coeffs
        return uouts, s, vins, u, v
    else
        return uouts, s, vins
    end
end


# V. Generalized Wigner-Smith operator (experimental beta)

"""
    detla_tmat(tmat1, tmat2, deltat)

Difference quotient of transmission matrices. BETA!
"""
function delta_tmat(
    tmat1::Matrix,
    tmat2::Matrix,
    deltat::Real
    )
    (tmat2 - tmat1)/deltat
end

"""
    gws_q(tmat0, tmat1, tmat2, deltat)

Generalized Wigner-Smith Q. BETA!
"""
function gws_q(
    tmat0::Matrix,
    tmat1::Matrix,
    tmat2::Matrix,
    deltat::Real
    )
    dtmat = delta_tmat(tmat1, tmat2, deltat)
    tdt = inv(tmat0) * dtmat
    tdthc = adjoint(tdt)
    -1im/2*(tdt - tdthc)
end

"""
    principal_modes(gws, lgs)

Principal modes. BETA!
"""
function principal_modes(
    gws::Matrix,
    lgs::Array
    )
    jlmax = size(lgs, 1)
    Nin = size(lgs, 2)
    thetas, vr = eigen(gws)
    pmodes = zeros(eltype(lgs), jlmax, Nin, Nin)
    for i in 1:jlmax
        for j in 1:jlmax
            pmodes[i, :, :] += vr[j, i] .* lgs[j, :, :]
        end
    end
    return thetas, pmodes
end

end