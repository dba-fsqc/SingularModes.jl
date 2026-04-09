# Statistics.jl
# This module contains all functions for the statistical evaluation of
# phase screens.
# Depends on `BasicFunctions.jl`, `ModeGeneration.jl`, and `Phasescreens.jl`.
# v3.0

module Statistics

export radial_avg, phase_struct_rd, zernike_covariance

using Random

using ..BasicFunctions
using ..ModeGeneration
using ..PhaseScreens


"""
    points_in_circ(radius, samples)

Unifromly distributed random points confined by a circle.
"""
function points_in_circ(
    radius::Number,
    samples::Integer
    )
    t = rand(samples)
    u = rand(samples)
    x = radius.*sqrt.(t).*cos.(2*pi*u)
    y = radius.*sqrt.(t).*sin.(2*pi*u)
    return x, y
end

"""
    radial_avg(a, d)

Average of a 2D array along centered circles.
"""
function radial_avg(
    a::Matrix,
    d::Number
    )
    n = size(a, 1)
    n2 = floor(Integer, n/2)  # assume square array
    ra = collect(0:n2-1) * d
    x = collect(-n2:n2-1)
    x, y = meshgrid(x, x)
    RR = sqrt.(x .^ 2 .+ y .^ 2)
    ma = zeros(eltype(a), n2)
    for r in 1:n2
        ma[r] = mean(a[(RR .>= r - 0.5) .& (RR .< r + 0.5)])
    end
    return ra, ma
end

"""
    chunk_avg(a, c)

Split 1D array into chunks of size c and compute their mean.
"""
function chunk_avg(
    a::Vector,
    c::Integer
    )
    nc = floor(Int, size(a, 1)/c)
    inds = collect(0:nc) .* c .+ 1
    mean.(
        getindex.(Ref(a), (:).(inds[1:end-1], inds[2:end].-1))
    )
end

"""
    phase_struct_ft(a, R)

Radial structure function of a 2D array based on Fourier Transform.
Needs zeropadding.
Not recommended to use, rather use [`phase_struct_rd`](@ref).
"""
function phase_struct_ft(
    ph::Matrix,
    R::Number
    )
    N = size(ph, 1)
    d = 2*R / N
    dph2 = abs.(struct_fcn2_ft(a, ones(N, N), ph))
    return radial_avg(dph2, d)
end

"""
    phase_struct_rd(ph, R; rs=0, samples=10^8)

Radial phase structure function of a 2D array based on random sampling.
"""
function phase_struct_rd(
    ph::Matrix,
    R::Number;
    rs::Integer=0,
    samples::Integer=10^8
    )
    Random.seed!(rs)
    N = size(ph, 1)
    radius = floor(Integer, size(ph, 1) / 2)
    x, y = points_in_circ(radius, samples)
    x = floor.(Integer, x) .+ radius .+ 1
    y = floor.(Integer, y) .+ radius .+ 1
    r = sqrt.(diff(x).^2 + diff(y).^2) * 2*R/N
    phs = getindex.(Ref(ph), x, y)
    rinds = sortperm(r)
    radii = r[rinds]
    dphs = (diff(phs).^2)[rinds]
    dradii = diff(radii)
    inds = vcat(0, findall(>(0), dradii), length(dradii))
    inds[1] = 1
    dph2m = mean.(
        getindex.(Ref(dphs), (:).(inds[1:end - 1], inds[2:end] .- 1))
        )
    uradii = unique(radii)
    c = trunc(Integer, samples / 10^5)
    return chunk_avg(uradii, c) , chunk_avg(dph2m, c)
end

"""
    zernike_covariance(
        pow_spec;
        nmaxScreen=9,
        nmaxZernike=10,
        samples=10,
        rs=0,
        N=1024,
        R=1,
        zks=zeros(1,1,1),
        kwargs...
        )

Zernike variance of screens for a given power spectrum.
See also [`SingularModes.PhaseScreens.zk_overlap`](@ref),
[`SingularModes.PhaseScreens.ph_x`](@ref).
"""
function zernike_covariance(
    pow_spec::Function;
    nmaxScreen::Integer=9,
    nmaxZernike::Integer=10,
    samples::Integer=10,
    rs::Integer=0,
    N::Integer=1024,
    R::Number=1,
    zks::Array=zeros(1,1,1),
    kwargs...
    )
    if samples <= 1
        samples = 2
    end
    Random.seed!(rs)
    rss = rand(1:10^9, samples)
    jzmaxScreen = nmax2j(nmaxScreen)
    jzmaxZernike = nmax2j(nmaxZernike)
    if zks == zeros(1,1,1)
        zksZernike = mode_pyramid_load(
            zk,
            nmaxZernike,
            N,
            R;
            kwargs...
            )
    end
    zksScreen = mode_pyramid_load(zk, nmaxScreen, N, R; kwargs...)
    as = zeros(samples, jzmaxZernike - 1)
    Threads.@threads for i in 1:samples
        screen = ph_x(
            pow_spec,
            N=N,
            R=R,
            rs=rss[i],
            zks=zksScreen,
            jzmax=jzmaxScreen;
            kwargs...
            )
        as[i, :] = zk_overlap(screen, zksZernike)
    end
    return cov(as)
end

end