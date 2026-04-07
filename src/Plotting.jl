# Plotting.jl
# This module contains all plotting routines for the `Plots` package.                                                        .
# v1.1

module Plotting

export s_plot, mode_plot, clean_mode_plot, ps_plot, tmat_plot

using Plots
using ColorSchemes

# """
#     minimal_default(x)

# Minimal defaults for plotting, if [`defaults`](@ref) from the module
# `FigureDefaults.jl` is not used.
# """
# function minimal_default(
#     tick_direction = :out,
#     thickness_scaling = 1.5,
#     background_color = :transparent,
#     foreground_color = :black,
#     dpi = 600
#     )
# end

# I. Plotting routines

"""
    s_plot(svs; m=-1)

Singularvalue plot with optional indicator at index `m`.
"""
function s_plot(
    svs::Vector;
    m::Int=-1
    )
    p = plot(1:size(svs, 1), svs, color=:red3)
    if m !== -1
        vline!([m], color=:gray, line=:dash)
        scatter!(
            [m], [svs[m]],
            color=:gray,
            markersize=1,
            markerstrokewidth=1.5
            )
        annotate!((0.85, 0.9), text(round(svs[m], digits=2)))
    end
    return p
end

"""
    phasemap(mode; cmap=hsv, yflip=false, square=true)

Translate (complex) matrix into HSV-RGBA image.
"""
function phasemap(
    mode::Matrix;
    cmap::ColorSchemes.ColorScheme=ColorSchemes.hsv,
    yflip::Bool=false,
    square::Bool=true,
    nonorm::Bool=false
    )

    Nx = size(mode, 1)
    Ny = size(mode, 2)
    if square
        alphas = abs2.(mode)
    else
        alphas = abs.(mode)
    end
    if ! nonorm
    alphas = alphas .- minimum(alphas)
    alphas = alphas ./maximum(alphas)
    else
        alphas = ones(Nx, Ny)
    end
    phasenorm = (angle.(mode) .+ π) ./ (2π)

    pic = zeros(RGBA{Float32}, Nx, Ny)
    for i in 1:Nx, j in 1:Ny
        if !yflip
            ip = Ny - i + 1
        else
            ip = i
        end
        pic[ip, j] = RGBA{Float32}(cmap[phasenorm[i, j]], alphas[i, j])
    end
    return pic
end

"""
    colorbar2d_plot(ratio=5)

2D colorbar for HSV-RGBA images. See also [`phasemap`](@ref).
"""
function colorbar2d_plot(
    ratio::Number=5;
    nogradient=false,
    cmap::ColorSchemes.ColorScheme=ColorSchemes.hsv
    )
    N = 100
    gradient = (0:N)' .* (exp.(-2π*1im.*(-N/2:N/2)/N))
    cpic = phasemap(gradient, cmap=cmap)
    xx = (0:N)./N
    yy = xx * 2π .- π
    if nogradient
        gradient = (ones(N+1))' .* (exp.(-2π*1im.*(-N/2:N/2)/N))
        cpic = phasemap(gradient, nonorm=true, cmap=cmap)
    end
    xticks=[0,1]
    if nogradient
        xticks=false
    end
    p2 = plot(
        xx,
        yy,
        cpic,
        ratio=ratio,
        ymirror=true,
        xlims=(0,1),
        ylims=(-π, π),
        grid=false,
        xticks=xticks,
        yticks = (
            (-π:π/2:π),
            ["+π", "+π/2", "0", "-π/2", "-π"]
            )
        )
    return p2
end

"""
    mode_plot(
        mode; R=1, kind=0, legend=true,
        boundary=true, circle=true, yflip=false, clean=false
        )

Plot intensity (`kind = 0`), phase (`kind = 1`) or combination (`kind = 2`)
of a mode.
"""
function mode_plot(
    mode::Matrix;
    R::Number=1,
    kind::Integer=2,
    legend::Bool=true,
    boundary::Bool=true,
    circle::Bool=true,
    yflip::Bool=false,
    clean::Bool=false,
    gridb::Bool=true,
    noclim=false,
    noabs=false,
    colormap::Symbol=:hsv,
    cmap::ColorSchemes.ColorScheme=ColorSchemes.hsv,
    boundarystyle=:solid,
    showaxis=true,
    kwargs...
    )
    N = size(mode, 1)
    d = 2R/N
    x = (-N/2:N/2)*d
    bf = 1.00
    if circle
        if kind !== 2
            mode = inner_circ(mode, nans=true)
        else
            mode = inner_circ(mode)
        end
    end
    if clean
        showaxis = false
        gridb = false
        legend = false
    end
    if kind == 0  # intensity
        if noabs
            pfunc = mode
        else
            pfunc = abs2.(mode)
        end
        p = heatmap(
            x, x, pfunc,
            yflip=yflip,
            ratio=:equal,
            xlims=[-bf*R, bf*R],
            ylims=[-bf*R, bf*R],
            color=:viridis,
            legend=legend,
            showaxis=showaxis,
            grid=gridb
            )
    elseif kind == 1  # phase
        if noclim
            p = heatmap(x, x, angle.(mode),
                yflip=yflip,
                ratio=:equal,
                xlims=[-bf*R, bf*R],
                ylims=[-bf*R, bf*R],
                color=colormap,
                legend=false,
                showaxis=showaxis,
                grid=gridb
                )
        else
            p = heatmap(x, x, angle.(mode),
                yflip=yflip,
                ratio=:equal,
                xlims=[-bf*R, bf*R],
                ylims=[-bf*R, bf*R],
                color=colormap,
                clims=(-π, π),
                legend=false,
                showaxis=showaxis,
                grid=gridb
                )
        end
        if legend
            p2 = colorbar2d_plot(3, nogradient=true, cmap=cmap)
            p = plot(
                p, p2,
                layout = grid(1, 2, widths=[0.95, 0.05])
                )
        end
    else  # composite
        pic = phasemap(mode)
        p = plot(
            x,
            x,
            pic,
            xlims=[-bf*R, bf*R],
            ylims=[-bf*R, bf*R],
            ratio=:equal,
            showaxis=showaxis,
            grid=gridb,
            yflip=false
            )
        if legend
            p2 = colorbar2d_plot(3)
            p = plot(
                p, p2,
                layout = grid(1, 2, widths=[0.95, 0.05])
                )
        end
    end
    if boundary
        pts = Plots.partialcircle(0, 2π, 100, 1.0 * R)
        x, y = Plots.unzip(pts)
        plot!(
            p,
            Shape(
                x,
                y
            ),
            fillalpha=0,
            linecolor=:black,
            linealpha=1,
            linewidth=1.2,
            legend=false,
            style=boundarystyle
            )
    end
    return p
end

"""
    clean_mode_plot(mode; kind=2)

Shorthand for clean mode plot. See also [`mode_plot`](@ref).
"""
function clean_mode_plot(
    mode::Matrix;
    kind::Int=2
    )
    mode_plot(
        mode,
        R=1,
        kind=kind,
        clean=true,
        boundary=false
    )
end

"""
    ps_plot(ph; R=1, r0=0, legend=true, boundary=true)

Phase screen plot with optional scale bar for `r0 != 0`.
"""
function ps_plot(
    ph::Matrix;
    R::Number=1,
    r0::Number=0,
    legend::Bool=true,
    boundary::Bool=true,
    colormap::Symbol=:hsv,
    grey::Bool=false,
    axes::Bool=false,
    cleanmode::Bool=false
    )
    if grey
        colormap=:Greys
    end
    if legend
        clean = false
    end
    if axes
        cleanmode = false
    end
    phcomplex = exp.(1im*ph)
    p = mode_plot(
        phcomplex,
        R=R,
        kind=1,
        legend=legend,
        boundary=boundary,
        colormap=:Greys,
        cmap=ColorSchemes.Greys,
        clean=cleanmode
        )
    if r0 !== 0
        ypos = -0.93
        xpos = 0.93
        plot!(
            p,
            [xpos - r0, xpos].*R,
            [ypos, ypos].*R,
            color=:black,
            marker=:vline,
            linecolor=:black,
            linewidth=1.5,
            markersize=3,
            legend=false
            )
        annotate!(
            p,
            (xpos - r0/2)*R,
            (ypos+0.1)*R,
            text("\$r_0\$", 10),
            legend=false
            )
    end
    return p
end

"""
    tmat_plot(tmat; grid=true, kind=0, legend=true, yflip=true, square=false)

Transmission matrix plot of different kinds. See also [`mode_plot`](@ref).
"""
function tmat_plot(
    tmat::Matrix;
    grid::Bool=true,
    kind::Integer=0,
    legend::Bool=true,
    yflip::Bool=true,
    square::Bool=false,
    noabs::Bool=false,
    kwargs...
    )
    (N, M) = size(tmat)
    square = false
    if N == M
        square = true
    end
    x = collect(0:N-1) .+ 0.5
    y = collect(0:M-1) .+ 0.5
    xt = x[end] <= 100 ?
        (vcat(0.5, collect(10:10:x[end]) .- 0.5),
        vcat(1, 10:10:floor(Int, x[end]))) :
        :auto
    yt = y[end] <= 100 ?
        (vcat(0.5, collect(10:10:y[end]) .- 0.5),
        vcat(1, 10:10:floor(Int, y[end]))) :
        :auto
    if noabs
        ptmat = tmat
    else
        ptmat = abs.(tmat)
    end
    if kind == 0  # here only absolute value *not abs2*
        p = heatmap(
            y,
            x,
            ptmat,
            yflip=yflip,
            ratio=square ? :equal : :auto,
            xlims=(0, M),
            ylims=(0, N),
            xticks=xt,
            yticks=yt,
            color=:viridis,
            legend=legend,
            grid=false;
            kwargs...
            )
    elseif kind == 1  # phase
        p = heatmap(
            y,
            x,
            xticks=xt,
            yticks=yt,
            angle.(tmat),
            yflip=yflip,
            ratio=square ? :equal : :auto,
            color=:hsv,
            clims=(-π, π),
            legend=legend,
            grid=false;
            kwargs...
            )
    else  # composite
        x = collect(0:N)
        y = collect(0:M)
        pic = phasemap(tmat, yflip=yflip, square=square)
        p = plot(
            y, x, pic,
            ratio=square ? :equal : :auto,
            xticks=xt,
            yticks=yt,
            xlims=(0, M + 0),
            ylims=(0 , N + 0),
            grid=false,
            kwargs...
            )
        if legend
            p2 = colorbar2d_plot(3)
            lay = @layout [a{0.95w} b{0.05w}]
            p = plot(p, p2, layout = lay)
        end
    end
    if grid
        gcol = :red
        offset = 1 # 0.5
        if kind == 2
            gcol = :black
            offset = 1
        end
        i = 2
        di = 2
        if N == M
            while i <= M
                vline!(p, [i-offset], color=gcol, linewidth=0.5, legend=false)
                hline!(p, [i-offset], color=gcol, linewidth=0.5, legend=false)
                i += di
                di += 1
                end
        else
            while i <= M
                vline!(p, [i-offset], color=gcol, linewidth=0.5, legend=false)
                i += di
                di += 1
            end
        end 
    end
    return p
end

end