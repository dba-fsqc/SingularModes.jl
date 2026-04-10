# FigureDefaults.jl
# This module specifies the defaults for figure plotting with the `Plots` package.
# v1.0

module FigureDefaults

using Plots
using LaTeXStrings
plot_font = "Computer Modern"

"""
    set_fig_defaults()

Set defaults for plotting figures with the `Plots` package.
"""
function set_fig_defaults()
    default()
    default(
        fontfamily = plot_font,
        linewidth = 2, 
        framestyle = :box, 
        background_color = :transparent,
        label = nothing, 
        grid = true,
        gridalpha = 1,
        background_color_legend = :white,
        foreground_color_grid = RGB(0.8,0.8,0.8),
        thickness_scaling = 1.5,
        tick_direction = :out
    )
end

end