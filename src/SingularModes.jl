module SingularModes

include("BasicFunctions.jl")
include("FigureDefaults.jl")
include("ModeGeneration.jl")
include("PhaseScreens.jl")
include("Plotting.jl")
include("PowerSpectra.jl")
include("Propagation.jl")
include("Statistics.jl")

using .BasicFunctions
export ft, ift, zeropad, inner_circ, rect, cart2pol, pol2cart
export meshgrid
export print_paras, msum, norm_mode, wopt, rayleigh_waist
export norm_coeff

using .ModeGeneration
export pl2j, nm2j, mn2j, j2pl, j2nm, j2mn, j2nmax, nmax2j, noll_pyramid
export zk, lg, hg, besselbeam, besselgaussian, points2grid
export lg_mode, hg_mode, bg_mode, zk_mode, mode_pyramid, zka_load
export mode_pyramid_load
export pcl, pcl_mod, pcl_arc, pcl_geo, pcl_lim, bcl
export rho_mean_bg, pcl_bg

using .PhaseScreens
export zk_covmat, ph_ft, ph_sh, ph_zk, ph_hd, ph_kl, ph_x
export phz_x, fourier_shift, wind_shift, get_phzT

using .Plotting
export s_plot, mode_plot, clean_mode_plot, ps_plot, tmat_plot

using .PowerSpectra
export atm_spec, oce_spec, ps2struct, atm_struct

using .Propagation
export fraunhofer_prop, lens_against_ft, ang_spec_prop
export rayleigh_length, sigma_R2, min_ns, multi_step_variable
export propagate, modeset_propagation, transmission_square
export svd_lgs_square, svd_square_wrapper, transmission_pixel
export svd_lgs_pixel

using .Statistics
export radial_avg, phase_struct_rd, zernike_covariance

end
