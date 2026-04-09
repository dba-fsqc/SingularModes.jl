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

using .ModeGeneration
#! format: off
export
    pl2j,
    nm2j,
    mn2j,
    j2pl,
    j2nm,
    j2mn,
    j2nmax,
    nmax2j,
    noll_pyramid,
    r_z,
    zk,
    laguerre,
    lg,
    hg_1d,
    hg,
    besselbeam,
    besselgaussian,
    make_grid,
    _default_N,
    points2grid,
    lg_mode,
    hg_mode,
    bg_mode,
    zk_mode,
    mode_pyramid,
    mode_pyramid_load,
    zka_load,
    pcl,
    pcl_mod,
    pcl_arc,
    pcl_geo,
    pcl_lim,
    bcl,
    rho_mean_bg,
    pcl_bg
#! format: on

using .PhaseScreens
#! format: off
export 
    pFq_reg,
	zk_cov_atm,
	zk_cov_num,
	zk_covmat,
	ph_ft,
	ph_sh,
	ph_zk,
	zk_overlap,
	ph_hd,
	ph_kl,
	ph_x,
	phz_x,
	screen_factor,
	fourier_shift,
	wind_shift,
	get_phzT
#! format: on

using .Plotting
#! format: off
export 
	s_plot,
	phasemap,
	colorbar2d_plot,
	mode_plot,
	clean_mode_plot,
	ps_plot,
	tmat_plot
#! format: on

using .PowerSpectra
#! format: off
export 
    atm_spec,
	_oce_delta,
	oce_spec,
	ps2struct,
	atm_struct
#! format: on

using .Propagation
#! format: off
export 
    fraunhofer_prop,
	lens_against_ft,
	ang_spec_prop,
	rayleigh_length,
	sigma_R2,
	sigma_R2_1,
	cn2,
	cn2sr,
	r0step,
	min_ns,
	multi_step_variable,
	propagate,
	modeset_propagation,
	transmission_square,
	svd_lgs_square,
	svd_square_wrapper,
	transmission_pixel,
	svd_lgs_pixel,
	delta_tmat,
	gws_q,
	principal_modes
#! format: on

using .Statistics
#! format: off
export
	points_in_circ,
	radial_avg,
	chunk_avg,
	phase_struct_ft,
	phase_struct_rd,
	zernike_covariance
#! format: on

end
