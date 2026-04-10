using SingularModes

# output directory for plots
datadir = replace(@__FILE__, ".jl" => "")
mkpath(datadir)

# load parameters from file
include("parameters.jl")

# generate a mode pyramid of input LG modes
lgs_in = mode_pyramid(lg, nmax, Nin, Rarray[1], w0=w0)

# generate a mode pyramid of output LG modes
lgs_vac = mode_pyramid(lg, nmax, Nin, Rout, z=zarr[end], w0=w0, wvl=wvl)

# generate hybrid phase screens
phz = phz_x(atm_spec, N/Nin*Rarray, rs=rs, r0tot=r0)

# singular modes and transmission matrix of the channel
uouts, s, vins, tmat = svd_square_wrapper(
    Nin,
    w0,
    wvl,
    nmax,
    Rarray,
    zarr,
    phz,
    lgs=lgs_in,
    lgs_vac=lgs_vac
)

# plot the highest transmitting in- and output modes
pin = mode_plot(vins[1, :, :], R=Rin)
savefig(pin, "$datadir/vins1.pdf")

pout = mode_plot(uouts[1, :, :], R=Rout)
savefig(pout, "$datadir/uouts1.pdf")

# plot the singular value distribution
psv = s_plot(s)
savefig(psv, "$datadir/singular-values.pdf")

# plot the transmission matrix
pt = tmat_plot(tmat)
savefig(pt, "$datadir/tmat.svg") # pdf forces too low resolution