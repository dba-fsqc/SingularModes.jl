using SingularModes

# output directory for plots
datadir = replace(@__FILE__, ".jl" => "")
mkpath(datadir)

# load parameters from file
include("parameters.jl")

# generate hybrid phase screens
phz = phz_x(spec, N / Nin * Rarray, N=N, kind="hd", rs=rs, jzmax=jzmax, r0tot=r0)

# plot first phase screen 
ps = ps_plot(phz[1, :, :], R=N / Nin * Rarray[2])
savefig(ps, "$datadir/phase-screen.pdf")

# generate an individual LG mode with p=0 and l=3
uin = lg_mode(p=0, l=3, w0=w0, R=Rin)

# plot the input mode
pin = mode_plot(uin, R=Rin)
savefig(pin, "$datadir/input-mode.pdf")

# propagate the mode via the split-step method
uout = propagate(uin, wvl, Rarray, zarr, phz)

# plot the output mode
pout = mode_plot(uout, R=Rarray[end])
savefig(pout, "$datadir/output-mode.pdf")
