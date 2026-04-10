using SingularModes

# output directory for plots
datadir = replace(@__FILE__, ".jl" => "")
mkpath(datadir)

# ===== Initial parameters =====
# --- Geometry ---
Rin   = 5e-2                              # input aperture radius / m
Rout  = 6.25e-2                           # output aperture radius / m
z     = 507                               # propagation distance / m

wvl   = 1550e-9                           # optical wavelength / m
w0    = 1.6e-2                            # beam waist / m
zR    = π * w0^2 / wvl                    # Rayleigh length

# --- Grid ---
N     = 1024                              # phase screen diameter / pixels
Nin   = 256                               # mode diameter / pixels
din   = 2*Rin / Nin                       # input pixel size
dout  = 2*Rout / Nin                      # output pixel size
outin = Rout / Rin                        # aperture size ratio

# --- Turbulence ---
spec  = atm_spec                          # atmospheric Kolmogorov spectrum
r0    = 25e-3                             # Fried parameter / m
jzmax = 21                                # max. incl. Zernike coeff. (int)
rs    = 1                                 # random seed (int)

# --- Propagation steps ---
ns    = min_ns(z, wvl, r0, 0.5)           # number of screens
sr2   = sigma_R2(wvl, z, r0)              # Rytov variance
hss   = collect(0:2*ns) ./ (2*ns)         # halfsteps
zarr  = hss .* z                          # zarray
sarr  = sqrt.(1 .+(outin^2-1)*zarr.^2/z^2)# scale array
darr  = din * sarr                        # pixel size array
Rarray= Rin * sarr                        # radii array

# ===== Propagation of a Lagurre-Gaussian (LG) beam =====
# generate hybrid phase screens
phz = phz_x(spec, N/Nin*Rarray, N=N, kind="hd", rs=rs, jzmax=jzmax, r0tot=r0)

# plot first phase screen 
ps = ps_plot(phz[1, :, :], R=N/Nin*Rarray[2])
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