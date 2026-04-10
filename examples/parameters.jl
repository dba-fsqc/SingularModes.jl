# parameters.jl
# This parameter file collects all inital parameters.
# Depends on `BasicFunctions.jl`, `ModeGeneration.jl` and `Propagation.jl`.
# v3.0

using SingularModes
using SpecialFunctions

set_fig_defaults()                        # set figure defaults for plotting 

# ===== MANUAL PARAMETERS =====

# --- Geometry ---
Rin   = 5e-2                              # input aperture radius / m
nmax  = 9                                 # highest mode family index
wvl   = 1550e-9                           # optical wavelength / m
zzR   = 1                                 # propagation distance / Rayleigh length

jlmax = nmax2j(nmax)                      # maximal LG superindex
w0    = Rin / sqrt(nmax + 1)              # beam width to fit into input aperture
zR    = π * w0^2 / wvl                    # Rayleigh length
z     = zzR * zR                          # total propagation distance
outin = 5 / 4                             # ouput aperture / input aperture size
Rout  = outin * Rin                       # output aperture radius
RoutD = Rin * sqrt(1 + (z / zR)^2)        # theo. output radius diffraction only

# --- Grid ---
N     = 1024                              # phase screen diameter / pixels
Nin   = 256                               # mode diameter / pixels

din   = 2*Rin / Nin                       # input pixel size
dout  = 2*Rout / Nin                      # output pixel size
nhigh = outin^2 / zzR^2 * (nmax + 1)^2    # theo. sum of singular values (SV)
ratio = nhigh / jlmax                     # sum SV / number of basis modes

# --- Turbulence ---
spect = "kol_spec_phi"                    # atmospheric Kolmogorov spectrum
wturb = 2                                 # distortion strength (Rin/r0)
lquot = Inf                               # Rin / inner scale
Lfact = Inf                               # outer scale / Rin
alpha = 5/3                               # (non)-Kolmogorov exponent
woce  = -4/5                              # oceanic mixing parameter

r0    = Rin / wturb                       # Fried parameter
if r0 == 0                                # no turbulence: r0 -> Inf
    r0 = Inf
end
l0    = Rin / lquot                       # inner scale of turbulence
L0    = Lfact * Rin                       # outer scale of turbulence

rs    = 1                                 # random seed (int)
jzmax = 21                                # max. incl. Zernike coeff. (int)
sh    = 7                                 # number of subharmonic corrections

# --- Wind ---
vmean = 60 * Rin                          # mean wind speed
vstd  = 20 * Rin                          # wind speed standard deviation

pretc = 2 / sqrt(24 / 4 * gamma(6 / 5))   # prefactor for tc (alpha=5/3 only)
tc    = pretc * r0 / vmean                # atmospheric coherence time (tc)
deltat= tc                                # time interval
Tarray= [0, 1, 2, 3, 4, 5, 8]*deltat/10   # time steps
Ts    = size(Tarray, 1)                   # number of time steps


# ===== DERIVED PARAMETERS =====

# --- Propagation steps ---
ns    = min_ns(z, wvl, r0, 0.5)           # number of screens
sr2   = sigma_R2(wvl, z, r0)              # Rytov variance
hss   = collect(0:2*ns) ./ (2*ns)         # halfsteps
zarr  = hss .* z                          # zarray
sarr  = sqrt.(1 .+(outin^2-1)*zarr.^2/z^2)# scale array
darr  = din * sarr                        # pixel size array
Rarray= Rin * sarr                        # radii array

# --- Power spectra ---
spec = atm_spec
spec_para = alpha
spec_para2 = 0
if contains(string(spect), "kol_spec_phi")
    spec_para = alpha
elseif contains(string(spect), "karman_spec_phi")
    spec_para = L0
elseif contains(string(spect), "tat_spec_phi")
    spec_para = L0
    spec_para2 = l0
elseif contains(string(spect), "oceanic_spec_phi")
    spec = oce_spec
    spec_para = -woce
end

# ===== Print parameters =====
print_paras(
    outin, zzR, N, Nin,
    Rin, Rout, RoutD, din, dout, wvl, z,
    w0, nmax, jlmax, zR, nhigh, ratio,
    spect, r0, l0, L0, sh, jzmax, rs, ns, sr2,
    vmean, vstd, tc, deltat, Ts
    )