# parameters.jl
# This parameter file collects all inital parameters.
# Depends on `BasicFunctions.jl`, `ModeGeneration.jl` and `Propagation.jl`.
# v2.0

using Random
using Printf
using SpecialFunctions
include("BasicFunctions.jl")
include("ModeGeneration.jl")
include("Propagation.jl")

#. ##### BASH PARAMETERS #####
#. 
#. sbtimea=(40:00:00)
#. sbnodes=1
#. sbmema=200000
#. sbmailtype=FAIL
#. sbmailuser=you@mail.com

#. r0a=(30)
#. rsa=$(seq 1 150)
# rsa=(0)

# #### MANUAL PARAMETERS #####

# --- Random seed ---
rs    = 1                                 # random number seed (integer)
Random.seed!(rs)


# --- Geometry ---
Rphys = 0.05                              # examplary physical input aperture radius / m 
nmax  = 9                                 # highest mode family index
Rin   = 1                                 # radius of the input aperture fixed to one
wvl   = 1550e-9 / Rphys                   # wavelength in units of input aperture radius
zzR   = 1                                 # ratio between propagation distance and Rayleigh length

jlmax = nmax2j(nmax)                      # resulting maximal j_lg index, source ModeGeneration
w0    = Rin / sqrt(nmax + 1)              # beam width (to fit into the input aperture)
zR    = π * w0^2 / wvl                    # Rayleigh length
z     = zzR * zR                          # total propagation distance
Rout  = Rin * sqrt(1 + (z / zR)^2)        # output radius of the mode after diffraction
outin = 5 / 4                             # ratio between output and input aperture size
RoutC = outin * Rin                       # output aperture radius


# --- Grid ---                            # (pixel numbers are best chosen to be powers of two)
N     = 1024                              # one dimensional number of pixels (integer)
Nin   = 256                               # pixel sidelength of apertures (integer)

din   = 2 * Rin / Nin                     # input pixel size
dout  = 2 * RoutC / Nin                   # output pixel size
nhigh = outin^2 / zzR^2 * (nmax + 1)^2    # theoretical sum of singular values
ratio = nhigh / jlmax                     # number of highly transmitting modes vs basis modes


# --- Turbulence ---
spec  = "kol_spec_phi"                    # power spectrum
r0fra = 2                                 # fraction of the input radius that determines r0
lquot = 10                                # fraction of the input radius that determines l0
Lfact = 10                                # factor of the input radius that determins L0
sh    = 7                                 # number of subharmonic corrections (integer)
jzmax = 21                                # maximal number of Zernike coefficients (integer)

alpha = 5/3                               # exponent of the (non-)Kolmogorov spectrum
woce  = -0.8                              # mixinf factor for the oceanic spectrum 


r0    = Rin / r0fra                       # Fried parameter
l0    = Rin / lquot                       # inner scale of turbulence
L0    = Lfact * Rin                       # outer scale of turbulence

spec_para = l0
spec_para2 = L0

if r0 == 0                                # no turbulence
    r0=Inf
end
ns    = min_ns(z, wvl, r0, 0.5)           # number of propagation steps, source Propagation
sr2   = sigma_R2(wvl, z, r0)              # Rytov variance, source Propagation


# --- Wind ---
vmean = 60 * Rin                          # mean wind speed
vstd  = 20 * Rin                          # wind speed standard deviation

pretc = 2 / sqrt(24 / 4 * gamma(6 / 5))   # prefactor for tc (alpha=5/3 only)
tc    = pretc * r0 / vmean                # atmospheric coherence time


deltat= tc                                # time interval for GWS derivative
Tarray= [0, 1, 2, 3, 4, 5, 8]*deltat/10   # time steps
Ts    = size(Tarray, 1)


# #### DERIVED PARAMETERS ####

# --- Power spectra ---
spec_para = alpha
spec_para2 = 0
if contains(string(spec), "kol_spec_phi")
    spec_para = alpha
elseif contains(string(spec), "karman_spec_phi")
    spec_para = L0
elseif contains(string(spec), "tat_spec_phi")
    spec_para = L0
    spec_para2 = l0
elseif contains(string(spec), "oceanic_spec_phi")
    spec_para = -woce
end

# --- Propagation steps ---
halfsteps = collect(0:2ns) ./(2ns)
zarray = halfsteps .* z
darray = din .* sqrt.(1 .+ (outin^2 - 1) .* zarray.^2 ./ z^2) 
Rarray = Rin .* sqrt.(1 .+ (outin^2 - 1) .* zarray.^2 ./ z^2) 
    
# --- Coordinate grids ---
x = -Nin/2:Nin/2 - 1
xx = x' .* ones(Nin)
yy = ones(Nin)' .* x
xin, yin = xx .* din, yy .* din
xoutC, youtC = xx .* dout, yy .* dout

# #### Print parameters ####
print_paras(
    datdir,
    Rphys, outin, zzR,
    N, Nin,
    Rin, RoutC, Rout, din, dout, wvl, z,
    w0, nmax, jlmax, zR, nhigh, ratio,
    spec, r0, l0, L0, sh, jzmax, rs, ns, sr2,
    vmean, vstd, tc, deltat, Ts)


