# PowerSpectra.jl
# This module contains the definitions of various power spectra and
# structure functions.
# v4.1

module PowerSpectra

export atm_spec, oce_spec, ps2struct, atm_struct

using SpecialFunctions
using QuadGK


# I. Power spectra

# all spectra and structure functions have r_0^{-5/3} factored out

# prefactor = 0.49
const prefactor = (2^(2/3)*gamma(11/6)^2) / (π^2) * (24/5 * gamma(6/5))^(5/6)

"""
    atm_spec(k; L0=Inf, l0=0, alpha=5/3)

Atmospheric phase power spectrum. E.g. for
- `L0=Inf, l0=0, alpha=5/3`: Kolmogorov
- `L0=Inf, l0=0, alpha=2`  : quadratic approximation
- `L0<Inf, l0=0, alpha=5/3`: von Karman
- `L0<Inf, l0>0, alpha=5/3`: Tatarskii.
"""
function atm_spec(
    k::Number;
    L0::Number=Inf,
    l0::Number=0,
    alpha::Number=5/3,
    kwargs...   # to dump potentially other spectral parameters
    )
    k0 = 2 * π / L0
    fac = (sqrt(3) * gamma(8/3) / 8 / π)^(-3/4)
    km = fac / l0
    term1 = (k^2 + k0^2)^((-alpha -2)/2) * exp(-k^2 / km^2)
    return prefactor * term1

end

"""
    _oce_delta(k, eta, C1)

Help function for oce_spec. See also [`oce_spec`](@ref).
"""
function _oce_delta(
    k::Number,
    eta::Number,
    C1::Number
    )
    return 1.5 * C1^2 * (k * eta)^(4/3) + C1^3 * (k * eta)^2
end

"""
    oce_spec(k; w=-0.5, eta=0.1)

Oceanic phase power spectrum based on Korotkova (2018).
"""
function oce_spec(
    k::Number;
    w::Number=-0.5,
    eta::Number=0.1,
    kwargs...
    )
    w = -abs(w) # ensure negativity of w
    C0 = 0.72  # Obukhov-Corrsin constant
    C1 = 2.35  # fitting parameter
    PrT = 7  # Prandtl number temperature
    PrS = 700  # Prandtl number salt
    AT = C0 * C1^(-2) / PrT
    AS = C0 * C1^(-2) / PrS
    ATS = 0.5 * C0 * C1^(-2) * (1 / PrT + 1 / PrS)
    aT = 2.6e-4
    aS = 1.7e-4
    #eta = 0.1  # innerscale / m
    epsilon = 1e-4  # 1e-2 (surface) to 1e-8 (deep) / m^2 s^-3
    chiT = 1e-6  # 1e-4 (surface) to 1e-10 (deep midwater) / K^2 s^-1
    
    chiS = aT^2 / w^2 / aS^2 * chiT
    chiTS = (chiT * chiS)^(1/2)
    chiN = aT^2 * chiT * (w - 1)^2 / w^2
    
    first = C0 * (4 * π)^(-1) * aT^2 * chiN * epsilon^(-1/3) * k^(-11/3)
    second = (1 + C1 * (k * eta)^(2/3))
    third1 = w^2 * exp(-AT * _oce_delta(k, eta, C1))
    third2 = exp(-AS * _oce_delta(k, eta, C1))
    third3 = 2 * w * exp(-ATS * _oce_delta(k, eta, C1))
    third = third1 + third2 - third3
    units = (1e6)^3  # numerical factor to pull it to larger numbers
    return first * second * third * units
end


# II. Structure functions

"""
    ps2struct(pow_spec, r; kwargs...)

Numerical conversion of a given phase power spectrum to its structure
function, the `kwargs...` are passed to `pow_spec`.
See also [`atm_spec`](@ref), [`oce_spec`](@ref).
"""
function ps2struct(
    pow_spec::Function,
    r::Vector;
    kwargs...
    )
    # numerical conversion from power spectrum to structure function
    function f(k::Number, r::Number)
        k * pow_spec(k; kwargs...) * (1 - besselj(0, k * r))
    end
    n = size(r, 1)
    out = zeros(n) 
    Threads.@threads for i in 1:n
        out[i] = 4 * π * quadgk(k -> f(k, r[i]), 0, Inf, rtol=1e-6)[1]
    end
    return out
end

"""
    atm_struct(r; L0=Inf)

Theoretical structure function of Kolmogorov (`L0=Inf`) or von Karman
(`L0 < Inf`) turbulence.
"""
function atm_struct(
    r::Number;
    L0::Number=Inf
    )
    if L0 == Inf
        return 2 * (24/5 * gamma(6/5))^(5/6) * r^(5/3)
    else
        k0 = 2 * π / L0
        first1 = 2 * gamma(11/6) / 2^(5/6) / π^(8/3)
        first2 = (24/5 * gamma(6/5))^(5/6)
        second = (2 * π / 1 / k0)^(5/3)
        third = gamma(5/6) / 2^(1/6) - (r * k0)^(5/6) * besselk(5/6, r * k0)
        return first1 * first2 * second * third
    end
end

end