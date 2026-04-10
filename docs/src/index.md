# SingularModes.jl

Software package written in [Julia](https://julialang.org) for the propagation of transverse laser modes through turbulent media with particular focus on singular modes.
By construction, singular modes are the optimal communication modes which are obtained via a singular value decomposition of a given channel's transmission matrix.

## Project status
This package emerged from a PhD project and is based on the following publications:

- Bachmann, D., Isoard, M., Shatokhin, V. N., Sorelli, G., Treps, N., and Buchleitner, A. (2023). Highly Transmitting Modes of Light in Dynamic Atmospheric Turbulence. *Physical Review Letters*, **130**(7):073801. [10.1103/physrevlett.130.073801](https://doi.org/10.1103/PhysRevLett.130.073801).
		
- Bachmann, D., Klug, A., Isoard, M., Shatokhin, V. N., Sorelli, G., Buchleitner, A., and Forbes, A. (2024). Universal Crosstalk of Twisted Light in Random Media. *Physical Review Letters,* **132**(6):063801. [10.1103/physrevlett.132.063801](https://doi.org/10.1103/physrevlett.132.063801).
		
- Bachmann, D., Isoard, M., Shatokhin, V. N., Sorelli, G., and Buchleitner, A. (2025). Accurate Zernike-corrected phase screens for arbitrary power spectra. *Optical Engineering,* **64**(5):058102. [10.1117/1.OE.64.5.058102](https://doi.org/10.1117/1.OE.64.5.058102).


## Installation

### Prerequisites: Install Julia

You need a running copy of [Julia](https://julialang.org) to use this package. We recommend to install it via the package manager of your system, e.g., `pacman` or `apt-get` on Arch or Debian-based Linux, respectively.

```bash
# Arch linux
sudo pacman -S julia
```

```bash
# Debian-based Linux
sudo apt-get install julia
```

On macOS, Julia can be installed via the package manager [Homebrew](https://brew.sh) which is recommended for easy updates later on. Install Homebrew and Julia:

```bash
# macOS
# Install Homebrew
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
# Install Julia
brew install julia
```

### Install the SingularModes.jl package

Once Julia is running on your system, the SingularModes package can straightforwardly be installed via Julia's internal package manager.

```bash
julia> ]
pkg> add SingularModes
```

## Package structure

This software package is structured as follows:

```
.
│
├── src                    [Julia modules]
│   ├── BasicFunctions.jl  [library of general functions and definitions]
│   ├── FigureDefaults.jl  [default settings for plotting]
│   ├── ModeGeneration.jl  [library for transverse mode generation]
│   ├── Phasescreens.jl    [library for phase screen generation]
│   ├── Plotting.jl        [library for plotting routines]
│   ├── PowerSpectra.jl    [library of power spectra for different media]
│   ├── Propagation.jl     [library for all functions related to propagation]
│   └── Statistics.jl      [library for the statistical analysis of phase screens]
│
├── docs                   [documentation of the package]
│   └── ...
│
├── examples               [code examples]
│   ├── 01-LG-beams... .jl [LG beam propagation with manual parameters]
│   ├── 02-LG-beams.jl     [LG beam propagation]
│   ├── 03-SVD-square.jl   [Singular modes with square transmission matrix]
│   ├── parameters.jl      [sample parameter file]
│   └── ...
│
├── tests                  [code tests]
│   └── ...
│
├── LICENSE                [GPL-3.0-or-later]
├── README.md              [Short documentation]
│
└── ...
```

The constituent Julia modules are documented below.

---

## BasicFunctions.jl

```@autodocs
Modules = [SingularModes.BasicFunctions]
Order   = [:function, :type]
```

## FigureDefaults.jl

```@autodocs
Modules = [SingularModes.FigureDefaults]
Order   = [:function, :type]
```

## ModeGeneration.jl

```@autodocs
Modules = [SingularModes.ModeGeneration]
Order   = [:function, :type]
```

## PhaseScreens.jl

```@autodocs
Modules = [SingularModes.PhaseScreens]
Order   = [:function, :type]
```

## Plotting.jl

```@autodocs
Modules = [SingularModes.Plotting]
Order   = [:function, :type]
```

## PowerSpectra.jl

```@autodocs
Modules = [SingularModes.PowerSpectra]
Order   = [:function, :type]
```

## Propagation.jl

```@autodocs
Modules = [SingularModes.Propagation]
Order   = [:function, :type]
```

## Statistics.jl

```@autodocs
Modules = [SingularModes.Statistics]
Order   = [:function, :type]
```

## Examples

TBC...