# SingularModes.jl

Software package written in [Julia](https://julialang.org) for the propagation of transverse laser modes through turbulent media.


## Installation

### Prerequisites: Install Julia

You need a running copy of [Julia](https://julialang.org) to use this package. We recommend to install it via the package manager of your system, e.g., `pacman` or `apt-get` on Arch or Debian-based Linux, respectively.

```bash
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
│   ├── static             [cached arrays of LG/ZK/KL modes for faster computation]
│   │
│   ├── parameters.jl      [sample parameter file]
│   ├── LG-Beams.jl        [sample code for the propagtion of Lagurre-Gaussian beams]
│   ├── Singularmodes.jl   [sample code for the computation of singular modes]
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