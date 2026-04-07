# # docs/make.jl
# using Documenter
# using SingularModes

# makedocs(
#     modules=[SingularModes],
#     sitename="SingularModes.jl documentation",
#     authors="David Ba",
#     format=Documenter.HTML(),
# )

using Documenter, SingularModes

makedocs(
    sitename="SingularModes.jl",
    authors="David Bachmann",
    modules=[
        SingularModes.BasicFunctions,
        # SingularModes.FigureDefaults,
        # SingularModes.ModeGeneration,
        # SingularModes.PhaseScreens,
        # SingularModes.Plotting,
        # SingularModes.PowerSpectra,
        # SingularModes.Propagation,
        # SingularModes.Statistics
        ],
    format=Documenter.HTML(
        edit_link="main",
    ),
    )

deploydocs(
    repo="github.com/dba-fsqc/SingularModes.jl.git",
    devbranch="main"
)
