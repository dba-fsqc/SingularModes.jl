using SingularModes
using Test
using FFTW

@testset "SingularModes.jl" begin
    ft([1.0, 2.0, 3.0], 1) ≈ fft([1.0, 2.0, 3.0], 1)
end