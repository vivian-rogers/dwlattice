# General imports.
using Plots
using Waveforms

"""
    Takes time domain function and returns aₙ,
    where aₙ is an amplitude of X(ω) in frequency domain.
    @params
    x::function
    t::Float64 (1 period)
    N::Int (no. of points)
    n::Int (harmonic index)
"""
function aₙ(x, t::Float64, N::Int, n::Float64)
    # fundamental freq.
    f = 1/t

    # range of points to iterate through.
    pts = LinRange(0.0, t, N)

    # final array
    final = zeros(ComplexF64, N)
    @. final = ℯ^(-im*n*2pi*f*pts) * x(pts)

    return sum(final) / t
end


# Magnitude for every nᵗʰ harmonics.
function Mag_Harmonics(x, t::Float64, N::Int, total::Float64)
    # Array of harmonic #s.
    range = LinRange(1.0, total, convert(Int, total))

    # A of every harmonic
    A = zeros(ComplexF64, convert(Int, total))
    @. A = aₙ(x, t, N, range)

    # Get magnitude of each harmonics.
    mag = zeros(convert(Int, total))
    @. mag = abs(A)

    return range, mag
end

squ_harm, squ_mag = Mag_Harmonics(squarewave, 2pi, 50, 11.0)

# Plot
plot(squ_harm, squ_mag, line=:stem)
title!("Magnitude of Harmonics")
xlabel!("n (ᵗʰ harmonic)")
ylabel!("Magnitude")
savefig("img/Harmonic_Mag.png")