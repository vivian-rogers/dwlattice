# General imports.
using Plots
#using Waveforms



function genFiring_neuron_DWpos(f::Float64)
    pos(t) = f*mod(t,1/f) -1/2
    return pos
end

"""
    Takes time domain function and returns aₙ,
    where aₙ is an amplitude of X(ω) in frequency domain.
    @params
    x::function
    t::Float64 (1 period)
    N::Int (no. of points)
    n::Int (harmonic index)
"""
function aₙ(x::Function, t::Float64, Npts::Int, n::Int)
    # fundamental freq.
    f = 1/t

    # range of points to iterate through.
    pts = LinRange(0.0, t, Npts)

    # final array
    final = zeros(ComplexF64, Npts)
    @. final = ℯ^(-im*n*2pi*f*pts) * x(pts)

    return sum(final) / t
end


# Magnitude for every nᵗʰ harmonics.
function Mag_Harmonics(x::Function, t::Float64, N::Int, nharmonics::Int)
    # Array of harmonic #s.
    range = [i for i = -nharmonics:1:nharmonics]
    #range = [-2,-1,0,1,2]
    #range = LinRange(1.0, total, convert(Int, total))

    # A of every harmonic
    A = zeros(ComplexF64, nharmonics*2+1)
    @. A = aₙ(x, t, N, range)

    # Get magnitude of each harmonics.
    #mag = zeros(nharmonics*2+1)
    mag = A
    #mag = real.(A)
    #@. mag = abs(A)

    return range, mag
end


function invFT(harmonic_indices::Vector{Int},aₙcoeffs::Vector{ComplexF64},f₀::Float64)
    function y(t::Float64)
        return sum([aₙcoeffs[i]*ℯ^(im*f₀*2*π*i*t) for i in harmonic_indices])
    end
    return y
end

f₀ = 1.0
squ_harm, squ_mag = Mag_Harmonics(genFiring_neuron_DWpos(f₀), 1/f₀, 50, 50)
#squ_harm, squ_mag = Mag_Harmonics(squarewave, 2pi, 50, 11.0)

function Plot_Spectrum(harmonics, magnitude)
    # Plot
    plot(harmonics, real.(magnitude), line=:stem,c="blue")
    plot!(harmonics, imag.(magnitude), line=:stem,c="red")
    title!("Magnitude of Harmonics")
    xlabel!("n (ᵗʰ harmonic)")
    ylabel!("Magnitude")
end

Plot_Spectrum(squ_harm, squ_mag)
savefig("vs_lattice/img/Harmonic_Mag.png")