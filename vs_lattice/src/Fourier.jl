# General imports.
using Plots
using Waveforms

#function genFiring_neuron_DWpos(f::Float64)
#    pos(t) = trianglewave(t*2*π*f)
#    #pos(t) = cos(2*π*f*t) + 1/2*sin(2*π*3*f*t) 
#    #pos(t) = f*mod(t,1/f) -1/2
#    return pos
#end

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
    pts = LinRange(0.0, t - t/Npts, Npts)

    # final array
    #final = ℯ.^(-im*n*2*π*f*pts) .* x.(pts)
    final = zeros(ComplexF64, Npts)
    @. final = ℯ^(-im*n*2pi*f*pts) * x(pts) / Npts

    return (sum(final) / t)
end


# Magnitude for every nᵗʰ harmonics.
function FourierSeries(x::Function, t::Float64, N::Int, nharmonics::Int)
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
        yval::ComplexF64 = 0
        for i in 1:size(harmonic_indices)[1]
            fᵢ = f₀*harmonic_indices[i] 
            yval += aₙcoeffs[i]*ℯ^(im*2pi*fᵢ*t)
        end
        return real(yval) 
            #sum([
            #    aₙcoeffs[i]*ℯ^(im*f₀*2*π*harmonic_indices[i]*t) for i in size(harmonic_indices)[1]
            #    ]))
    end
    return y
end


function Plot_Spectrum(harmonics, magnitude, title="")
    # Plot
    plot(harmonics, real.(magnitude), line=:stem,c="blue")
    plot!(harmonics, imag.(magnitude), line=:stem,c="red")
    title!("Magnitude of Harmonics")
    xlabel!("n (ᵗʰ harmonic)")
    ylabel!("Magnitude")
    savefig("vs_lattice/img/"*title*" FourierSeries.png")
end

#f₀ = 1.0
#x = genFiring_neuron_DWpos(f₀)

#nharmonics = 10; npts = 50
#harmonics, coeffs = FourierSeries(x, 1/f₀, npts, nharmonics)

#approx_x = invFT(harmonics,coeffs,f₀)

#Plot_IO(x,approx_x,2*1/f₀)
#squ_harm, squ_mag = Mag_Harmonics(squarewave, 2pi, 50, 11.0)
#Plot_Spectrum(harmonics, coeffs)
# plotIO(harmonics, coeffs)