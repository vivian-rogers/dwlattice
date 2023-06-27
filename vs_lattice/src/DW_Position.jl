# Specific Imports.
include("DW_Lattice.jl")
include("Fourier.jl")

# units
μm = 10^-6; GHz = 1; 

testSystem = DWLattice(2, [0, 2.5]*μm, [1,-1], 5*μm, [0.1,0.1]*GHz, [10,15]*GHz, true, 10^8, 2.5*μm);
H_AFM_racetrack = constructHamiltonian(testSystem,2)

AFM2R_lattice = init(n_racetracks=2, racetrack_positions=[0, 2.5]*μm, orientations=[1,-1], 
    a=5*μm, γ=0.1*[1,1]*GHz, ω₀=0.0*[9,8]*GHz, PBC=true, C=10.0*GHz, R₀=2.5*μm);
AFM3R_lattice = init(n_racetracks=3, racetrack_positions=[0, 2.8, 5.0]*μm, orientations=[1,-1,1], 
    a=7.5*μm, γ=0.5*[1,0.5,1]*GHz, ω₀=10*[0,2.0,2.1]*GHz, PBC=true, C=400.0*GHz, R₀=2.5*μm);
FM1R_lattice = init(n_racetracks=1, racetrack_positions=[0]*μm, orientations=[1], 
    a=5*μm, γ=0.5*[1]*GHz, ω₀=10*[2.0]*GHz, PBC=true, C=80.0*GHz, R₀=2.5*μm);

"""
Returns the displacement of DW depending on input signal.

@params
x (input signal function)
t (period)
N (# of points collected during 1 period)
total (# harmonics)
T (Transfer Function)
"""
function genFiring_neuron_DWpos(f::Float64)
    pos(t) = trianglewave(t*2*π*f)
    #pos(t) = cos(2*π*f*t) + 1/2*sin(2*π*3*f*t) 
    #pos(t) = f*mod(t,1/f) -1/2
    return pos
end
#start_f = 0.0000001*GHz
#stop_f = 15.0*GHz
#function genFiring_neuron_DWpos(f::Float64)
#    pos(t) = f*mod(t,1/f) -1/2
#    return pos
#end

ΔT = 2*π/(6.66*GHz); f₀ = 1/ΔT
x = genFiring_neuron_DWpos(1/ΔT)
ϵ = 0.001

function DWPosition(system::DWLattice, n_lattice::Int, race_i::Int, race_j::Int, 
    x::Function, f₀::Float64, nharmonics::Int=10, NNs::Int=100, npts::Int=150)

    Δt = 1/f₀
    T = transferFunc(system, n_lattice, race_i, race_j, 0.0,0.0, 1, NNs, false)

    # Array of harmonic #s and corresponding Magnitude 
    harmonic_indices, aₙvals = FourierSeries(x, Δt, npts, nharmonics)
    #Plot_Spectrum(harmonic_indices,aₙvals,"anvals")
    ωvals = (2π/Δt).*harmonic_indices .+ ϵ*GHz
    bₙvals = aₙvals.*T.(ωvals)
    #BodePlot(T,-nharmonics*2*π/Δt,nharmonics*2*π/Δt,200)
    #Plot_Spectrum(harmonic_indices,bₙvals,"bnvals")
    y = invFT(harmonic_indices,bₙvals,f₀)
    return y
end

npts = 150; nharmonics = 3
function plotIO(x::Function, f₀)
    t_range = LinRange(0, ΔT, npts)
    plot(t_range, DWPosition(x, ΔT, npts, nharmonics, T))
    plot!(t_range,x)
    title!("y(2,1)")
    xlabel!("Time (s)")
    ylabel!("Position (m)")
    savefig("img/vs_lattice/DW_Position.png")
end