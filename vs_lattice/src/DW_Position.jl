# Specific Imports.
include("DW_Lattice.jl")
include("Fourier.jl")

# units
μm = 10^-6; GHz = 10^9; 

testSystem = DWLattice(2, [0, 2.5]*μm, [1,-1], 5*μm, [0.1,0.1]*GHz, [10,15]*GHz, true, 10^8, 2.5*μm);
H_AFM_racetrack = constructHamiltonian(testSystem,2)

AFM2R_lattice = init(n_racetracks=2, racetrack_positions=[0, 2.6]*μm, orientations=[1,-1], 
    a=5*μm, γ=0.5*[1,1]*GHz, ω₀=1.0*[9,8]*GHz, PBC=true, C=10.0*GHz, R₀=2.5*μm);
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
start_f = 0.0000001*GHz
stop_f = 15.0*GHz

function genFiring_neuron_DWpos(f::Float64)
    pos(t) = f*mod(t,1/f) -1/2
    return pos
end

ΔT = 1/(5*GHz)
x = genFiring_neuron_DWpos(1/ΔT)

T, ω_val, mag, ϕ = transferFunc(AFM2R_lattice, 5, 1, 2, start_f, stop_f, 200)

function DWPosition(x::Function, Δt::Float64, N::Int, total::Float64, T::Function)
    # Array of harmonic #s and corresponding Magnitude 
    f₀
    range, aₙvals = Mag_Harmonics(x,Δt, N, total)
    ωvals = (2π/Δt).*range
    bₙvals = aₙvals.*T.(ωvals)
    y = invFT(range,bₙvals,f₀)
    # range in frequency.
    #ω_val = zeros(convert(Int, total))
    #@. ω_val = range /Δt
    # Array of Transfer Function magnitudes for corresponding freqs.
    #T_mag = zeros(convert(Int, total))
    #@. T_mag = (abs∘T).(ω_val)
    # Property: Inverse Fourier is sum of Fourier * Fourier (Magnitude)
    #y = zeros(convert(Int, total))
    #@. y = aₙ * T_mag
    
    return y
end

t_range = LinRange(0, ΔT, 150)
plot(t_range, DWPosition(x, ΔT, 50, 150.0, T))
plot!(t_range,x)
title!("Position of DW of racetrack 10 when oscillating racetrack 1")
xlabel!("Time (s)")
ylabel!("Position (m)")
savefig("img/DW_Position.png")
