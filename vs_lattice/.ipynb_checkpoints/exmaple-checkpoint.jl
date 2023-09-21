using Plots
using FFTW
using DSP
include("./src/DW_Lattice.jl")

# Number of points 
N = 100
# Sample period
Ts = 1 / N
# Start time 
t0 = 0 
tmax = t0 + N * Ts
# time coordinate
# t = t0:Ts:tmax
t = range(t0, stop=tmax, step=Ts)
# cycles/sec
f = 1

# signal 
signal = cos.(2π * f .* t) # sin (2π f t) 

# Fourier Transform of it 
F = fft(signal) |> fftshift
freqs = fftfreq(length(t), 1.0/Ts) |> fftshift

# Lattice structure
testSystem = DWLattice(2, [0, 2.5]*μm, [1,-1], 5*μm, [0.1,0.1]*GHz, [10,15]*GHz, true, 10^8, 2.5*μm);
H_AFM_racetrack = constructHamiltonian(testSystem,2)

AFM2R_lattice = init(n_racetracks=2, racetrack_positions=[0, 2.6]*μm, orientations=[1,-1], 
    a=5*μm, γ=0.5*[1,1]*GHz, ω₀=1.0*[9,8]*GHz, PBC=true, C=10.0*GHz, R₀=2.5*μm);
AFM3R_lattice = init(n_racetracks=3, racetrack_positions=[0, 2.8, 5.0]*μm, orientations=[1,-1,1], 
    a=7.5*μm, γ=0.5*[1,0.5,1]*GHz, ω₀=10*[0,2.0,2.1]*GHz, PBC=true, C=400.0*GHz, R₀=2.5*μm);
FM1R_lattice = init(n_racetracks=1, racetrack_positions=[0]*μm, orientations=[1], 
    a=5*μm, γ=0.5*[1]*GHz, ω₀=10*[2.0]*GHz, PBC=true, C=80.0*GHz, R₀=2.5*μm);

T, ω_val, mag, ϕ = transferFunc(AFM2R_lattice, 5, 1, 10, 0.0000001*GHz, 15.0*GHz, 101)

# basic plots 
# time_domain = plot(t, signal, title = "Signal")
# freq_mag = plot(ω_val, abs.(F), title = "Mang. Spectrum", xlim=(0, 20)) 
# plot(time_domain, freq_mag, layout=2)
# savefig("Wave.png")

# With Transfer Function
y = ifft(F.*(T.(ω_val)))
mag_sp = plot(ω_val, abs.(F), title = "Mag. Spectrum")
transfer = plot(ω_val, abs.(T.(ω_val)), title = "Transfer Function")
position = plot(ω_val, abs.(y), title = "Position Change")
plot(mag_sp, transfer, position, layout = 3)
savefig("Position.png")