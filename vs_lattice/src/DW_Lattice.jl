# General imports.
using LinearAlgebra
using Plots
using Measures
using DelimitedFiles
#using FFTW
gr()

# Specific imports.
include("Greens_Func.jl")

#using Kronecker
⊗(a,b) = kron(a,b)

# units
μm = 10^-6; GHz = 1; 

"""
DW Lattice struct w/ parameters
Made racetrack_positions and PBC mutable.
racetrack_positions is 2D vector, so that we can duplicate set of racetracks
"""
struct DWLattice
    n_racetracks::Ref{Int} # number of racetracks in the unit cell
    racetrack_positions::Ref{Vector{Float64}}
    orientations::Ref{Vector{Int}} # relative orientation of each domain wall;
    a::Float64 # lattice constant
    γ::Ref{Vector{Float64}} # damping constant for each oscillator
    ω₀::Ref{Vector{Float64}} # resonant frequency of each DW
    PBC::Ref{Bool} # periodic boundary conditions? for k-dependent Hamiltonian?
    C::Float64 # Stray-field coupling constant
    R₀::Float64 # distance between racetracks at which C was calculated
end

function init(; n_racetracks::Int, racetrack_positions::Vector{Float64}, orientations::Vector{Int},
        a::Float64, γ::Vector{Float64}, ω₀::Vector{Float64}, PBC::Bool, C::Float64, R₀::Float64)
    return DWLattice(n_racetracks, racetrack_positions, orientations, a, γ, ω₀, PBC, C, R₀);
end

testSystem = DWLattice(2, [0, 2.5]*μm, [1,-1], 5*μm, [0.1,0.1]*GHz, [10,15]*GHz, true, 10^8, 2.5*μm);

"""
Hamiltonian
"""
function HO_Hamiltonian(ω₀::Float64,γ::Float64)
    H₀ = [0 1; -ω₀^2 -2*γ]
    #display(H₀)
    return H₀
end

function Coupling_Hamiltonian(C::Union{ComplexF64, Float64})
    H = [0 0; C 0]; # adds a term dₜ²x₁ = -ω₀²x₁ - 2γdₜx₁ + C*x₂
    return H
end

function HO_Hamiltonian(ω₀::Float64,γ::Float64)
    H₀ = [0 1; -ω₀^2 -2*γ]
    #display(H₀)
    return H₀
end

function Coupling_Hamiltonian(C::Union{ComplexF64, Float64})
    H = [0 0; C 0]; # adds a term dₜ²x₁ = -ω₀²x₁ - 2γdₜx₁ + C*x₂
    return H
end

function constructHamiltonian(system::DWLattice, NNs::Int)
    function H(k::Union{ComplexF64,Float64})
        H₀ = zeros(ComplexF64,2*system.n_racetracks[],2*system.n_racetracks[])
        # add oscillators on diagonal
        for i = 1:system.n_racetracks[]
            ivec = zeros(system.n_racetracks[]); ivec[i] = 1;
            H₀ += ivec⊗(ivec')⊗HO_Hamiltonian(system.ω₀[][i],system.γ[][i]) 
        end
        # add coupling between oscillators
        for i = 1:system.n_racetracks[]
            for j = ((i-NNs):(i-1))∪((i+1):(i+NNs))
                if(j != i)
                    j_index = mod(j-1, system.n_racetracks[])+1
                    dR = (floor((j-1)/system.n_racetracks[]))*system.a # distance of nth unit cell away
                    ΔR = abs(system.racetrack_positions[][i] - (dR + system.racetrack_positions[][j_index])) # distance of nth racetrack
                    # make huge approximation here, say that Cij decays with 1/r² from calculated
                    # may need minus sign
                    Cij = system.orientations[][i][]*system.orientations[][j_index][]*system.C*(system.R₀/ΔR)^(2)
                    α_index = zeros(system.n_racetracks[],system.n_racetracks[]);
                    β_index = zeros(system.n_racetracks[],system.n_racetracks[]);
                    β_index[i,j_index] = -1; α_index[i,i] = 1
                    #println("i = $(i), j = $(j), J_index = $(j_index)")
                    if(j < 1 || j > system.n_racetracks[])
                        if(system.PBC[])
                            H₀ += β_index⊗Coupling_Hamiltonian(Cij*exp(im*k*dR))
                            H₀ += α_index⊗Coupling_Hamiltonian(Cij)
                            #println("k⋅dR/(2π) = $(k*dR/(2*π))")
                        end
                    else
                        H₀ += α_index⊗Coupling_Hamiltonian(Cij)
                        H₀ += β_index⊗Coupling_Hamiltonian(Cij)
                    end
                end
            end
        end
        #display(H₀)
        return H₀
    end
    return H
end


"""
Bands
"""
H_AFM_racetrack = constructHamiltonian(testSystem,2)

function plot1DBands(H::Function,system::DWLattice,nk::Int,Broadening::Bool=false, nE::Int = 200, DOS=true)
    a = system.a
    kvals = LinRange(-π/a,π/a,nk)
    ys = []
    xs = []
    maxE = maximum(imag.(eigvals(H(π/a)) ∪ eigvals(H(0.0))))
    if(Broadening)
        # get some reasonable maximum for the energies
        Evals = LinRange(0,maxE*1.1,nE)
        n = 2*system.n_racetracks[]
        bands = zeros(nk,nE)
        plottingBands = zeros(nk,nE)
        @Threads.threads for ik in eachindex(kvals)
            k = kvals[ik]
            Hatk = H(k)
            DOS_k(E) = (-1/π)imag.(tr(inv(E*I(n) - im*Hatk)))
            DOS_k_E = DOS_k.(Evals)
            bands[ik,:] = DOS_k_E
        end
        bands = (1/maximum(bands))*bands
        fig = heatmap(kvals*(a/π), Evals, bands', clims=(0,maximum(bands)), xlabel="k (π/a)", ylabel="Frequency (GHz)")
    end
    @Threads.threads for k in kvals
        Es = imag(eigvals(H(k)))
        append!(ys,Es/GHz)
        append!(xs,k*ones(size(Es)))
    end
    fig = scatter!(xs*(a/π),ys, legend=false,ylims=(-0.001,maxE*1.1), xlims=(-1,1), c="white", markersize=2.0, markerstrokewidth=0)
    if(DOS)
        DOS_tot = [sum(bands[:,iE]) for iE = 1:nE]; DOS_tot = (1/maximum(DOS_tot))*DOS_tot
        fDOS = plot(DOS_tot,Evals,legend=false,xlabel="DOS(ω)", ylabel="Frequency (GHz)",ylims=(-0.0001,maxE*1.1),xlims=(0,1))
    end
    fig2 = plot(fig,fDOS,layout=grid(1,2, widths=(5/8,3/8)), size=(800,300),margin=5mm)
    return fig2
end

function getBands(system::DWLattice,NNs::Int,broadening::Bool=false)
    H = constructHamiltonian(system,NNs)
    fig = plot1DBands(H,system,200,broadening)
    #display(fig)
end

"""
Making supercell.
Modify: n_racetracks, racetrack_positions, orientation, γ, ωₒ
(Ones that are vectors)
"""
function realSpaceArray(unit_cell::DWLattice, num::Int)
    # Shallow copy, creating alias
    orig_racetrack = unit_cell.racetrack_positions[]
    orig_orientations = unit_cell.orientations[]
    orig_γ = unit_cell.γ[]
    orig_ω₀ = unit_cell.ω₀[]

    # Duplicate unit_cell to make supercell.
    dup_lattice = deepcopy(unit_cell)
    dup_lattice.PBC[] = false

    for i = 1:(num-1)
        for j = 1:unit_cell.n_racetracks[]
            new_race = i*unit_cell.a + orig_racetrack[j]
            push!(dup_lattice.racetrack_positions[], new_race)
            push!(dup_lattice.orientations[], orig_orientations[j])
            push!(dup_lattice.γ[], orig_γ[j])
            push!(dup_lattice.ω₀[], orig_ω₀[j])
        end
    end

    # Increase n_racetracks accordingly
    dup_lattice.n_racetracks[] = num * unit_cell.n_racetracks[]

    return dup_lattice
end

"""
Returns transfer function
    @params
    Gʳ (Green's Function)
    race_i (racetrack being perturbed)
    race_j (racetrack being affected)
    startP: starting point of ω
    endP: ending point of ω
    step: # of points inbetween.

Includes the creation of supercell
"""
function transferFunc(system::DWLattice, n_lattice::Int, race_i::Int, race_j::Int, startP::Float64, endP::Float64, step::Int, NNs::Int=500, eval::Bool=true)
    # Configure the lattice
    new_sys = realSpaceArray(system, n_lattice)
    Gʳ = genGʳ_nok(constructHamiltonian(new_sys, NNs))

    # local function that returns transfer function
    function T(ω::Float64)
        Gʳ_val = Gʳ(ω)
        #writedlm("Greens_Func.txt", Gʳ_val)
        return Gʳ_val[2*race_i-1, 2*race_j-1] / Gʳ_val[2*race_i-1, 2*race_i-1]
    end
    if(eval==false)
        return T
    end
    # ω and corresponding transfer function value
    ω_total = LinRange(startP, endP, step)
    mag::Vector = (abs∘T).(ω_total)
    ϕ::Vector = (angle∘T).(ω_total)

    return T, ω_total, mag, ϕ
end