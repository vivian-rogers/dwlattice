# General imports.
using LinearAlgebra
using Plots
using Measures
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
mutable struct DWLattice
    n_racetracks::Int # number of racetracks in the unit cell
    racetrack_positions::Vector{Float64}
    orientations::Vector{Int} # relative orientation of each domain wall;
    a::Float64 # lattice constant
    γ::Vector{Float64} # damping constant for each oscillator
    ω₀::Vector{Float64} # resonant frequency of each DW
    PBC::Bool # periodic boundary conditions? for k-dependent Hamiltonian?
    C::Float64 # Stray-field coupling constant. Is Fext/(m*Hext)
    #R₀::Float64 # distance between racetracks at which C was calculated
    w_RT::Float64 # width of racetrack in meters
    Ms::Float64 # magnetic saturation value of material
    t_RT::Float64 # thickness of racetrack in meters
    α::Float64 # gilbert damping
    m::Float64 # effective mass of domain wall
end

nm=10^-9
μ₀ = 9.274E-24 # Bohr Magneton in Joules/Tesla
γ₀ = 1.7609E11 # electron gyromagnetic ratio in rad/(second*Tesla)
fudge = 1E0
"""
Creates a new system struct 
"""
function init(; n_racetracks::Int, racetrack_positions::Vector{Float64}, orientations::Vector{Int},
        a::Float64, ω₀::Vector{Float64}, PBC::Bool, w_RT::Float64, Ms::Float64, t::Float64, α::Float64, w_DW::Float64=40*nm)
        NzminNy = 1;
        γ = α*Ms*γ₀*NzminNy/(2*(1+α^2))*fudge
        m = 2*(1+α^2)*μ₀*(t*w_RT)/(γ₀*w_DW*NzminNy)
        C = 2*μ₀*w_RT*t*Ms/m
        return DWLattice(n_racetracks, racetrack_positions, orientations, a, γ*ones(n_racetracks), ω₀, PBC, C, w_RT, Ms, t, α, m);
end

#function init(; n_racetracks::Int, racetrack_positions::Vector{Float64}, orientations::Vector{Int},
#        a::Float64, γ::Vector{Float64}, ω₀::Vector{Float64}, PBC::Bool, C::Float64, R₀::Float64)
#    return DWLattice(n_racetracks, racetrack_positions, orientations, a, γ, ω₀, PBC, C, R₀);
#end

#testSystem = DWLattice(2, [0, 2.5]*μm, [1,-1], 5*μm, [0.1,0.1]*GHz, [10,15]*GHz, true, 10^8, 2.5*μm);

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

function ∂Hz_∂x(system::DWLattice,ΔR::Float64) # taylor series in x approximation for stray fields
    Λplus = 2*(system.w_RT-2*ΔR)/(system.t_RT*√( system.t_RT^2 + system.w_RT^2 - 4*system.w_RT*ΔR + 4*ΔR^2 ) )
    Λmin = 2*(system.w_RT+2*ΔR)/(system.t_RT*√( system.t_RT^2 + system.w_RT^2 + 4*system.w_RT*ΔR + 4*ΔR^2 ) )
    return -(system.Ms/π)*(Λmin + Λplus)
end

function constructHamiltonian(system::DWLattice, NNs::Int)
    function H(k::Union{ComplexF64,Float64})
        H₀ = zeros(ComplexF64,2*system.n_racetracks,2*system.n_racetracks)
        # add oscillators on diagonal
        for i = 1:system.n_racetracks
            ivec = zeros(system.n_racetracks); ivec[i] = 1;
            H₀ += ivec⊗(ivec')⊗HO_Hamiltonian(system.ω₀[i],system.γ[i]) 
        end
        # add coupling between oscillators
        for i = 1:system.n_racetracks
            for j = ((i-NNs):(i-1))∪((i+1):(i+NNs))
                if(j != i)
                    j_index = mod(j-1, system.n_racetracks)+1
                    dR = (floor((j-1)/system.n_racetracks))*system.a # distance of nth unit cell away
                    ΔR = abs(system.racetrack_positions[i] - (dR + system.racetrack_positions[j_index])) # distance of nth racetrack
                    # make huge approximation here, say that Cij decays with 1/r² from calculated
                    # may need minus sign
                    Cij = system.orientations[i]*system.orientations[j_index]*system.C*∂Hz_∂x(system,ΔR)
                    α_index = zeros(system.n_racetracks,system.n_racetracks);
                    β_index = zeros(system.n_racetracks,system.n_racetracks);
                    β_index[i,j_index] = -1; α_index[i,i] = 1
                    #println("i = $(i), j = $(j), J_index = $(j_index)")
                    if(j < 1 || j > system.n_racetracks)
                        if(system.PBC)
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
        return (1/GHz)*H₀
    end
    return H
end

"""
Making supercell.
Modify: n_racetracks, racetrack_positions, orientation, γ, ωₒ
(Ones that are vectors)
"""
function realSpaceArray(unit_cell::DWLattice, num::Int)
    # Shallow copy, creating alias
    orig_racetrack = unit_cell.racetrack_positions
    orig_orientations = unit_cell.orientations
    orig_γ = unit_cell.γ
    orig_ω₀ = unit_cell.ω₀

    # Duplicate unit_cell to make supercell.
    dup_lattice = deepcopy(unit_cell)
    dup_lattice.PBC = false

    for i = 1:(num-1)
        for j = 1:unit_cell.n_racetracks
            new_race = i*unit_cell.a + orig_racetrack[j]
            push!(dup_lattice.racetrack_positions, new_race)
            push!(dup_lattice.orientations, orig_orientations[j])
            push!(dup_lattice.γ, orig_γ[j])
            push!(dup_lattice.ω₀, orig_ω₀[j])
        end
    end

    # Increase n_racetracks accordingly
    dup_lattice.n_racetracks = num * unit_cell.n_racetracks

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
function transferFunc(system::DWLattice, Ncells::Int, race_i::Int, race_j::Int, startP::Float64, endP::Float64, step::Int; NNs::Int=500, eval::Bool=true)
    # Configure the lattice
    new_sys = realSpaceArray(system, Ncells)
    Gʳ = genGʳ_nok(constructHamiltonian(new_sys, NNs))

    # local function that returns transfer function
    function T(ω::Float64)
        Gʳ_val = Gʳ(ω)
        if(ω≈0)
            return 0
        end
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