


# two ways of calculating the retarded green's function in code! 
# (I like the phrase "time-retarded" but the internet may be more helpful without saying "time")

# 1: Immediately spit out Gʳ::Matrix{ComplexF64} values of a specific Hamiltonian::Matrix{ComplexF64} at a specific k vector and energy
# 2: Or, after defining a H::Function (k::Vector{Float64}) in a given scope, can define Gʳ::Function as a function of ω and k. 

# The latter is more useful if you get in the "functional" mindset. 
# it makes more sense in my head but this might be infurating idk


# APPROACH 1:
function GʳofH(H::Union{Matrix{ComplexF64},Matrix{Float64}}, ω::Union{Float64, ComplexF64})
	return inv(ω*I - im*H)
end

# APPROACH 2:
# Okay this is a cool feature of julia
# We will define a function that spits out a function!
# Ultimately, this will return a function Gʳ(k,ω) "the retarded green's function"
function genGʳ(H::Function)
	function Gʳ(ω::Union{Float64, ComplexF64}, k::Vector{Float64})
		return inv(ω*I - im*H(k))
	end
	return Gʳ
end
