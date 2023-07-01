using LinearAlgebra
using SparseArrays
using Plots


nm = 1E-9; ε₀ = 8.854E-12; q=1.602E-19 
#nm = 1.0; ε₀ = 1.0;
function getPos(p::NamedTuple, ivec::Vector{Int})
    return (ivec.*[p.Lx/p.Nx,p.Lz/p.Nz])
end

function bounds(p::NamedTuple, ivec::Vector{Int})
    ix = ivec[1]; iz = ivec[2]
    return [mod(ix, p.Nx), mod(iz,p.Nz)]
end

function index(p::NamedTuple,ivec::Vector{Int})
    return p.Nz*ivec[1] + ivec[2] +1
end

function mat_at_index(p::NamedTuple,ivec::Vector{Int})
    R = getPos(p,ivec)
    x = R[1]; z = R[2]
    if z < p.z_substrate
        return "substrate"
    elseif z < (p.z_substrate + p.t_HM)
        return "HM"
    elseif z < (p.z_substrate + p.t_HM + p.t_racetrack)
        return "RT"
    elseif abs(x-p.Lx/2) > p.w_gate/2
        return "ins"
    elseif z < (p.z_substrate + p.t_HM + p.t_racetrack + p.t_gate)
        return "gatedielectric"
    else
        return "gatecontact"
    end
end
function ε_at_index(p::NamedTuple,ivec::Vector{Int})
    R = getPos(p,ivec)
    x = R[1]; z = R[2]
    if z < p.z_substrate
        return p.ε_substrate
    elseif z < (p.z_substrate + p.t_HM)
        return p.ε_HM
    elseif z < (p.z_substrate + p.t_HM + p.t_racetrack)
        return p.ε_RT
    elseif abs(x-p.Lx/2) > p.w_gate/2
        return p.ε_ins
    elseif z < (p.z_substrate + p.t_HM + p.t_racetrack + p.t_gate)
        return p.ε_gate
    else
        return p.ε_RT
    end
end

function inlattice(p::NamedTuple,ivec::Vector{Int})
    if(ivec[1] < p.Nx && ivec[1] >= 0 && ivec[2] < p.Nz && ivec[2] >= 0)
        return true
    else
        return false
    end
end

function ε_vec(p::NamedTuple)
    εvals = zeros(p.Nx*p.Nz)
    for ix_site = 0:(p.Nx-1)
        for iz_site = 0:(p.Nz-1)
            ivec = [ix_site,iz_site]
            εvals[index(p,ivec)] = ε_at_index(p,ivec)
        end
    end
    return εvals
end

function define_∇(p::NamedTuple)
    M = zeros(p.Nx*p.Nz,p.Nx*p.Nz)
    δx = p.Lx/p.Nx; δz = p.Lz/p.Nz; 
    for ix_site = 0:(p.Nx-1)
        for iz_site = 0:(p.Nz-1)
            iAvec = [ix_site,iz_site]
            #ε = ε_at_index(p,iAvec)
            iBvec = deepcopy(iAvec)
            M[index(p,iAvec),index(p,iBvec)] = -2/δx^2 -2/δz^2 
            for Δivec ∈ [[1,0],[-1,0]]
                iBvec = iAvec + Δivec
                iBvec = bounds(p,iBvec)
                #if(inlattice(p,iBvec))
                M[index(p,iAvec),index(p,iBvec)] = 1/δx^2
                #end
            end
            for Δivec ∈ [[0,-1],[0,1]]
                iBvec = iAvec + Δivec
                if(inlattice(p,iBvec))
                    M[index(p,iAvec),index(p,iBvec)] = 1/δz^2
                end
            end

        end
    end
    return M
end

function fixedVoltages(p::NamedTuple,ΔV::Float64)
    fixedVs = zeros(p.Nx*p.Nz)
    for ix_site = 0:(p.Nx-1)
        for iz_site = 0:(p.Nz-1)
            ivec = [ix_site,iz_site]
            material = mat_at_index(p,ivec)
            if(material == "RT" || material == "HM")
                fixedVs[index(p,ivec)] = -ΔV/2
            elseif(material =="gatecontact")
                fixedVs[index(p,ivec)] = ΔV/2
            end
        end
    end
    # bit complicated. Basically we want to discard rows where we say 0*V = 0
    fix_flag = diagm((!).(iszero.(fixedVs))) # vector of flags denoting fixed voltages 
    extendA = fix_flag[vec(mapslices(col -> any(col .!= 0), fix_flag, dims = 2)), :]
    filter!(V->V≠0,fixedVs)
    return extendA, fixedVs
end


function fixedCharges(p::NamedTuple)
    fix_flag = zeros(p.Nx*p.Nz)
    n_fixed = 0
    for ix_site = 0:(p.Nx-1)
        for iz_site = 0:(p.Nz-1)
            ivec = [ix_site,iz_site]
            material = mat_at_index(p,ivec)
            fix_flag[index(p,ivec)] = 0; n_fixed +=0
            #if(material == "substrate" || material == "ins" || material == "gatedielectric")
            #    fix_flag[index(p,ivec)] = 0
            #    n_fixed += 1
            #end
        end
    end
    # get a vector with true if the element is not zero
    # this is perhaps a bit jank. Should be diagonal, with the presence of 
    # a "1" designating that the charge density will be fixed to 0. 
    fixedρs = zeros(n_fixed)
    extendA = diagm(fix_flag)
    extendA = extendA[vec(mapslices(col -> any(col .!= 0), extendA, dims = 2)), :]
    return extendA, fixedρs
end

function performPoisson(p::NamedTuple,ΔV::Float64)
    # wish to set up problem as 
    # A*x = b. First define A
    ∇² = define_∇(p) # the laplacian times dielectric const
    εmat = diagm(ε_vec(p)) #diagonal with permittivities
    extendA, fixedVs = fixedVoltages(p,ΔV) # vector of voltage values which will be fixed
    A = sparse(vcat(-εmat*∇²,extendA)) # set up laplacian with boundary conditions
    #V_to_fix = diagm((!).(iszero.(fixedVs))) # vector of flags denoting fixed voltages 
    #heatmap((log10∘abs).(A.+1E-10)',dpi=600)
    #savefig("./micromagnetics/matrix.png")
    # then define b, the charge densities and fixed voltages
    ρvals = zeros(p.Nx*p.Nz)
    b = vcat(ρvals,fixedVs)

    # perfect. Now we have a bad guess for the voltage values. 
    Vvals = A \ b
    # let's now fixed the charge density boundary conditions
    inv_∇² = inv(∇²)
    extendA_charges, fixedρs = fixedCharges(p)
    A_charges = sparse(vcat(-inv_∇²*inv(εmat),extendA_charges))
    for ip = 1:p.n_iter
        # and let's just loop this for a while
        Vprev = deepcopy(Vvals)
        b_charges = vcat(Vvals, fixedρs)
        ρvals = A_charges \ b_charges
        b = vcat(ρvals,fixedVs)
        Vvals = A \ b
        error = norm(Vprev - Vvals)/norm(Vvals)
        println("iteration $(ip), error = $(error)")
    end
    return Vvals, ρvals
end

function plotVoltageCharge(p::NamedTuple,Vvals::Vector{Float64},ρvals::Vector{Float64})
    Vmat = reshape(Vvals,(p.Nz,p.Nx))
    ρmat = reshape(ρvals,(p.Nz,p.Nx))
    zvals = collect(LinRange(0,p.Lz,p.Nz))
    xvals = collect(LinRange(0,p.Lx,p.Nx))
    fig1 =  heatmap(xvals,zvals,Vmat,cmap=:bwr)
    savefig(fig1,"./micromagnetics/Voltage.png")
    maxρ = maximum(ρmat)
    fig2 =  heatmap(xvals,zvals,ρmat,cmap=:bwr,clim=(-maxρ,maxρ))
    savefig(fig2,"./micromagnetics/Charge_density.png")
end


params = (n_iter = 40, Nx = 50, Nz = 50, Lz = 12*nm, Lx = 100*nm,
    z_substrate=2*nm, t_HM=1*nm, t_racetrack=2*nm, t_gate=3*nm, w_gate=40*nm, 
    ε_gate = 25*ε₀, ε_ins = 9*ε₀, ε_substrate=5*ε₀, ε_RT = 1E4*ε₀, ε_HM =1E3*ε₀)

Vvals, ρvals = performPoisson(params,0.25)
plotVoltageCharge(params,Vvals,ρvals)