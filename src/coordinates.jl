abstract type AbstractOrbitCoordinate{T} end

#Enable Broadcasting
Base.broadcastable(x::T) where T <: AbstractOrbitCoordinate = (x,)

struct EPRCoordinate{T} <: AbstractOrbitCoordinate{T}
    energy::T # Kinetic Energy
    pitch::T
    r::T
    z::T
    t::T
    m::Float64
    q::Int
end

function EPRCoordinate(T=Float64; amu=H2_amu, q = 1)
    z = zero(T)
    return EPRCoordinate(z, z, z, z, z, amu*mass_u, q)
end

function EPRCoordinate(energy, pitch, r, z; t = zero(z), amu=H2_amu, q = 1)
    return EPRCoordinate(promote(energy, pitch, r, z, t)..., amu*mass_u, q)
end

function EPRCoordinate(M::AbstractEquilibrium, energy, pitch, R ; amu=H2_amu, q=1, dz=0.2)
    rmaxis, zaxis = magnetic_axis(M)
    zmax = zaxis + dz
    zmin = zaxis - dz
    if (M.psi[2]-M.psi[1]) > 0.0
        res = optimize(x->M(R,x), zmin, zmax) # Minimize
    else
        res = optimize(x->-M(R,x), zmin, zmax) # Maximize
    end
    Z = Optim.minimizer(res)
    (Z == zmax || Z == zmin) && error(@sprintf("Unable to find starting Z value with dz = %.2f. Increase dz",dz))
    return EPRCoordinate(promote(energy, pitch, R, Z, 0)..., amu*mass_u, q)
end

function GCParticle(c::EPRCoordinate)
    GCParticle(c.energy,c.pitch,c.r,c.z,c.m,c.q)
end

function EPRParticle(oc::EPRCoordinate; t=oc.t)
    EPRParticle(oc.energy,oc.pitch,oc.r,oc.z,t,oc.m,oc.q)
end

function Base.show(io::IO, c::EPRCoordinate)
    println(io,typeof(c))
    @printf(io," E = %.3f keV\n",c.energy)
    @printf(io," pitch = %.3f\n",c.pitch)
    @printf(io," Rmax = %.3f m",c.r)
end

#This will fill out orbits in the PSGrid, but not the orbit grid
#I want it to just be 
struct GCEPRCoordinate
    energy::Float64
    pitch::Float64
    r::Float64
    z::Float64
    pitch_m::Float64
    r_m::Float64
    z_m::Float64
    t::Float64
    class::Char
    tau_p::Float64
    tau_t::Float64
    jacdet::Float64
    gcvalid::Bool
    m::Float64
    q::Int
end

function GCEPRCoordinate(args...; jacdet::Float64=0.0, gcvalid::Bool=false, m::Float64=H2_amu*mass_u, q::Int=1) 
    GCEPRCoordinate(promote(args)..., jacdet, gcvalid, m, q)
end

function GCEPRCoordinate(gcp::GCParticle,class::Char='i'; gcvalid::Bool=false) 
    return GCEPRCoordinate(gcp.energy,gcp.pitch,gcp.r,gcp.z,0.0,0.0,0.0,0.0,class,0.0,0.0,0.0,gcvalid,gcp.m,gcp.q)
end

struct HamiltonianCoordinate{T} <: AbstractOrbitCoordinate{T}
    energy::T #Kinetic + Potential Energy
    mu::T
    p_phi::T
    m::Float64
    q::Int
end

function HamiltonianCoordinate(T=Float64; amu=H2_amu, q=1)
    z = zero(T)
    return HamiltonianCoordinate(z, z, z, amu*mass_u, q)
end

function HamiltonianCoordinate(energy, mu, p_phi; amu=H2_amu, q=1)
    return HamiltonianCoordinate(promote(energy, mu, p_phi)..., amu*mass_u, q)
end

function HamiltonianCoordinate(M::AbstractEquilibrium, KE, pitch, r, z, m, q)
    psi = M(r,z)
    babs = norm(Bfield(M,r,z))
    g = poloidal_current(M,psi)
    pot = electric_potential(M,psi)
    PE = 1e-3*pot

    mc2 = m*c0^2
    KE_j = e0*KE*1e3
    p_rel2 = ((KE_j + mc2)^2 - mc2^2)/(c0*c0)
    p_para = sqrt(p_rel2)*pitch*B0Ip_sign(M)
    p_perp2 = p_rel2*(1-pitch^2)

    mu = p_perp2/(2*m*babs)
    Pphi = p_para*g/babs + q*e0*psi

    return HamiltonianCoordinate(promote(KE+PE,mu,Pphi)...,m,q)
end

function HamiltonianCoordinate(M::AbstractEquilibrium, c::EPRCoordinate)
    HamiltonianCoordinate(M,c.energy,c.pitch,c.r,c.z,c.m,c.q)
end

function HamiltonianCoordinate(M::AbstractEquilibrium, c::HamiltonianCoordinate)
    return c
end

function HamiltonianCoordinate(M::AbstractEquilibrium, KE, p, r, z; amu=H2_amu, q=1)
    HamiltonianCoordinate(M, KE, p, r, z, amu*mass_u, q)
end

function HamiltonianCoordinate(M::AbstractEquilibrium, c::GCParticle)
    HamiltonianCoordinate(M,c.energy,c.pitch,c.r,c.z,c.m,c.q)
end

function Base.show(io::IO, c::HamiltonianCoordinate)
    println(io,typeof(c))
    @printf(io," E = %.3f keV\n",c.energy)
    @printf(io," μ₀ = %.3e\n",c.mu)
    @printf(io," Pᵩ = %.3e",c.p_phi)
end
