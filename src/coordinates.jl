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
    return EPRCoordinate(energy, pitch, r, z, t, amu*mass_u, q)
end

function EPRCoordinate(M::AxisymmetricEquilibrium, energy, pitch, R ; amu=H2_amu, q=1, dz=0.2)
    zaxis = M.axis[2]
    zmax = zaxis + dz
    zmin = zaxis - dz
    res = optimize(x->M.psi_rz(R,x), zmin, zmax)
    Z = Optim.minimizer(res)
    (Z == zmax || Z == zmin) && error(@sprintf("Unable to find starting Z value with dz = %.2f. Increase dz",dz))
    return EPRCoordinate(energy, pitch, R, Z, zero(Z), amu*mass_u, q)
end

function Base.show(io::IO, c::EPRCoordinate)
    println(io,typeof(c))
    @printf(io," E = %.3f keV\n",c.energy)
    @printf(io," pitch = %.3f\n",c.pitch)
    @printf(io," Rmax = %.3f m",c.r)
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
    return HamiltonianCoordinate(energy, mu, p_phi, amu*mass_u, q)
end

function HamiltonianCoordinate(M::AxisymmetricEquilibrium, KE, pitch, r, z, m, q)
    psi = M.psi_rz(r,z)
    babs = M.b(r,z)
    g = M.g(psi)
    pot = M.phi(psi)
    PE = 1e-3*pot

    mc2 = m*c0^2
    KE_j = e0*KE*1e3
    p_rel2 = ((KE_j + mc2)^2 - mc2^2)/(c0*c0)
    p_para = sqrt(p_rel2)*pitch*M.sigma
    p_perp2 = p_rel2*(1-pitch^2)

    mu = p_perp2/(2*m*babs)
    Pphi = -p_para*g/babs + q*e0*psi

    return HamiltonianCoordinate(KE+PE,mu,Pphi,m,q)
end

function HamiltonianCoordinate(M::AxisymmetricEquilibrium, c::EPRCoordinate)
    HamiltonianCoordinate(M,c.energy,c.pitch,c.r,c.z,c.m,c.q)
end

function HamiltonianCoordinate(M::AxisymmetricEquilibrium, c::HamiltonianCoordinate)
    return c
end

function HamiltonianCoordinate(M::AxisymmetricEquilibrium, KE, p, r, z; amu=H2_amu, q=1)
    HamiltonianCoordinate(M, KE, p, r, z, amu*mass_u, q)
end

function HamiltonianCoordinate(M::AxisymmetricEquilibrium, c::GCParticle)
    HamiltonianCoordinate(M,c.energy,c.pitch,c.r,c.z,c.m,c.q)
end

function Base.show(io::IO, c::HamiltonianCoordinate)
    println(io,typeof(c))
    @printf(io," E = %.3f keV\n",c.energy)
    @printf(io," μ₀ = %.3e\n",c.mu)
    @printf(io," Pᵩ = %.3e",c.p_phi)
end
