abstract type AbstractParticle{T} end

struct Particle{T} <: AbstractParticle{T}
    r::T
    phi::T
    z::T
    vr::T
    vt::T
    vz::T
    m::Float64
    q::Int
end

function Electron(r,phi,z,vr,vt,vz)
    Particle(r,phi,z,vr,vt,vz,e_amu*mass_u,-1)
end

function Proton(r,phi,z,vr,vt,vz)
    Particle(r,phi,z,vr,vt,vz,H1_amu*mass_u,1)
end

function Deuteron(r,phi,z,vr,vt,vz)
    Particle(r,phi,z,vr,vt,vz,H2_amu*mass_u,1)
end

function Triton(r,phi,z,vr,vt,vz)
    Particle(r,phi,z,vr,vt,vz,H3_amu*mass_u,1)
end

function Triton(r,phi,z,vr,vt,vz)
    Particle(r,phi,z,vr,vt,vz,He3_amu*mass_u,2)
end

struct GCParticle{T} <: AbstractParticle{T}
    energy::T
    pitch::T
    r::T
    z::T
    m::Float64
    q::Int
end

function GCParticle(energy,pitch,r,z; amu=H2_amu, q=1)
    GCParticle(energy, pitch, r, z, mass_u*amu, q)
end

function GCParticle(c::EPRCoordinate)
    GCParticle(c.energy,c.pitch,c.r,c.z,c.m,c.q)
end

function GCElectron(energy,pitch,r,z)
    GCParticle(energy,pitch,r,z,e_amu*mass_u,-1)
end

function GCProton(energy,pitch,r,z)
    GCParticle(energy,pitch,r,z,H1_amu*mass_u,1)
end

function GCDeuteron(energy,pitch,r,z)
    GCParticle(energy,pitch,r,z,H2_amu*mass_u,1)
end

function GCTriton(energy,pitch,r,z)
    GCParticle(energy,pitch,r,z,H3_amu*mass_u,1)
end

function GCAlpha(energy,pitch,r,z)
    GCParticle(energy,pitch,r,z,He3_amu*mass_u,2)
end

struct EPRParticle{T} <: AbstractParticle{T}
    energy_c::T
    pitch_c::T
    r_c::T
    z_c::T
    t::T
    m::Float64
    q::Int
end

function EPRParticle(energy, pitch, r, z, t; amu=H2_amu, q=1)
    EPRParticle(energy, pitch, r, z, t, mass_u*amu, q)
end

function EPRParticle(oc::EPRCoordinate; t=oc.t)
    EPRParticle(oc.energy,oc.pitch,oc.r,oc.z,t,oc.m,oc.q)
end

#Enable Broadcasting
Base.broadcastable(x::T) where T <: AbstractParticle = (x,)
