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

function Electron(args...)
    Particle(promote(args...)...,e_amu*mass_u,-1)
end

function Proton(args...)
    Particle(promote(args...)...,H1_amu*mass_u,1)
end

function Deuteron(args...)
    Particle(promote(args...)...,H2_amu*mass_u,1)
end

function Triton(args...)
    Particle(promote(args...)...,H3_amu*mass_u,1)
end

function Alpha(args...)
    Particle(promote(args...)...,He3_amu*mass_u,2)
end

struct GCParticle{T} <: AbstractParticle{T}
    energy::T
    pitch::T
    r::T
    z::T
    m::Float64
    q::Int
end

function GCElectron(args...)
    GCParticle(promote(args...)...,e_amu*mass_u,-1)
end

function GCProton(args...)
    GCParticle(promote(args...)...,H1_amu*mass_u,1)
end

function GCDeuteron(args...)
    GCParticle(promote(args...)...,H2_amu*mass_u,1)
end

function GCTriton(args...)
    GCParticle(promote(args...)...,H3_amu*mass_u,1)
end

function GCAlpha(args...)
    GCParticle(promote(args...)...,He3_amu*mass_u,2)
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

function EPRParticle(args...; amu=H2_amu, q=1)
    EPRParticle(promote(args)..., mass_u*amu, q)
end

#Enable Broadcasting
Base.broadcastable(x::T) where T <: AbstractParticle = (x,)
