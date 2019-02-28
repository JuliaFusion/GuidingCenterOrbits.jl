struct RelativisticOrbitPath{T}
    ppar::Vector{T}
    pperp::Vector{T}
    r::Vector{T}
    z::Vector{T}
    phi::Vector{T}
    dt::Union{T,Vector{T}}
end
Base.length(op::RelativisticOrbitPath) = length(op.r)

mutable struct RelativisticOrbitStatus
    errcode::Int
end
RelativisticOrbitStatus() = RelativisticOrbitStatus(1)

struct RelativisticCoordinate
    r
    z
    ppar
    pperp
    mu
end

struct RelativisticGCParticle
    r
    z
    ppar
    pperp
    mu
end

function relativistic_gc_velocity(M::AxisymmetricEquilibrium)

end
