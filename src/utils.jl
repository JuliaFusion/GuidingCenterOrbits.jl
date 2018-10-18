function get_pitch(M::AxisymmetricEquilibrium, c::HamiltonianCoordinate, r::T, z::T) where {T<:Real}
    psi = M.psi_rz(r,z)
    g = M.g(psi)
    babs = M.b(r,z)
    phi = M.phi(psi)
    KE = c.energy - 1e-3*phi
    f = -babs/(sqrt(2e3*e0*KE*c.m)*g*M.sigma)
    pitch = f*(c.p_phi - c.q*e0*psi)
    #pitchabs = sqrt(max(1.0-(c.mu*babs/(1e3*e0*KE)), 0.0))
    #if !isapprox(abs(pitch), pitchabs, atol=1.e-1)
    #    warn("abs(pitch) != abspitch: ",pitchabs," ",pitch)
    #end
    return clamp(pitch,-1.0,1.0)
end

function get_pitch(M::AxisymmetricEquilibrium, c::T, r::Vector{S}, z::Vector{S}) where {T<:AbstractOrbitCoordinate, S<:Real}
    n = length(r)
    pitch = zeros(S,n)
    @inbounds for i=1:n
        pitch[i] = get_pitch(M, c, r[i], z[i])
    end
    return pitch
end

function get_pitch(M::AxisymmetricEquilibrium, c::T, path::OrbitPath) where {T<:AbstractOrbitCoordinate}
    return get_pitch(M, c, path.r, path.z)
end

function get_pitch(M::AxisymmetricEquilibrium, c::T, r::S, z::S) where {T<:AbstractOrbitCoordinate, S<:Real}
    hc = HamiltonianCoordinate(M, c)
    return get_pitch(M, hc, r, z)
end

function get_pitch(M::AxisymmetricEquilibrium, o::Orbit)
    path = o.path
    n = length(path.r)
    pitch = zeros(eltype(path.r),n)
    @inbounds for i=1:n
        pitch[i] = get_pitch(M, o.coordinate, path.r[i], path.z[i])
    end
    return pitch
end

function get_kinetic_energy(M::AxisymmetricEquilibrium, c::HamiltonianCoordinate, r::T, z::T) where {T<:Real}
    psi = M.psi_rz(r,z)
    return c.energy - 1e-3*M.phi(psi)
end

function get_kinetic_energy(M::AxisymmetricEquilibrium, c::T, r::Vector{S}, z::Vector{S}) where {S, T<:AbstractOrbitCoordinate}
    n = length(r)
    energy = zeros(S,n)
    @inbounds for i=1:n
        energy[i] = get_kinetic_energy(M, c, r[i], z[i])
    end
    return energy
end

function get_kinetic_energy(M::AxisymmetricEquilibrium, c::T, path::OrbitPath) where {T<:AbstractOrbitCoordinate}
    return get_kinetic_energy(M, c, path.r, path.z)
end

function get_kinetic_energy(M::AxisymmetricEquilibrium, c::T, r::S, z::S) where {T<:AbstractOrbitCoordinate, S<:Real}
    hc = HamiltonianCoordinate(M, c)
    return get_kinetic_energy(M, hc, r, z)
end

function get_kinetic_energy(M::AxisymmetricEquilibrium, o::Orbit)
    path = o.path
    c = o.coordinate
    return get_kinetic_energy(M, o.coordinate, o.path)
end

function classify(path::OrbitPath, axis; n=length(path))

    pp = take(zip(path.r,path.z),n)
    op = Limiter()
    for p in pp
        push!(op.vertices,p)
    end
    push!(op.vertices,first(pp))

    sign_p = sign.(path.pitch)
    if in_vessel(op, axis)
        if all(sign_p .== sign_p[1])
            if sign_p[1] > 0.0
                class = :co_passing
            else
                class = :ctr_passing
            end
        else
            class = :potato
        end
    else
        if all(sign_p .== sign_p[1])
            class = :stagnation
        else
            class = :trapped
        end
    end

    return class
end

function hits_wall_path(path::OrbitPath, wall::Limiter)

    not_in_vessel = [~in_vessel(wall,p) for p in zip(path.r,path.z)]
    ind = findfirst(not_in_vessel)
    if ind != nothing
        return false, path
    else
        ind = ind-1
        return true, OrbitPath(path.r[1:ind], path.z[1:ind] ,path.phi[1:ind],
                               path.pitch[1:ind], path.energy[1:ind],
                               path.dt[1:ind], path.dl[1:ind])
    end
end

function hits_wall(path::OrbitPath, wall::Limiter)
    not_in_vessel = [~in_vessel(wall,p) for p in zip(path.r,path.z)]
    return any(not_in_vessel)
end

function ion_cyclotron_period(M,p::AbstractParticle)
    return 2*pi*(p.m)/(M.b(p.r,p.z)*e0) #Ion Cyclotron Period
end
