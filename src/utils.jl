function get_pitch{T<:Real}(M::AxisymmetricEquilibrium, c::HamiltonianCoordinate, r::T, z::T)
    psi = M.psi_rz[r,z]
    g = M.g[psi]
    babs = M.b[r,z]
    phi = M.phi[psi]
    KE = c.energy - 1e-3*phi
    f = -babs/(sqrt(2e3*e0*KE*mass_u*c.amu)*g*M.sigma)
    pitch = f*(c.p_phi - c.q*e0*psi)
    #pitchabs = sqrt(max(1.0-(c.mu*babs/(1e3*e0*KE)), 0.0))
    #if !isapprox(abs(pitch), pitchabs, atol=1.e-1)
    #    warn("abs(pitch) != abspitch: ",pitchabs," ",pitch)
    #end
    return clamp(pitch,-1.0,1.0)
end

function get_pitch{S,T<:AbstractOrbitCoordinate}(M::AxisymmetricEquilibrium, c::T, r::Vector{S}, z::Vector{S})
    n = length(r)
    pitch = zeros(S,n)
    @inbounds for i=1:n
        pitch[i] = get_pitch(M, c, r[i], z[i])
    end
    return pitch
end

function get_pitch{T<:AbstractOrbitCoordinate}(M::AxisymmetricEquilibrium, c::T, path::OrbitPath)
    return get_pitch(M, c, path.r, path.z)
end

function get_pitch{T<:AbstractOrbitCoordinate,S<:Real}(M::AxisymmetricEquilibrium, c::T, r::S, z::S)
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

function get_kinetic_energy{T<:Real}(M::AxisymmetricEquilibrium, c::HamiltonianCoordinate, r::T, z::T)
    psi = M.psi_rz[r,z]
    return c.energy - 1e-3*M.phi[psi]
end

function get_kinetic_energy{S,T<:AbstractOrbitCoordinate}(M::AxisymmetricEquilibrium, c::T, r::Vector{S}, z::Vector{S})
    n = length(r)
    energy = zeros(S,n)
    @inbounds for i=1:n
        energy[i] = get_kinetic_energy(M, c, r[i], z[i])
    end
    return energy
end

function get_kinetic_energy{T<:AbstractOrbitCoordinate}(M::AxisymmetricEquilibrium, c::T, path::OrbitPath)
    return get_kinetic_energy(M, c, path.r, path.z)
end

function get_kinetic_energy{T<:AbstractOrbitCoordinate,S<:Real}(M::AxisymmetricEquilibrium, c::T, r::S, z::S)
    hc = HamiltonianCoordinate(M, c)
    return get_kinetic_energy(M, hc, r, z)
end

function get_kinetic_energy(M::AxisymmetricEquilibrium, o::Orbit)
    path = o.path
    c = o.coordinate
    return get_kinetic_energy(M, o.coordinate, o.path)
end

function classify(path::OrbitPath, pitch, axis; n=length(path))

    pp = take(zip(path.r,path.z),n)
    op = Limiter()
    for p in pp
        push!(op.vertices,p)
    end
    push!(op.vertices,first(pp))

    if in_vessel(op, axis)
        if all(sign.(pitch) .== sign(pitch[1]))
            if sign(pitch[1]) > 0.0
                class = Symbol("co_passing")
            else
                class = Symbol("ctr_passing")
            end
        else
            class = Symbol("potato")
        end
    else
        if all(sign.(pitch) .== sign(pitch[1]))
            class = Symbol("stagnation")
        else
            class = Symbol("trapped")
        end
    end

    return class
end

function hits_wall(path::OrbitPath, wall::Limiter)
    not_in_vessel = [~in_vessel(wall,p) for p in zip(path.r,path.z)]
    return any(not_in_vessel)
end

function hits_wall(o::Orbit, wall::Limiter)
    hits_wall(o.path, wall)
end
