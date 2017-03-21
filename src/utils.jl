function get_pitch{T<:Real}(M::AxisymmetricEquilibrium, c::HamiltonianCoordinate, r::T, z::T)
    rz = [r,z]
    psi = M.psi(rz)
    g = M.g([psi])
    babs = M.b(rz)

    f = -babs/(sqrt(2e3*e0*c.energy*mass_u*c.amu)*g*M.sigma)
    pitch = f*(c.p_phi - c.q*e0*psi)
    pitchabs = sqrt(max(1.0-(c.mu*babs/(1e3*e0*c.energy)), 0.0))
    if !isapprox(abs(pitch), pitchabs, atol=1.e-1)
        warn("abs(pitch) != abspitch: ",pitchabs," ",pitch)
    end
    return clamp(pitch,-1.0,1.0)
end

function get_pitch{T<:AbstractOrbitCoordinate}(M::AxisymmetricEquilibrium, c::T, path::OrbitPath)
    n = length(path.r)
    pitch = zeros(eltype(path.r),n)
    @inbounds for i=1:n
        pitch[i] = get_pitch(M, c, path.r[i], path.z[i])
    end
    return pitch
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

function classify(path::OrbitPath, pitch, axis)

    pp = zip(path.r,path.z)
    op = Polygon()
    for p in pp
        push!(op.vertices,p)
    end
    push!(op.vertices,first(pp))

    if in_polygon(axis, op)
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

function hits_wall(o::Orbit, wall::Polygon)
    not_in_poly = [~in_polygon(p, wall) for p in zip(o.path.r,o.path.z)]
    return any(not_in_poly)
end
