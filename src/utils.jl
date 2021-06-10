function get_pitch(M::AbstractEquilibrium, c::HamiltonianCoordinate, r::T, z::T) where {T<:Number}
    psi = M(r,z)
    g = poloidal_current(M,psi)
    babs = norm(Bfield(M,r,z))
    phi = electric_potential(M,psi)
    KE = c.energy - 1e-3*phi
    f = -babs/(sqrt(2e3*e0*KE*c.m)*g*B0Ip_sign(M))
    pitch = f*(c.p_phi - c.q*e0*psi)
    #pitchabs = sqrt(max(1.0-(c.mu*babs/(1e3*e0*KE)), 0.0))
    #if !isapprox(abs(pitch), pitchabs, atol=1.e-1)
    #    warn("abs(pitch) != abspitch: ",pitchabs," ",pitch)
    #end
    return clamp(pitch,-1.0,1.0)
end

function get_pitch(M::AbstractEquilibrium, gcp::GCParticle, p_para::T, mu::T, r::T, z::T) where {T<:Number}
    Babs = norm(Bfield(M,r,z))
    m = gcp.m
    p = sqrt(p_para^2 + 2*m*Babs*mu)
    pitch = B0Ip_sign(M)*p_para/p
    return pitch
end

function get_pitch(M::AbstractEquilibrium, c::T, path::OrbitPath) where {T<:AbstractOrbitCoordinate}
    return get_pitch.(M, c, path.r, path.z)
end

function get_pitch(M::AbstractEquilibrium, c::T, r::S, z::S) where {T<:AbstractOrbitCoordinate, S<:Number}
    hc = HamiltonianCoordinate(M, c)
    return get_pitch(M, hc, r, z)
end

function get_pitch(M::AbstractEquilibrium, o::Orbit)
    return get_pitch.(M, o.coordinate, o.path)
end

function get_kinetic_energy(M::AbstractEquilibrium, c::HamiltonianCoordinate, r::T, z::T) where {T<:Number}
    psi = M(r,z)
    return c.energy - 1e-3*electric_potential(M,psi)
end

function get_kinetic_energy(M::AbstractEquilibrium, gcp::GCParticle, p_para::T, mu::T, r::T, z::T) where {T<:Number}
    Babs = norm(Bfield(M,r,z))
    m = gcp.m
    mc2 = m*c0^2
    p = sqrt(p_para^2 + 2*m*Babs*mu)
    KE = 1e-3*(hypot(p*c0, mc2) - mc2)/e0 #keV
end

function get_kinetic_energy(M::AbstractEquilibrium, c::T, path::OrbitPath) where {T<:AbstractOrbitCoordinate}
    return get_kinetic_energy.(M, c, path.r, path.z)
end

function get_kinetic_energy(M::AbstractEquilibrium, c::T, r::S, z::S) where {T<:AbstractOrbitCoordinate, S<:Number}
    hc = HamiltonianCoordinate(M, c)
    return get_kinetic_energy(M, hc, r, z)
end

function get_kinetic_energy(M::AbstractEquilibrium, o::Orbit)
    path = o.path
    c = o.coordinate
    return get_kinetic_energy(M, o.coordinate, o.path)
end

function classify(r, z, pitch, axis; n=length(r))

    op = Boundary(r,z)

    sign_p = sign.(pitch)
    if in_boundary(op, axis)
        if all(sign_p .== sign_p[1])
            if sign_p[1] > 0.0
                class = :co_passing
            else
                class = :ctr_passing
            end
        else
            if sign_p[argmax(r)] >= 0.0
                class = :potato
            else
                class = :ctr_passing #really ctr_potato
            end
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

function classify(path::OrbitPath, axis)
    return classify(path.r,path.z,path.pitch, axis, n=length(path))
end

"""

gcde_check(M, o; verbose=false)

For each point along the orbit path, compute the gyro-radius ρ⊥, the magnetic field
amplitude B and the square root of the maximum eigenvalue of the matrix ̂M. Use these quantities
to determine whether the usage of guiding-center drift equations (GCDE) could be expected to produce 
a significant discrepancy compared to full-orbit equations (FOE). The matrix ̂M is the matrix
found in equation (2) of D. Pfefferlé et al (2015) (https://doi.org/10.1088/0741-3335/57/5/054017)
and the condition to be evaluated is the condition found in equation (3) of the same paper.
In section 3 of the paper, the abritrary threshold of 0.073 is used. Here, we instead use 0.04.
If the condition is satisfied for a threshold of 0.04, no significant discrepancy between GCDE 
and FOE is to be expected. The equations to compute ̂M are from Appendix B of the same paper (B.2).

Return true if the criterion is fulfilled.
Return false if the criterion is not fulfilled.
    
Note! This function would have been a callback, had it not been for ForwardDiff.jacobian() being so slow.
"""
function gcde_check(M::AbstractEquilibrium, o::Orbit; verbose=false)
    m = o.coordinate.m # Mass of particle. kg
    KE = o.coordinate.energy # Kinetic energy. keV
    mc2 = m*c0*c0 # Rest energy. Joule
    KE_j = e0*KE*1e3 # Kinetic energy. Joule
    p_rel2 = ((KE_j + mc2)^2 - mc2^2)/(c0*c0) # Square of relativistic momentum
    q = o.coordinate.q*e0 # Particle charge. Coulomb

    for i=1:length(o.path.r)
        verbose && println("Evaluating gcde criterion for orbit path position $(i) of $(length(o.path.r))... ")
        R = o.path.r[i] # Major radius position of particle
        φ = o.path.phi[i] # Toroidal angle position of particle
        Z = o.path.z[i] # Vertical position of particle
        pitch = o.path.pitch[i] # Pitch of particle

        p_perp2 = p_rel2*(1-pitch^2) # Square of relativistic perpendicular momentum
        B = Bfield(M, R, Z) # Magnetic field vector at particle position.
        Babs = norm(B) # Magnetic field magnitude. Tesla

        r_g = sqrt(p_perp2) / (abs(q)*Babs) # Gyroradius at particle position. Meter

        jacB = ForwardDiff.jacobian(x->Bfield(M,x[1],x[3]),SVector{3}(R,zero(R),Z)) # Jacobian of magnetic field
        #verbose && println("Jacobian(B): $(jacB)")
        chrB = [0.0 -R*B[2] 0.0;B[2]/R B[1]/R 0.0;0.0 0.0 0.0] # The Christoffel symbol component of the covariant derivative
        Bij = jacB .+ chrB

        gij = [1.0 0.0 0.0;0.0 R*R 0.0;0.0 0.0 1.0]

        Pjk = inv(gij) .- ((1/(Babs^2)) .*[B[1]*B[1] B[1]*B[2] B[1]*B[3];B[2]*B[1] B[2]*B[2] B[2]*B[3];B[3]*B[1] B[3]*B[2] B[3]*B[3]])

        dxmdrk = [cos(φ) -R*sin(φ) 0.0;sin(φ) R*cos(φ) 0.0;0.0 0.0 1.0]

        Vim = Bij*Pjk*dxmdrk
        Mhat = Vim*inv(gij)*transpose(Vim)
        #Mhat = Vim*gij*transpose(Vim) # Exact formulation in D. Pfefferlé et al (2015). Should be inv(gij)?

        λmax = maximum(eigvals(Mhat)) # Could use λmax = tr(Mhat) instead to speed up (trace of matrix), but less accurate.

        criterion = sqrt(λmax)*r_g / Babs

        if criterion > 0.04 # criterion violated
            verbose && println("criterion violated!")
            verbose && println("criterion: $(criterion)")
            verbose && println("r_g: $(r_g) m")
            verbose && println("|B|: $(Babs) T")
            verbose && println("√(λmax): $(sqrt(λmax))")
            return false
        end
    end

    return true
end

function hits_wall_path(path::OrbitPath, wall::Wall)

    not_in_vessel = [~in_vessel(wall,p) for p in zip(path.r,path.z)]
    ind = findfirst(not_in_vessel)
    if ind != nothing
        return false, path
    else
        ind = ind-1
        return true, OrbitPath(path.vacuum, path.drift,
                               path.energy[1:ind], path.pitch[1:ind],
                               path.r[1:ind], path.z[1:ind], path.phi[1:ind],
                               path.dt[1:ind])
    end
end

function hits_wall(path::OrbitPath, wall::Wall)
    not_in_vessel = [~in_vessel(wall,p) for p in zip(path.r,path.z)]
    return any(not_in_vessel)
end

function lorentz_factor(K_keV,m)
    KE_j = e0*1e3*K_keV
    mc2 = m*c0^2
    p_rel2 = ((KE_j + mc2)^2 - mc2^2)/(c0*c0)
    return sqrt(1 + p_rel2/((m*c0)^2))
end

function lorentz_factor(p::GCParticle)
    return lorentz_factor(p.energy, p.m)
end

function lorentz_factor(p::Particle)
    v = hypot(p.vr,p.vt,p.vz)
    beta = v/c0
    return inv(sqrt(1 - beta^2))
end

function cyclotron_frequency(M,p::T) where T <: AbstractParticle
    gamma = lorentz_factor(p)
    Babs = norm(Bfield(M, p.r, p.z))
    return abs(p.q*e0)*Babs/(gamma*p.m)
end

function cyclotron_period(M,p::T) where T <: AbstractParticle
    return 2*pi/cyclotron_frequency(M,p)
end

function gyroradius(M,p::GCParticle)
    γ = lorentz_factor(p)
    Babs = norm(Bfield(M, p.r, p.z)) # Magnetic field magnitude. Tesla
    m = p.m # Mass of particle
    KE = p.energy # Kinetic energy. keV
    mc2 = m*c0*c0 # Rest energy. Joule
    KE_j = e0*KE*1e3 # Kinetic energy. Joule
    p_rel2 = ((KE_j + mc2)^2 - mc2^2)/(c0*c0) # Relativistic momentum
    p_perp2 = p_rel2*(1-p.pitch^2) # Square of relativistic perpendicular momentum

    return sqrt(p_perp2) / (abs(p.q)*B)
end

function gyroradius(M,p::Particle)
    B = Bfield(M, p.r, p.z)
    vpara, vperp = parallel_perpendicular([p.vr p.vt p.vz], B)
    γ = lorentz_factor(p)
    return γ*(p.m)*vperp / (abs(p.q)*norm(B))
end

function normalize(M::AbstractEquilibrium, hc::HamiltonianCoordinate)
    E = hc.energy
    B0 = norm(Bfield(M, magnetic_axis(M)...))
    mu = (abs(B0)/(E*e0*1e3))*hc.mu
    flux = abs(-(psi_limits(M)...))
    Pphi = (B0Ip_sign(M)/(e0*flux))*hc.p_phi

    return E, Pphi, mu
end

function normalize(M::AbstractEquilibrium, KE, p, r, z; amu=H2_amu, q=1)
    normalized(M,HamiltonianCoordinate(M,KE,p,r,z,amu*mass_u,q))
end

"""
    v_para, v_perp = parallel_perpendicular(v, B)

Returns parallel velocity component (`v_para`) w.r.t to B and the perpendicular vector (`v_perp`)
"""
function parallel_perpendicular(v, B)
    Bhat = B/norm(B)
    v_para = dot(v,Bhat)
    v_perp = cross(-Bhat,cross(Bhat,v))
    return v_para, v_perp
end

"""
    a, c = perpendicular_vectors(b)

Calculates normalized vectors that are perpendicular to `b` such that `a x c = b_norm`
"""
function perpendicular_vectors(B::T) where T
    Bhat = B/norm(B)

    if abs(Bhat[3]) == 1
        a = convert(T, [1, 0, 0])
        c = convert(T, [0, 1, 0])
    else
        if Bhat[3] == 0
            a = convert(T, [0, 0, 1])
            c = convert(T, [Bhat[2], -Bhat[1], 0])
        else
            a = convert(T, [Bhat[2], -Bhat[1], 0]/norm(Bhat[1:2]))
            c = convert(T, -[a[2] , -a[1] , (a[1]*Bhat[2] - a[2]*Bhat[1])/Bhat[3]])
            c = c/norm(c)
            if Bhat[3] < 0
                c = -c
            end
        end
    end
    return a, c
end