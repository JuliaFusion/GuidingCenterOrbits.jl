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

For each point along the orbit path, compute the gyro-radius r_g, the magnetic field
amplitude B and the square root of the maximum eigenvalue of the matrix ̂M. Use these quantities
to determine whether the usage of guiding-center drift equations (GCDE) could be expected to produce 
a significant discrepancy compared to full-orbit equations (FOE). The matrix ̂M is the matrix
found in equation (2) of D. Pfefferlé et al (2015) (https://doi.org/10.1088/0741-3335/57/5/054017)
and the condition to be evaluated is the condition found in equation (3) of the same paper.
In section 3 of the paper, the abritrary threshold of 0.073 is used. The default is the same here.
If the condition is satisfied for the threshold, no significant discrepancy between GCDE 
and FOE is to be expected. The equations to compute ̂M are from Appendix B of the same paper (B.2).

Return true if the criterion is fulfilled, and it is ok to use GCDE.
Return false if the criterion is not fulfilled, and FOE should be used instead.
"""
function gcde_check(M::AbstractEquilibrium, o::Orbit; threshold=0.073, verbose=false)

    g_rz = (x) -> poloidal_current(M, M(x[1],x[2])) # Poloidal current function (F=R*Bt) as a function of R,Z
    m = o.coordinate.m # Mass of particle. kg
    KE = o.coordinate.energy # Kinetic energy. keV
    mc2 = m*c0*c0 # Rest energy. Joule
    KE_j = e0*KE*1e3 # Kinetic energy. Joule
    p_rel2 = ((KE_j + mc2)^2 - mc2^2)/(c0*c0) # Square of relativistic momentum
    q = o.coordinate.q*e0 # Particle charge. Coulomb

    maxcrit = 0.0
    for i=1:length(o.path.r)
        verbose && println("Evaluating gcde criterion for orbit path position $(i) of $(length(o.path.r))... ")
        R = o.path.r[i] # Major radius position of particle
        φ = o.path.phi[i] # Toroidal angle position of particle
        Z = o.path.z[i] # Vertical position of particle
        pitch = o.path.pitch[i] # Pitch of particle

        p_perp2 = p_rel2*(1-pitch^2) # Square of relativistic perpendicular momentum
        B = Equilibrium.Bfield(M, R, Z) # Magnetic field vector at particle position.
        Babs = norm(B) # Magnetic field magnitude. Tesla

        r_g = sqrt(p_perp2) / (abs(q)*Babs) # Gyroradius at particle position. Meter

        psi = M(R,Z) # Poloidal flux in Weber/rad
        F = poloidal_current(M,psi) # Poloidal current function in meter*Tesla
        grad_psi = SVector{2,Float64}(psi_gradient(M,R,Z)) # [dψ/dR,dψ/dZ]. SVector for efficiency
        cc = cocos(M) # Cocos factor
        cocos_factor = cc.sigma_RpZ*cc.sigma_Bp/((2pi)^cc.exp_Bp)
        
        grad_F = SVector{2,Float64}(ForwardDiff.gradient(g_rz,SVector{2}(R,Z))) # [∂F/∂r,∂F/∂z] SVector for efficiency
        J_grad_psi = SMatrix{2,2}(ForwardDiff.jacobian(x->psi_gradient(M,x[1],x[2]),SVector{2}(R,Z))) # [d^2ψ/dr^2 d^2ψ/dzdr;d^2ψ/drdz d^2ψ/dz^2] SMatrix for efficiency

        D = SMatrix{3,3}(cocos_factor*inv(R)*(J_grad_psi[2,1]-grad_psi[2]*inv(R)),grad_F[1]-F*inv(R),cocos_factor*inv(R)*(-J_grad_psi[1,1]+inv(R)*grad_psi[1]),-F*inv(R),cocos_factor*grad_psi[2],0.0,cocos_factor*inv(R)*J_grad_psi[2,2],grad_F[2],-cocos_factor*inv(R)*J_grad_psi[1,2]) # Matrix D as in equation (2) in D. Pfefferlé et al (2015). Note D_{ij}=B_{i;j} in appendix B.

        G = SMatrix{3,3}([1.0 0.0 0.0;0.0 inv(R*R) 0.0;0.0 0.0 1.0])

        P = SMatrix{3,3}(G .- ((1/(Babs^2)) .*[grad_psi[2]*grad_psi[2]*inv(R^2) cocos_factor*grad_psi[2]*F*inv(R^3) -grad_psi[1]*grad_psi[2]*inv(R^2);cocos_factor*grad_psi[2]*F*inv(R^3) (F^2)*inv(R^4) -cocos_factor*F*grad_psi[1]*inv(R^3);-grad_psi[1]*grad_psi[2]*inv(R^2) -cocos_factor*F*grad_psi[1]*inv(R^3) grad_psi[1]*grad_psi[1]*inv(R^2)]))

        Λ  = SMatrix{3,3}([cos(φ) -R*sin(φ) 0.0;sin(φ) R*cos(φ) 0.0;0.0 0.0 1.0])

        Mhat = transpose(Λ)*transpose(P)*transpose(D)*G*D*P*Λ

        λmax = maximum(eigvals(Array(Mhat))) # Could use λmax = tr(Mhat) instead to speed up. But less accurate.

        criterion = sqrt(λmax)*r_g / Babs

        if criterion > maxcrit
            maxcrit = criterion
        end

        if criterion > threshold # criterion violated
            verbose && println("--- Criterion violated! ---")
            verbose && println("Criterion > Threshold: $(criterion) > $(threshold)")
            verbose && println("- Gyroradius (r_g) at bad position: $(r_g) m")
            verbose && println("- |B| at bad position: $(Babs) T")
            verbose && println("- √(λmax) at bad position: $(sqrt(λmax))")
            verbose && println("√(λmax) * r_g / |B| > Threshold")
            return false
        end
    end
    verbose && println("- Criterion was not violated at any point along the orbit path.")
    verbose && println("maximum(criterion) < threshold: $(maxcrit) < $(threshold)")
    verbose && println("- gcde ok to use!")
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

"""
    velocity(gcp::GCParticle, gamma)

Returns the velocity of the particle in cartesian coordinates. Assumes toroidal angle = 0
"""
function velocity(M::AbstractEquilibrium, gcp::GCParticle, gamma)
    mc2 = gcp.m*c0^2

    p0 = sqrt(((1e3*e0*gcp.energy + mc2)^2 - mc2^2)/(c0*c0)) # The initial particle momentum
    if abs(gcp.pitch) == 1.0
        pitch0 = sign(gcp.pitch)*prevfloat(abs(gcp.pitch))
    else
        pitch0 = gcp.pitch
    end

    p_para0 = p0*pitch0*B0Ip_sign(M) # The initial parallel momentum
    p_perp0 = p0*sqrt(1.0 - pitch0^2) # The initial perpendicular momentum

    B = Bfield(M,gcp.r,0.0,gcp.z)
    Babs = norm(B)
    b = B/Babs
    a, c = perpendicular_vectors(B)

    pvec = p_perp0*cos(gamma)*a .+ p_perp0*sin(gamma)*c .+ p_para0*b

    v = pvec/(gcp.m*lorentz_factor(gcp))
    return v
end

"""
    gyro_step(M::AbstractEquilibrium, gcp, gamma) -> r_gyro

Returns gyro-step such that the r_p = r_gc - r_gyro in cartesian coordinates. Assumes toroidal angle=0
"""
function gyro_step(M::AbstractEquilibrium, gcp::GCParticle, gamma)

    v = velocity(M, gcp, gamma)

    B = Bfield(M,gcp.r,0.0,gcp.z)
    Babs = norm(B)
    b = B/Babs

    # First order
    Ω_c = cyclotron_frequency(M,gcp)
    r_gyro = cross(v,b)/Ω_c

    # Second order
    # Belova, E. V., N. N. Gorelenkov, and C. Z. Cheng. "Self-consistent equilibrium model of low aspect-
    # ratio toroidal plasma with energetic beam ions." Physics of Plasmas (1994-present) 10.8 (2003):
    vpara = dot(v,b)
    cB = curlB(M,gcp.r,0.0,gcp.z)/Babs
    term1 = vpara*dot(b,cB)/Ω_c
    gB = gradB(M,gcp.r,0.0,gcp.z)
    term2 = -dot(r_gyro,gB)/(2*Babs)

    correction = (1 - term1 - term2)
    if (correction <= 0.0) || (correction >= 2.0)
        @warn "Gyro correction results in negative distances or too large of a shift: $(correction)"
    end

    r_gyro = r_gyro*correction

    return r_gyro
end

"""
    larmor_radius(M::AbstractEquilibrium, gcp::GCParticle)

Returns the average larmor radius
"""
function larmor_radius(M::AbstractEquilibrium, gcp::GCParticle; N=10)
    return mean(norm(gyro_step(M, gcp, g)) for g in range(0,2pi,length=N))
end