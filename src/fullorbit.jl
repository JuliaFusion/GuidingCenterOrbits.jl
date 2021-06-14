struct FullOrbitPath{T<:Number}
    dt::Union{T,Vector{T}}
    r::Vector{T}
    phi::Vector{T}
    z::Vector{T}
    vr::Vector{T}
    vphi::Vector{T}
    vz::Vector{T}
end

mutable struct FullOrbitStatus
    errcode::Int
end
FullOrbitStatus() = FullOrbitStatus(1)

function FullOrbitPath(T::DataType=Float64)
    FullOrbitPath(zero(T),T[],T[],T[],T[],T[],T[])
end

Base.length(op::FullOrbitPath) = length(op.r)

function borispush(M, m, v, u, dt)
        x = u[1]
        y = u[2]
        z = u[3]

        q_m_half_dt = (0.5*dt*e0/m)

        B, E = fields(M,x,y,z)

        t = q_m_half_dt*B

        half_acc = q_m_half_dt*E
        v_minus = v .+ half_acc
        v_minus_x_t = cross(v_minus,t)
        v_prime = v_minus + v_minus_x_t

        s = (2.0/(1+dot(t,t)))*t
        v_prime_x_s = cross(v_prime,s)
        v_plus = v_minus + v_prime_x_s

        v = v_plus .+ half_acc

        return v
end

"""
    integrate(M, pc::Particle; dt=0.01, tmax=1000) -> FullOrbitPath, stat

Integrate the particle up to tmax normalized to the cyclotron_period with time step dt is a fraction of a cyclotron period. Default tmax is 1000
"""
function integrate(M::AbstractEquilibrium, pc::Particle; dt= 0.01, tmax = 1000)

    if lorentz_factor(pc) > 1.05
        @warn "Relativistic Full Orbit has not been implemented: Lorentz factor > 1.05"
    end

    T_c = cyclotron_period(M,pc)
    dt_sec = dt*T_c
    tmax_sec = tmax*T_c
    t = 0:dt:tmax

    nstep = length(t)
    x_arr = zeros(nstep)
    y_arr = zeros(nstep)
    z_arr = zeros(nstep)
    vx_arr = zeros(nstep+1)
    vy_arr = zeros(nstep+1)
    vz_arr = zeros(nstep+1)

    sp, cp = sincos(pc.phi)
    u = SVector{3}([pc.r*cp,pc.r*sp,pc.z])
    v = SVector{3}([pc.vr*cp - pc.vt*sp,pc.vr*sp + pc.vt*cp,pc.vz])
    x_arr[1] = u[1]
    y_arr[1] = u[2]
    z_arr[1] = u[3]
    v = borispush(M,pc.m,v,u,-0.5*dt_sec)
    vx_arr[1] = v[1]
    vy_arr[1] = v[2]
    vz_arr[1] = v[3]
    for i=2:nstep
        v = borispush(M,pc.m,v,u,dt_sec)
        x_arr[i] = x_arr[i-1] + v[1]*dt_sec
        y_arr[i] = y_arr[i-1] + v[2]*dt_sec
        z_arr[i] = z_arr[i-1] + v[3]*dt_sec
        u  = SVector{3}(x_arr[i],y_arr[i],z_arr[i])
        vx_arr[i] = v[1]
        vy_arr[i] = v[2]
        vz_arr[i] = v[3]
    end
    v = borispush(M,pc.m,v,u,0.5*dt_sec)
    vx_arr[end] = v[1]
    vy_arr[end] = v[2]
    vz_arr[end] = v[3]

    r_arr = sqrt.(x_arr.^2 .+ y_arr.^2)
    phi_arr = atan.(y_arr,x_arr)
    sp = sin.(phi_arr)
    cp = cos.(phi_arr)
    #velocity is 0.5 step behind position
    vr_arr = [0.5*(vx_arr[i]+vx_arr[i+1])*cp[i] .+ 0.5*(vy_arr[i] + vy_arr[i+1])*sp[i] for i=1:nstep]
    vphi_arr = [-0.5*(vx_arr[i]+vx_arr[i+1])*sp[i] .+ 0.5*(vy_arr[i] + vy_arr[i+1])*cp[i] for i=1:nstep]
    vz_arr = [0.5(vz_arr[i] + vz_arr[i+1]) for i=1:nstep]

    return FullOrbitPath(dt_sec, r_arr, phi_arr,z_arr,vr_arr,vphi_arr,vz_arr), FullOrbitStatus(0)
end

"""
    get_full_orbit(M, pc::Particle; dt=T_c/100, tmax=1000*T_c) -> FullOrbitPath, stat

Integrate the full orbit up to tmax (μs) with time step dt (μs). Default tmax is 1000*cyclotron_period.
"""
function get_full_orbit(M::AbstractEquilibrium, pc::Particle; kwargs...)
    return integrate(M, pc; kwargs...)
end

"""
    get_full_orbit(M, gcp::GCParticle; dt=T_c/100, tmax=1000*T_c) -> FullOrbitPath, stat

Integrate the full orbit up to tmax (μs) with time step dt (μs). Default tmax is 1000*cyclotron_period.
"""
function get_full_orbit(M::AbstractEquilibrium, gcp::GCParticle; gamma = 0.0, verbose=false, kwargs...)

    # Turn GCParticle into a Particle
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
    b = B/norm(B)
    a, c = perpendicular_vectors(B)

    pvec = p_perp0*cos(gamma)*a .+ p_perp0*sin(gamma)*c .+ p_para0*b

    v = pvec/(gcp.m*lorentz_factor(gcp))

    Ω_c = cyclotron_frequency(M,gcp)
    r_gyro = cross(v,b)/Ω_c

    verbose && println("|r_gyro|: $(norm(r_gyro))")
    
    r_p = SVector{3}(gcp.r, zero(gcp.r), gcp.z) .- r_gyro

    r = norm(r_p[1:2])
    phi = atan(r_p[2],r_p[1])
    z = r_p[3]
    sp,cp = sincos(phi)
    vr = v[1]*cp + v[2]*sp
    vphi = -v[1]*sp + v[2]*cp
    vz = v[3]

    pc = Particle(r,phi,z,vr,vphi,vz,gcp.m,gcp.q)

    return integrate(M, pc::Particle; kwargs...)
end


function hits_wall(M, pc::Particle, wall; dt=1e-3, tmax=100)

    dt_sec = dt*1e-6
    t = 0:dt:tmax
    nstep = length(t)

    sp,cp = sincos(pc.phi)
    u = SVector{3}([pc.r*cp,pc.r*sp,z])
    v = SVector{3}([pc.vr*cp - pc.vt*sp, pc.vr*sp + pc.vt*cp,pc.vz])
    hit = false
    v = borispush(M,pc.m,v,u,-0.5*dt_sec)
    T = 0.0
    for i=1:nstep
        v = borispush(M,pc.m,v,u,dt_sec)
        u = u .+ v*dt_sec
        rr = sqrt(u[1]^2 + u[2]^2)
        zz = u[3]
        if ~in_vessel(wall,(rr,zz))
            hit = true
            break
        end
        T += dt_sec
    end

    return hit, T
end
