struct FullOrbitPath{T<:Real}
    dt::Union{T,Vector{T}}
    r::Vector{T}
    phi::Vector{T}
    z::Vector{T}
    vr::Vector{T}
    vphi::Vector{T}
    vz::Vector{T}
end

function FullOrbitPath(T::DataType=Float64)
    FullOrbitPath(zero(T),T[],T[],T[],T[],T[],T[])
end

Base.length(op::FullOrbitPath) = length(op.r)

function borispush(M, v, u, dt)
        x = u[1]
        y = u[2]
        r = sqrt(x*x + y*y)
        phi = atan(y,x)
        z = u[3]

        q_m_half_dt = (0.5*dt*e0/(mass_u*H2_amu))

        F = fields(M,r,z)
        cp = cos(phi)
        sp = sin(phi)
        B = SVector{3}(F.B[1]*cp - F.B[2]*sp,F.B[1]*sp + F.B[2]*cp,F.B[3])
        E = SVector{3}(F.E[1]*cp - F.E[2]*sp,F.E[1]*sp + F.E[2]*cp,F.E[3])

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

function get_full_orbit(M, r, phi, z, vr, vphi, vz; dt= 1e-9, tmax = 100*1e-6)

    tau_c = 2*pi*(mass_u*H2_amu)/(M.b(r,z)*e0)

    t = 0:dt:tmax
    nstep = length(t)
    x_arr = zeros(nstep)
    y_arr = zeros(nstep)
    z_arr = zeros(nstep)
    vx_arr = zeros(nstep+1)
    vy_arr = zeros(nstep+1)
    vz_arr = zeros(nstep+1)

    sp = sin(phi)
    cp = cos(phi)
    u = SVector{3}([r*cp,r*sp,z])
    v = SVector{3}([vr*cp - vphi*sp,vr*sp + vphi*cp,vz])
    x_arr[1] = u[1]
    y_arr[1] = u[2]
    z_arr[1] = u[3]
    v = borispush(M,v,u,-0.5*dt)
    vx_arr[1] = v[1]
    vy_arr[1] = v[2]
    vz_arr[1] = v[3]
    for i=2:nstep
        v = borispush(M,v,u,dt)
        x_arr[i] = x_arr[i-1] + v[1]*dt
        y_arr[i] = y_arr[i-1] + v[2]*dt
        z_arr[i] = z_arr[i-1] + v[3]*dt
        u  = SVector{3}(x_arr[i],y_arr[i],z_arr[i])
        vx_arr[i] = v[1]
        vy_arr[i] = v[2]
        vz_arr[i] = v[3]
    end
    v = borispush(M,v,u,0.5*dt)
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

    return FullOrbitPath(dt, r_arr,phi_arr,z_arr,vr_arr,vphi_arr,vz_arr)
end

function hits_wall(M, r, phi, z, vr, vphi, vz, wall; dt=1e-9, tmax=100e-6)
    t = 0:dt:tmax
    nstep = length(t)

    sp = sin(phi)
    cp = cos(phi)
    u = SVector{3}([r*cp,r*sp,z])
    v = SVector{3}([vr*cp - vphi*sp,vr*sp + vphi*cp,vz])
    hit = false
    v = borispush(M,v,u,-0.5*dt)
    T = 0.0
    for i=1:nstep
        v = borispush(M,v,u,dt)
        u = u .+ v*dt
        rr = sqrt(u[1]^2 + u[2]^2)
        zz = u[3]
        if ~in_vessel(wall,(rr,zz))
            hit = true
            break
        end
        T += dt
    end

    return hit, T
end
