immutable OrbitPath{T}
    r::Vector{T}
    z::Vector{T}
    phi::Vector{T}
    dt::Vector{T}
end

function OrbitPath(T::DataType=Float64)
    OrbitPath(T[],T[],T[],T[])
end

Base.length(op::OrbitPath) = length(op.r)

immutable Orbit{T,S<:AbstractOrbitCoordinate{Float64}}
    coordinate::S
    class::Symbol
    omega_p::T
    omega_t::T
    path::OrbitPath{T}
end

function Orbit{T}(c::AbstractOrbitCoordinate{T})
    return Orbit(c, :incomplete, zero(T), zero(T), OrbitPath(T))
end

function Orbit{T}(c::AbstractOrbitCoordinate{T}, op::OrbitPath{T})
    return Orbit(c, :incomplete, zero(T), zero(T), op)
end

function make_gc_ode{T<:AbstractOrbitCoordinate}(M::AxisymmetricEquilibrium, c::T)
    oc = HamiltonianCoordinate(M, c)
    res = DiffBase.GradientResult(rand(2))
    function vgc(t, y, ydot)
        ri = [y[1],y[3]]
        psi = M.psi(ri)
        g = M.g(psi)
        B = M.B(ri)
        ForwardDiff.gradient!(res, M.b, ri)
        gradB = DiffBase.gradient(res)
        babs = DiffBase.value(res)

        Wperp = oc.mu*babs
        Wpara = 1e3*e0*oc.energy - Wperp
        vpara = -babs*(oc.p_phi - oc.q*e0*psi)/(oc.amu*mass_u*g)

        vd = (1/(oc.q*e0))*(Wperp + 2*Wpara)*cross(B,[gradB[1],0.0,gradB[2]])/(babs^3)
        vg = vpara*B/babs + vd
        ydot[1] = vg[1]
        ydot[2] = vg[2]/ri[1]
        ydot[3] = vg[3]
        return Sundials.CV_SUCCESS
    end
end

function get_orbit(M::AxisymmetricEquilibrium, E, pitch_i, ri, zi, amu, q::Int, nstep::Int, tmax)

    hc = HamiltonianCoordinate(M, E, pitch_i, ri, zi, amu=amu, q=q)

    #Function to check if in function
    hits_boundary = [false]
    in_boundary = function ib(x)
        if (M.r_domain[1] < x[1] < M.r_domain[2]) && (M.z_domain[1] < x[3] < M.z_domain[2])
            return true
        else
            hits_boundary[1] = true
            return false
        end
    end

    #Function to check if orbit is complete
    r0 = [ri, 0.0, zi]
    r1 = copy(r0)
    initial_dir = [0]
    npol = [0]
    orbit_complete = [false]
    is_complete = function ic(r2)
        if initial_dir[1] == 0
            initial_dir[1] = sign(r2[3]-r0[3])
            r1[:] = convert(Vector, r2)
            return false
        end

        dr10 = r1 - r0
        dr20 = r2 - r0
        dr21 = r2 - r1

        if sign(dr10[3]*dr20[3]) < 0
            npol[1] = npol[1] + 1
            if norm(dr10[1:2:3]) < norm(dr21[1:2:3]) && sign(dr20[3]) == initial_dir[1]
                orbit_complete[1] = npol[1] == 2 || npol[1] == 4
                return true
            end
        end
        r1[:] = convert(Vector, r2)
        return false
    end

    #Function to call between timesteps
    cb = function callback(mem, t, x)
        !is_complete(x) && in_boundary(x)
    end

    if !in_boundary(r0)
        error("Starting point outside boundary")
    end

    f = make_gc_ode(M, hc)
    t = 1e-6*collect(linspace(0.0,tmax,nstep))

    res = Sundials.cvode(f, r0, t, reltol=1e-8,abstol=1e-12, callback = cb)

    r = res[:,1]
    phi = res[:,2]
    z = res[:,3]

    n = length(r)
    dt = fill(t[2]-t[1],n)

    #correct for last timestep
    ydot = zeros(3)
    flag = f(0.0, r0, ydot)
    dt[end] = sqrt((r[end] - r[1])^2 + (z[end] - z[1])^2)/norm(ydot[1:2:3])

    path = OrbitPath(r,z,phi,dt)
    pitch = get_pitch(M, hc, path)
    rmax, ind = findmax(r)
    pitch_rmax = pitch[ind]

    c = EPRCoordinate(E, pitch_rmax, rmax, z[ind], hc.amu, hc.q)

    if !orbit_complete[1] || hits_boundary[1]
        #warn("Incomplete orbit")
        return Orbit(c, path)
    end

    T_p = sum(dt)
    omega_p = 2pi/T_p
    omega_t = abs((phi[end]-phi[1]))/T_p

    class = classify(path, pitch, M.axis)

    return Orbit(c, class, omega_p, omega_t, path)
end

function get_orbit(M::AxisymmetricEquilibrium, E, p, r, z; amu=H2_amu, q=1, nstep=3000, tmax=500.0)
    return get_orbit(M, E, p, r, z, amu, q, nstep, tmax)
end

function get_orbit(M::AxisymmetricEquilibrium, c::EPRCoordinate; nstep=3000, tmax=500.0)
    o = get_orbit(M, c.energy, c.pitch, c.r, c.z, c.amu, c.q, nstep, tmax)
    if o.class != :incomplete
        maximum(o.path.r) > c.r && return Orbit(c)
    end
    return o
end

function Base.show(io::IO, orbit::Orbit)
    class_str = Dict(:trapped=>"Trapped ",:co_passing=>"Co-passing ",:ctr_passing=>"Counter-passing ",
                     :stagnation=>"Stagnation ",:potato=>"Potato ",:incomplete=>"Incomplete ")
    println(io, typeof(orbit))
    println(io, class_str[orbit.class],"Orbit:")
    @printf(io, " ωₚ = %.3f rad/ms\n", orbit.omega_p*1e-3)
    @printf(io, " ωₜ = %.3f rad/ms\n", orbit.omega_t*1e-3)
end

function plot_orbit(o::Orbit;rlim=(1.0,2.5),zlim=(-1.25,1.25),xlim = (-2.3,2.3),ylim=(-2.3,2.3))
#    pyplot()
    t = cumsum(o.path.dt)*1e6
    n = length(t)
    r = o.path.r
    z = o.path.z
    phi = o.path.phi
    x = r.*cos.(phi)
    y = r.*sin.(phi)
    l = @layout [a{0.38w} b]
    p1 = plot(layout=l)

    rr = vec(hcat(r[1:end-1],r[2:end],fill(NaN,n-1))')
    zz = vec(hcat(z[1:end-1],z[2:end],fill(NaN,n-1))')
    plot!(p1,rr,zz,line_z=t,c=:bluesreds,label=o.class,
          colorbar=false,subplot=1,xlim=rlim,
          ylim=zlim,xlabel="R [m]",ylabel="Z [m]",aspect_ratio=1.0)

    xx = vec(hcat(x[1:end-1],x[2:end],fill(NaN,n-1))')
    yy = vec(hcat(y[1:end-1],y[2:end],fill(NaN,n-1))')
    plot!(p1,xx,yy,line_z=t,c=:bluesreds,label=o.class,subplot=2,
          xlim=xlim,ylim=ylim,xlabel="X [m]",
          ylabel="Y [m]",aspect_ratio=1.0)

    return p1
end
