struct OrbitPath{T}
    r::Vector{T}
    z::Vector{T}
    phi::Vector{T}
    pitch::Vector{T}
    energy::Vector{T}
    dt::Vector{T}
#    dl::Vector{T}
end

function OrbitPath(T::DataType=Float64)
    OrbitPath(T[],T[],T[],T[],T[],T[])
end

Base.length(op::OrbitPath) = length(op.r)

struct Orbit{T,S<:AbstractOrbitCoordinate{Float64}}
    coordinate::S
    class::Symbol
    tau_p::T
    tau_t::T
    path::OrbitPath{T}
end

function Orbit(c::AbstractOrbitCoordinate{T},class=:incomplete) where {T}
    return Orbit(c, class, zero(T), zero(T), OrbitPath(T))
end

mutable struct GCStatus{T<:Number}
    errcode::Int
    ri::SArray{Tuple{3},T,1,3}
    vi::SArray{Tuple{3},T,1,3}
    nr::Int
    rm::T
    zm::T
    pm::T
    tm::T
    tau_p::T
    tau_t::T
    poloidal_complete::Bool
    hits_boundary::Bool
    class::Symbol
end

function GCStatus(T=Float64)
    z = zero(T)
    return GCStatus(1,SVector{3}(z,z,z),SVector{3}(z,z,z),0,z,z,z,z,z,z,false,false,:incomplete)
end

function make_gc_ode(M::AxisymmetricEquilibrium, c::T, stat::GCStatus) where {T<:AbstractOrbitCoordinate}
    oc = HamiltonianCoordinate(M, c)
    ode = function f(y,p::Bool,t)
        stat
        r = y[1]
        z = y[3]

        F = fields(M,r,z)
        psi = F.psi
        g = F.g
        B = F.B
        E = F.E

        babs = norm(B)
        gradB = Interpolations.gradient(M.b,r,z)
        gradB = SVector{3}(gradB[1],zero(gradB[1]),gradB[2])
        Wperp = oc.mu*babs
        vpara = -babs*(oc.p_phi - oc.q*e0*psi)/(oc.m*g)
        Wpara = 0.5*oc.m*vpara^2

        # ExB Drift
        v_exb = cross(E, B)/babs^2

        # GradB Drift
        b1 = cross(B,gradB)/(babs^3)
        v_grad = Wperp*b1/(oc.q*e0)

        # Curvature Drift
        v_curv = 2*Wpara*b1/(oc.q*e0)

        # Guiding Center Velocity
        v_gc = (vpara*B/babs + v_exb + v_grad + v_curv)
        return SVector{3}(v_gc[1], v_gc[2]/r, v_gc[3])
    end
    return ode
end

function integrate(M::AxisymmetricEquilibrium, gcp::GCParticle,
                   dt, tmax, integrator, wall::Union{Nothing,Limiter}, interp_dt, classify_orbit::Bool,
                   one_transit::Bool, store_path::Bool, maxiter::Int, adaptive::Bool, autodiff::Bool,verbose::Bool)

    stat = GCStatus(typeof(gcp.r))

    r0 = @SVector [gcp.r,zero(typeof(gcp.r)),gcp.z]
    if !((M.r[1] < gcp.r < M.r[end]) && (M.z[1] < gcp.z < M.z[end]))
        @warn "Starting point outside boundary: " r0
        return OrbitPath(), stat
    end

    stat.ri = r0
    hc = HamiltonianCoordinate(M, gcp)
    gc_ode = make_gc_ode(M,hc,stat)
    stat.vi = gc_ode(r0,false,0.0)

    tspan = (zero(gcp.r),one(gcp.r)*tmax*1e-6)
    ode_prob = ODEProblem(gc_ode,r0,tspan,one_transit)

    if wall != nothing
        wall_cb = wall_callback(wall)
    end
    if one_transit
        cb = transit_callback
        if wall != nothing
            cb = CallbackSet(cb.continuous_callbacks..., wall_cb, cb.discrete_callbacks...)
        end
    else
        cb = oob_cb
        if wall != nothing
            cb = CallbackSet(wall_cb,oob_cb)
        end
    end

    dts = dt*1e-6
    success = false
    retcode = :TotalFailure
    try
        sol = solve(ode_prob, integrator, dt=dts, reltol=1e-8, abstol=1e-12, verbose=false,
                    callback=cb,adaptive=adaptive,save_everystep=store_path)
        success = sol.retcode == :Success || sol.retcode == :Terminated
        retcode = sol.retcode
    catch err
        verbose && println(err)
        if isa(err,InterruptException)
            throw(err)
        end
    end

    if !success && adaptive #Try non-adaptive
        try
            sol = solve(ode_prob, integrator, dt=dts, reltol=1e-8, abstol=1e-12, verbose=false,
                        callback=cb,adaptive=false,save_everystep=store_path)
            success = sol.retcode == :Success || sol.retcode == :Terminated
            retcode = sol.retcode
        catch err
            verbose && println(err)
            if isa(err,InterruptException)
                throw(err)
            end
        end
    end

    for i=1:maxiter #Try progressivly smaller time step with symplectic integrator
        success && break
        try
            sol = solve(ode_prob, ImplicitMidpoint(autodiff=autodiff), dt=dts, reltol=1e-8, abstol=1e-12, verbose=false,
                        callback=cb,save_everystep=store_path)
            success = sol.retcode == :Success || sol.retcode == :Terminated
            retcode = sol.retcode
        catch err
            verbose && println(err)
            if isa(err,InterruptException)
                throw(err)
            end
        end
        dts = dts/10
    end

    if !success
        @warn "Unable to find Guiding Center Orbit" gcp retcode
        stat.class = :incomplete
        return OrbitPath(), stat
    end
    stat.errcode=0

    if one_transit && stat.class != :lost && !stat.poloidal_complete #Try one more time
        try
            sol = solve(ode_prob, ImplicitMidpoint(autodiff=autodiff), dt=dts/10, reltol=1e-8, abstol=1e-12, verbose=false,
                        callback=cb)
            success = sol.retcode == :Success || sol.retcode == :Terminated
            retcode = sol.retcode
        catch err
            verbose && println(err)
            if isa(err,InterruptException)
                throw(err)
            end
        end
    end

    if one_transit && stat.class != :lost && !stat.poloidal_complete
        @warn "Orbit did not complete one transit in allotted time" gcp tmax retcode
    end

    if interp_dt > 0.0 && store_path
        n = floor(Int,sol.t[end]/(interp_dt*1e-6))
        if n > 10
            sol = sol(range(0.0,stop=sol.t[end],length=n))
        end
    else
        n = length(sol)
    end

    #Needed for classification
    r = getindex.(sol.u,1)
    z = getindex.(sol.u,3)
    pitch = get_pitch(M, hc, r, z)

    if stat.class == :unknown && store_path && classify_orbit && one_transit
        stat.class = classify(r,z,pitch,M.axis)
    end

    if !store_path
        return OrbitPath(), stat
    end

    # Get everything else in path
    phi = getindex.(sol.u,2)
    dt = eltype(r)[(sol.t[min(i+1,n)] - sol.t[i]) for i=1:n]
    energy = get_kinetic_energy(M, hc, r, z)

    # P_rz
#    prz = zeros(n)
#    dl = zeros(n)
#    @inbounds for i=1:n
#        v_gc = gc_ode([r[i],phi[i],z[i]],false,0.0)
#        v = norm(v_gc[1:2:3])
##        prz[i] = 1/(T_p*v)
#        dl[i] = v*dt[i]
#    end

    path = OrbitPath(r,z,phi,pitch,energy,dt)

    if stat.class == :lost
        stat.tau_p = zero(stat.tau_p)
        stat.tau_t = zero(stat.tau_t)
        return path, stat
    end

    return path, stat
end

function integrate(M::AxisymmetricEquilibrium, gcp::GCParticle; dt=0.1, tmax=1000.0, integrator=Tsit5(),
                   wall=nothing, interp_dt = 0.1, classify_orbit=true,one_transit=false, store_path=true,
                   maxiter=3, adaptive=true, autodiff=true, verbose=false)

    path, stat = integrate(M, gcp, dt, tmax, integrator, wall, interp_dt,
                           classify_orbit, one_transit, store_path, maxiter,
                           adaptive, autodiff,verbose)

    return path, stat
end

function integrate(M::AxisymmetricEquilibrium, c::EPRCoordinate; kwargs...)
    gcp = GCParticle(c)
    return integrate(M, gcp; kwargs...)
end

function get_orbit(M::AxisymmetricEquilibrium, gcp::GCParticle; kwargs...)
    path, stat = integrate(M, gcp; one_transit=true, kwargs...)

    if stat.class == :incomplete || stat.class == :lost
        return Orbit(EPRCoordinate(typeof(gcp.r)),stat.class,stat.tau_p,stat.tau_t,path)
    end
    hc = HamiltonianCoordinate(M,gcp)
    KE = get_kinetic_energy(M,hc,stat.rm,stat.zm)
    c = EPRCoordinate(KE,stat.pm,stat.rm,stat.zm)

    return Orbit(c,stat.class,stat.tau_p,stat.tau_t,path)
end

function get_orbit(M::AxisymmetricEquilibrium, c::EPRCoordinate; kwargs...)
    #work around integrate 1ms to get away from rmax
    gcp = GCParticle(c.energy,c.pitch,c.r,c.z,c.m,c.q)
    path, stat = integrate(M, gcp; tmax=0.2,interp_dt=0.0)

    gcp = GCParticle(path.energy[end],path.pitch[end],path.r[end],path.z[end],c.m,c.q)
    path, stat = integrate(M, gcp; one_transit=true, kwargs...)
    if stat.class != :incomplete && stat.class != :lost
        stat.rm > c.r && !isapprox(stat.rm,c.r,rtol=1e-4) && (stat.class = :degenerate)
    else
        stat.tau_p=zero(stat.tau_p)
        stat.tau_t=zero(stat.tau_t)
    end
    return Orbit(c,stat.class,stat.tau_p,stat.tau_t,path)
end

function Base.show(io::IO, orbit::Orbit)
    classes = Dict(:trapped=>"Trapped ",:co_passing=>"Co-passing ",:ctr_passing=>"Counter-passing ",
                   :stagnation=>"Stagnation ",:potato=>"Potato ",:incomplete=>"Incomplete ",
                   :degenerate=>"Degenerate ",:meta=>"Meta ",:lost=>"Lost ",:unknown=>"Unknown ")
    class_str = orbit.class in keys(classes) ? classes[orbit.class] : string(orbit.class)

    println(io, class_str*"Orbit:")
    @printf(io, " τₚ = %.3f μs\n", orbit.tau_p*1e6)
    @printf(io, " τₜ = %.3f μs", orbit.tau_t*1e6)
end
