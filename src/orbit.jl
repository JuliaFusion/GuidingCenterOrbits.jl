struct OrbitPath{T}
    energy::Vector{T}
    pitch::Vector{T}
    r::Vector{T}
    z::Vector{T}
    phi::Vector{T}
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

Base.length(o::Orbit) = length(o.path.r)

mutable struct GCStatus{T<:Number}
    errcode::Int
    ri::SArray{Tuple{3},T,1,3}
    vi::SArray{Tuple{3},T,1,3}
    initial_dir::Int
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
    return GCStatus(1,SVector{3}(z,z,z),SVector{3}(z,z,z),0,0,z,z,z,z,z,z,false,false,:incomplete)
end

@inline function gc_velocity(M::AxisymmetricEquilibrium, hc::HamiltonianCoordinate, r, z)
    F = fields(M,r,z)
    psi = F.psi
    g = F.g
    B = F.B
    E = F.E

    babs = norm(B)
    gradB = Interpolations.gradient(M.b,r,z)
    gradB = SVector{3}(gradB[1],zero(gradB[1]),gradB[2])
    Wperp = hc.mu*babs
    vpara = -babs*(hc.p_phi - hc.q*e0*psi)/(hc.m*g)
    Wpara = 0.5*hc.m*vpara^2

    # ExB Drift
    v_exb = cross(E, B)/babs^2

    # GradB Drift
    b1 = cross(B,gradB)/(babs^3)
    v_grad = Wperp*b1/(hc.q*e0)

    # Curvature Drift
    v_curv = 2*Wpara*b1/(hc.q*e0)

    # Guiding Center Velocity
    v_gc = (vpara*B/babs + v_exb + v_grad + v_curv)
    return SVector{3}(v_gc[1], v_gc[2]/r, v_gc[3])
end

function make_gc_ode(M::AxisymmetricEquilibrium, c::T, stat::GCStatus) where {T<:AbstractOrbitCoordinate}
    oc = HamiltonianCoordinate(M, c)
    ode = function f(y,p::Bool,t)
        stat
        v_gc = gc_velocity(M, oc, y[1], y[3])
    end
    return ode
end

function integrate(M::AxisymmetricEquilibrium, gcp::GCParticle, phi0,
                   dt, tmin, tmax, integrator, wall::Union{Nothing,Limiter}, interp_dt, classify_orbit::Bool,
                   one_transit::Bool, store_path::Bool, maxiter::Int, adaptive::Bool, autodiff::Bool,
                   r_callback::Bool,verbose::Bool)

    stat = GCStatus(typeof(gcp.r))

    r0 = @SVector [gcp.r,one(gcp.r)*phi0,gcp.z]
    if !((M.r[1] < gcp.r < M.r[end]) && (M.z[1] < gcp.z < M.z[end]))
        verbose && @warn "Starting point outside boundary: " r0
        return OrbitPath(), stat
    end

    stat.ri = r0
    hc = HamiltonianCoordinate(M, gcp)
    gc_ode = make_gc_ode(M,hc,stat)
    stat.vi = gc_ode(r0,false,0.0)
    stat.initial_dir = abs(stat.vi[1]) > abs(stat.vi[3]) ? 1 : 3

    tspan = (one(gcp.r)*tmin*1e-6,one(gcp.r)*tmax*1e-6)
    ode_prob = ODEProblem(gc_ode,r0,tspan,one_transit)

    if wall != nothing
        wall_cb = wall_callback(wall)
    end
    if one_transit
        if r_callback
            cb = transit_callback
        else
            stat.nr = 2
            stat.rm = r0[1]
            stat.zm = r0[3]
            cb = CallbackSet(pol_cb,oob_cb)
        end
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
        sol = solve(ode_prob, integrator, dt=dts, reltol=1e-8, abstol=1e-12, verbose=false, force_dtmin=true,
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
            sol = solve(ode_prob, integrator, dt=dts, reltol=1e-8, abstol=1e-12, verbose=false, force_dtmin=true,
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
                        callback=cb,save_everystep=store_path, force_dtmin=true)
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
        verbose && @warn "Unable to find Guiding Center Orbit" gcp retcode
        stat.class = :incomplete
        return OrbitPath(), stat
    end
    stat.errcode=0

    if one_transit && stat.class != :lost && !stat.poloidal_complete #Try one more time
        try
            sol = solve(ode_prob, ImplicitMidpoint(autodiff=autodiff), dt=dts/10, reltol=1e-8, abstol=1e-12, verbose=false,
                        callback=cb, force_dtmin=true)
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
        verbose && @warn "Orbit did not complete one transit in allotted time" gcp tmax retcode
    end

    if interp_dt > 0.0 && store_path
        n = floor(Int,abs(sol.t[end]/(interp_dt*1e-6)))
        if n > 10
            sol = sol(range(tmin,sol.t[end],length=n))
        end
    end
    n = length(sol)

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

    path = OrbitPath(energy,pitch,r,z,phi,dt)

    if stat.class == :lost
        stat.tau_p = zero(stat.tau_p)
        stat.tau_t = zero(stat.tau_t)
        return path, stat
    end

    return path, stat
end

function integrate(M::AxisymmetricEquilibrium, gcp::GCParticle; phi0=0.0, dt=0.1, tmin=0.0,tmax=1000.0,
                   integrator=Tsit5(), wall=nothing, interp_dt = 0.1, classify_orbit=true,
                   one_transit=false, store_path=true,maxiter=3,adaptive=true, autodiff=true,
                   r_callback=true, verbose=false)

    path, stat = integrate(M, gcp, phi0, dt, tmin, tmax, integrator, wall, interp_dt,
                           classify_orbit, one_transit, store_path, maxiter,
                           adaptive, autodiff, r_callback, verbose)
    return path, stat
end

function integrate(M::AxisymmetricEquilibrium, c::EPRCoordinate; kwargs...)
    gcp = GCParticle(c)
    return integrate(M, gcp; r_callback=false, kwargs...)
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
    gcp = GCParticle(c)
    path, stat = integrate(M, gcp; one_transit=true, r_callback=false, kwargs...)
    rmax = maximum(path.r)
    if stat.class != :incomplete && stat.class != :lost
        if rmax > c.r && !isapprox(rmax,c.r,rtol=1e-4)
            stat.class = :degenerate
        end
    else
        stat.tau_p=zero(stat.tau_p)
        stat.tau_t=zero(stat.tau_t)
    end
    return Orbit(c,stat.class,stat.tau_p,stat.tau_t,path)
end

function Base.show(io::IO, orbit::Orbit)
    classes = Dict(:trapped=>"Trapped ",:co_passing=>"Co-passing ",:ctr_passing=>"Counter-passing ",
                   :stagnation=>"Stagnation ",:potato=>"Potato ",:incomplete=>"Incomplete ",
                   :degenerate=>"Degenerate ",:meta=>"Meta ",:lost=>"Lost ")
    class_str = orbit.class in keys(classes) ? classes[orbit.class] : string(orbit.class)

    println(io, class_str*"Orbit:")
    @printf(io, " τₚ = %.3f μs\n", orbit.tau_p*1e6)
    @printf(io, " τₜ = %.3f μs", orbit.tau_t*1e6)
end
