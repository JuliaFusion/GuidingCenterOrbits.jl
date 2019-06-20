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
    ri::SArray{Tuple{5},T,1,5}
    vi::SArray{Tuple{5},T,1,5}
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
    return GCStatus(1,SVector{5}(z,z,z,z,z),SVector{5}(z,z,z,z,z),0,0,z,z,z,z,z,z,false,false,:incomplete)
end

@inline function gc_velocity(M::AxisymmetricEquilibrium, gcp::GCParticle, r, z, p_para, mu)
    # Phys. Plasmas 14, 092107 (2007); https://doi.org/10.1063/1.2773702
    # NOTE: Equations 13 and 15 are incorrect. Remove c factor to correct
    q = gcp.q
    m = gcp.m
    mc2 = m*c0*c0

    F = fields(M,r,z)
    Babs = norm(F.B)
    bhat = F.B/Babs

    gB = gradB(M,r,z)
    cB = curlB(M,r,z) #Jfield(M,r,z)*4pi*10^-7
    bhatXgradB = cross(bhat,gB)
    cbhat = (cB + bhatXgradB)/Babs

    Bstar = F.B + (p_para/q)*cbhat
    Bstar_para = dot(Bstar,bhat)

    gamma = sqrt(1 + (2*mu*Babs)/mc2 + (p_para/(m*c0))^2)
    gradg = (mu/(gamma*mc2))*gB
    Estar = F.E - mc2*gradg/q

    Xdot = ((p_para/(gamma*m))*Bstar + cross(Estar,bhat))/Bstar_para

    p_para_dot = dot(q*Estar,Bstar/Bstar_para)

    return SVector{5}(Xdot[1], Xdot[2]/r, Xdot[3], p_para_dot, zero(p_para_dot))
end

function make_gc_ode(M::AxisymmetricEquilibrium, gcp::GCParticle, stat::GCStatus)
    ode = function f(y,p::Bool,t)
        stat
        v_gc = gc_velocity(M, gcp, y[1], y[3], y[4], y[5])
    end
    return ode
end

function integrate(M::AxisymmetricEquilibrium, gcp::GCParticle, phi0,
                   dt, tmin, tmax, integrator, wall::Union{Nothing,Limiter}, interp_dt, classify_orbit::Bool,
                   one_transit::Bool, store_path::Bool, maxiter::Int, adaptive::Bool, autodiff::Bool,
                   r_callback::Bool,verbose::Bool)

    stat = GCStatus(typeof(gcp.r))

    if !((M.r[1] < gcp.r < M.r[end]) && (M.z[1] < gcp.z < M.z[end]))
        verbose && @warn "Starting point outside boundary: " r0 = (gcp.r, gcp.z)
        return OrbitPath(), stat
    end

    # Initial Conditions
    mc2 = gcp.m*c0^2
    p0 = sqrt(((gcp.energy + mc2)^2 - mc2^2)/(c0*c0))
    p_para0 = p0*gcp.pitch*M.sigma
    p_perp0 = p0*sqrt(1-gcp.pitch^2)
    mu_0 = (p_perp0^2)/(2*gcp.m*M.b(gcp.r,gcp.z))
    r0 = SVector{5}(gcp.r, one(gcp.r)*phi0, gcp.z, one(gcp.r)*p_para0, one(gcp.r)*mu_0)

    stat.ri = r0

    gc_ode = make_gc_ode(M,gcp,stat)
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
                    callback=cb,adaptive=adaptive)
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
                        callback=cb,adaptive=false)
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
                        callback=cb, force_dtmin=true)
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

    if interp_dt > 0.0
        n = floor(Int,abs(sol.t[end]/(interp_dt*1e-6)))
        if n > 10
            sol = sol(range(tmin,sol.t[end],length=n))
        end
    end
    n = length(sol)

    #Needed for classification
    r = getindex.(sol.u,1)
    phi = getindex.(sol.u,2)
    z = getindex.(sol.u,3)
    ppara = getindex.(sol.u,4)
    mu = getindex.(sol.u,5)
    pitch = get_pitch.(M, gcp, ppara, mu, r, z)

    if !r_callback
        ind = argmax(r)
        if r[ind] > stat.rm
            stat.rm = r[ind]
            stat.zm = z[ind]
            stat.pm = pitch[ind]
        end
    end

    if stat.class == :unknown && classify_orbit && one_transit
        stat.class = classify(r,z,pitch,M.axis)
    end

    if !store_path
        return OrbitPath(), stat
    end

    # Get everything else in path
    dt = eltype(r)[(sol.t[min(i+1,n)] - sol.t[i]) for i=1:n]
    energy = get_kinetic_energy.(M, gcp, ppara, mu, r, z)

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
    tm = (stat.tau_p - stat.tm)/stat.tau_p
    c = EPRCoordinate(KE,stat.pm,stat.rm,stat.zm,t=tm)

    return Orbit(c,stat.class,stat.tau_p,stat.tau_t,path)
end

function get_orbit(M::AxisymmetricEquilibrium, c::EPRCoordinate; hard = false, rtol=1e-4, kwargs...)
    gcp = GCParticle(c)
    path, stat = integrate(M, gcp; one_transit=true, r_callback=false, kwargs...)
    rmax = stat.rm
    if stat.class != :incomplete && stat.class != :lost
        if rmax > c.r && (hard || !isapprox(rmax,c.r,rtol=rtol))
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
