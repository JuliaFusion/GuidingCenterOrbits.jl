"""
A struct that fully describes the guiding-centre motion.

The fields are:\\
`vacuum` - Determines if the path was computed assuming a vacuum or not\\
`drift` - If true, then drift motion was included when computing path\\
`energy` - The energy array of the path\\
`pitch` - The pitch array of the path\\
`r` - The r array of the path\\
`z` - The z array of the path\\
`phi` - The phi array of the path\\
`dt` - The incremental times (array) of the path
"""
struct OrbitPath{T}
    vacuum::Bool
    drift::Bool
    energy::Vector{T}
    pitch::Vector{T}
    r::Vector{T}
    z::Vector{T}
    phi::Vector{T}
    dt::Vector{T}
#    dl::Vector{T}
end

function OrbitPath(T::DataType=Float64; vacuum=false,drift=false)
    OrbitPath(vacuum,drift,T[],T[],T[],T[],T[],T[])
end

Base.length(op::OrbitPath) = length(op.r)

function Base.show(io::IO, path::OrbitPath)
    print(io, typeof(path),"\n")
    print(io, " length = $(length(path))\n")
    print(io, " vacuum = $(path.vacuum)\n")
    print(io, " drift = $(path.vacuum)")
end

"""
A struct that describes a tokamak orbit.

The fields are:\\
`coordinate` - The coordinate of the orbit. Please see coordinates.jl.\\
`class` - The class of the orbit (co-passing, trapped, potato etc).\\
`tau_p` - The poloidal transit time of the orbit.\\
`tau_t` - The toroidal transit time of the orbit.\\
`path` - The OrbitPath of the orbit.
"""
struct Orbit{T,S<:AbstractOrbitCoordinate{Float64}}
    coordinate::S
    class::Symbol
    tau_p::T
    tau_t::T
    path::OrbitPath{T}
    gcvalid::Bool
end

function Orbit(c::AbstractOrbitCoordinate{T},class=:incomplete) where {T}
    return Orbit(c, class, zero(T), zero(T), OrbitPath(T), false)
end

Base.length(o::Orbit) = length(o.path.r)

"""
A struct that fully describes the integration status of the guiding-centre path/motion.

The fields are:\\
`errcode` - The error code of the integration\\
`ri` - The initial guiding-centre vector (r0,phi0,z0,p_para0,mu0)\\
`vi` - The initial guiding-centre vector (dr,dphi,dz,dp_para,dmu)\\
`initial_dir` - The initial direction of the guiding-centre motion\\
`nr` - A process number (please see integrate() and callbacks.jl)\\
`rm` - The rm coordinate\\
`zm` - The zm coordinate\\
`pm` - The pm coordinate\\
`tm` - The time at the rm,zm coordinate\\
`tau_p` - The poloidal transit time\\
`tau_t` - The toroidal transit time\\
`poloidal_complete` - True if the guiding-centre particle has completed one orbit poloidally\\
`hits_boundary` - True if the guiding-centre particle hit the boundaries set by M (AbstractEquilibrum)\\
`class` - The type of orbit (:co-passing, :trapped, :potato etc)
"""
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

function reset!(stat::GCStatus{T}, one_transit, r_callback) where T
    stat.errcode = 1
    stat.nr = 0
    stat.rm = zero(T)
    stat.zm = zero(T)
    stat.tm = zero(T)
    stat.tau_p = zero(T)
    stat.tau_t = zero(T)
    stat.poloidal_complete = false
    stat.hits_boundary = false
    stat.class = :incomplete
    if one_transit && !r_callback
        stat.nr = 2
        stat.rm = stat.ri[1]
        stat.zm = stat.ri[3]
    end
end

"""
    gc_velocity(M, gcp, r, z, p_para, mu, vacuum, drift)

Calculate the velocity of the guiding-centre particle gcp, given the inputs. The equations
are from https://doi.org/10.1063/1.2773702, with a correction for a factor c.
"""
@inline function gc_velocity(M::AbstractEquilibrium, gcp::GCParticle, r, z, p_para, mu, vacuum::Bool, drift::Bool)
    # Phys. Plasmas 14, 092107 (2007); https://doi.org/10.1063/1.2773702
    # NOTE: Equations 13 and 15 are incorrect. Remove c factor to correct
    q = gcp.q*e0 # gcp.q is given in units of integer number of elementary charges. e0 is the elementary charge.
    m = gcp.m
    mc2 = m*c0*c0

    B, E = fields(M,r,z) # Please check Equilibrium.jl/equil.jl/fields()
    Babs = norm(B)
    bhat = B/Babs

    if drift
        gB = gradB(M,r,z) # Please check Equilibrium.jl/equil.jl/gradB()
        if !vacuum
            cB = curlB(M,r,z) #Jfield(M,r,z)*4pi*10^-7
        else
            cB = SVector{3}(0.0,0.0,0.0)
        end
        bhatXgradB = cross(bhat,gB)
        cbhat = (cB + bhatXgradB)/Babs

        Bstar = B + (p_para/q)*cbhat
    else
        Bstar = B
    end
    Bstar_para = dot(Bstar,bhat)

    gamma = sqrt(1 + (2*mu*Babs)/mc2 + (p_para/(m*c0))^2) # The relativistic gamma factor
    if drift
        gradg = (mu/(gamma*mc2))*gB
        Estar = E - mc2*gradg/q
        Xdot = ((p_para/(gamma*m))*Bstar + cross(Estar,bhat))/Bstar_para
    else
        Estar = E
        Xdot = ((p_para/(gamma*m))*Bstar)/Bstar_para
    end

    p_para_dot = dot(q*Estar,Bstar/Bstar_para)

    Xdot = cylindrical_cocos(cocos(M), Xdot[1], Xdot[2]/r, Xdot[3])
    return SVector{5}(Xdot[1], Xdot[2], Xdot[3], p_para_dot, zero(p_para_dot))
end

"""
    make_gc_ode(M, gcp, stat, vacuum, drift)

Make and return an ordinary differential equation function for the guiding-centre motion.
"""
function make_gc_ode(M::AbstractEquilibrium, gcp::GCParticle, stat::GCStatus,vacuum::Bool,drift::Bool)
    ir::Int, iphi::Int, iz::Int = cylindrical_cocos_indices(cocos(M))
    ode = function f(y,p::Bool,t)
        stat
        v_gc = gc_velocity(M, gcp, y[ir], y[iz], y[4], y[5],vacuum,drift)
    end
    return ode
end

"""
    integrate(M, gcp, phi0, ..., maxphi, debug)

Integrate the guiding-centre particle motion given the axisymmetric equilibrium M
and lots of input.

Inputs:\\
`M` - The AbstractEquilibrium struct containing all info about the tokamak equilibrium. Please see Equilibrium.jl/equil.jl.\\
`gcp` - The guiding-centre particle whose motion is to be integrated. Please see GuidingCenterParticles.jl/particles.jl.\\
`phi0` - The initial toroidal angle of the guiding-centre particle\\
`dt` - The (initial) incremental time step of the integration (s)\\
`tmin` - The starting time for the integration (s)\\
`tmax` - The stopping time for the integration (s)\\
`integrator` - The type of ODE integrator to be used. For example Tsit5(), BS3(), Vern9() etc. See DifferentialEquations.jl for full list\\
`wall` - Boundary object defined in Equilibrium.jl. The (R,z) coordinates of the wall of the tokamak (in meters)\\
`interp_dt` - Interpolate the resulting integrated path onto a path evenly spaced in time, with time step size inter_dt. If 0.0, interpolate onto max_length time points.\\
`classify_orbit` - If true, let the algorithm classify the resulting orbit\\
`one_transit` - If true, stop the integration after particles completes one poloidal transit\\
`store_path` - If true, then the orbit path will be returned. Otherwise, an empty path object is returned with the status object.\\
`max_length` - The solve() function returns path arrays that are arbitrarily long. max_length will be the length of an interpolated path that is returned instead.\\
`maxiter` - If the first adaptive and non-adaptive integration attempt fail, the algorithm will re-try with progresively smaller time steps this number of times.\\
`toa` - Stands for 'try only adaptive'. If true, non-adaptive integration will not be attempted (safe but possibly incomplete).\\
`maxiters` - The solve() function has a default maximum number of iterations of 1e5. Use maxiters to increase this number.\\
`autodiff` - Deprecated.\\
`r_callback` - If true, then the integration will be terminated when the ratio between the r direction speed and the total speed has gone to zero 20 times (prevents infinite loop)\\
`verbose` - If true, lots of output messages will be printed during execution.\\
`vacuum` - If true, then vacuum will be assumed\\
`drift` - If true, then drift effects will be included\\
`limit_phi` - If true, then the integration will be terminated when phi direction reaches the value maxphi\\
`maxphi` - Please see 'limit_phi'\\
`debug` - If true, then de-bugging mode is activated. Function will terminate and return after first adaptive integration (if failed).
"""
function integrate(M::AbstractEquilibrium, gcp::GCParticle, phi0,
                   dt, tmin, tmax, integrator, wall::Union{Nothing,Wall}, interp_dt, classify_orbit::Bool,
                   one_transit::Bool, store_path::Bool, max_length::Int, maxiter::Int, toa::Bool, maxiters::Int, autodiff::Bool,
                   r_callback::Bool,verbose::Bool, vacuum::Bool, drift::Bool, limit_phi::Bool, maxphi, debug::Bool)

    stat = GCStatus(typeof(gcp.r))
    rlims, zlims = limits(M)
    if !((rlims[1] < gcp.r < rlims[end]) && (zlims[1] < gcp.z < zlims[end]))
        verbose && @warn "Starting point outside boundary: " r0 = (gcp.r, gcp.z)
        return OrbitPath(;vacuum=vacuum,drift=drift), stat
    end

    # Initial Conditions
    mc2 = gcp.m*c0^2
    p0 = sqrt(((1e3*e0*gcp.energy + mc2)^2 - mc2^2)/(c0*c0)) # The initial particle momentum
    # move abs(pitch) away from one to avoid issues when calculating jacobian;
    # don't try to be clever by using clamp it doesn't work with the autodiff
    if abs(gcp.pitch) == 1.0
        pitch0 = sign(gcp.pitch)*prevfloat(abs(gcp.pitch))
    else
        pitch0 = gcp.pitch
    end
    p_para0 = p0*pitch0*B0Ip_sign(M) # The initial parallel momentum
    p_perp0 = p0*sqrt(1.0 - pitch0^2) # The initial perpendicular momentum
    mu_0 = (p_perp0^2)/(2*gcp.m*norm(Bfield(M,gcp.r,gcp.z))) # The initial (and presumably constant) magnetic moment

    r0 = SVector{5}(cylindrical_cocos(cocos(M), gcp.r, one(gcp.r)*phi0, gcp.z)...,
                    one(gcp.r)*p_para0, one(gcp.r)*mu_0) # The initial guiding-centre element vector

    stat.ri = r0 # The initial status vector is the initial vector

    ir, iphi, iz = cylindrical_cocos_indices(cocos(M))
    gc_ode = make_gc_ode(M,gcp,stat,vacuum,drift) # Make the ordinary differential equation
    stat.vi = gc_ode(r0,false,tmin) # Obtain the initial velocity

    r_callback = r_callback && !iszero(stat.vi[ir])

    stat.initial_dir = abs(stat.vi[ir]) > abs(stat.vi[iz]) ? ir : iz # Is it mostly in the R- or z-direction? Element number 1 or 3?

    tspan = (one(gcp.r)*tmin,one(gcp.r)*tmax)
    ode_prob = ODEProblem(gc_ode,r0,tspan,one_transit) # Define the ODE problem

    # The callbacks. Please see callbacks.jl for specification.
    if wall != nothing
        wall_cb = wall_callback(wall)
    end
    if one_transit
        if r_callback
            cb = transit_callback
        else
            stat.nr = 2
            stat.rm = r0[ir]
            stat.zm = r0[iz]
            cb = CallbackSet(pol_cb,oob_cb, brr_cb)
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
    if limit_phi
        phi_cb = phi_callback(maxphi)
        cb = CallbackSet(cb.continuous_callbacks..., phi_cb, cb.discrete_callbacks...)
    end

    if verbose
        println(" ")
        println("---- Initiating orbit integration ---- ")
        println("- gcp: $(gcp)")
        println("- dt: $(round(dt/cyclotron_period(M,gcp),sigdigits=4)) times a cyclotron period")
        println("- (tmin,tmax)= ($(tmin),$(round(tmax/dt))*dt)")
        println("- Integrator: $(integrator)")
        println("- interp_dt: $(interp_dt)")
        println("- classify_orbit: $(classify_orbit)")
        println("- one_transit: $(one_transit)")
        println("- Maximum orbit length (will otherwise be interpolated onto this length): $(max_length)")
        println("- The first integration will use $(maxiters) numerical iterations. ")
        r_callback && println("- r_callback will be used.")
        vacuum && println("- Vacuum is assumed.")
        drift && println("- Particle drifts are included.")
        limit_phi && println("- The integration will be toroidally limited to $(round(maxphi/(2*pi))) turns.")
        println("--------------------------------")
    end

    # Always try adaptive first
    dts = dt
    success = false
    retcode = :TotalFailure # I think we should change this to ReturnCode.Failure (or maybe ReturnCode.Default)
    try
        sol = solve(ode_prob, integrator, dt=dts, reltol=1e-8, abstol=1e-12, verbose=false, force_dtmin=true,
                    callback=cb,adaptive=true,maxiters=maxiters)
        if sol.t[end] >= tspan[2] && one_transit
            reset!(stat,one_transit,r_callback)
            ode_prob = remake(ode_prob; tspan=(tspan[1],100*tspan[2]))
            sol = solve(ode_prob, integrator, dt=dts, reltol=1e-8, abstol=1e-12, verbose=false, force_dtmin=true,
                        callback=cb,adaptive=true,maxiters=Int(maxiters*10))
        end
        success = SciMLBase.successful_retcode(sol) &&
                  ((stat.poloidal_complete || !one_transit) || limit_phi) || stat.hits_boundary
        retcode = sol.retcode

    catch err
        verbose && (println("Adaptive"); println(err))
        if isa(err,InterruptException)
            rethrow(err)
        end
    end

    if !success && debug
        verbose && println("Adaptive failed. You wanted to fix bugs. Here is some extra information that might help.")
        verbose && println("Success: $(success)")
        verbose && println("sol.retcode: $(retcode)")
        verbose && println("stat.poloidal_complete: $(stat.poloidal_complete)")
        verbose && println("!one_transit: $(!one_transit)")
        verbose && println("Length(sol): $(length(sol))")
        verbose && println(gcp)
        verbose && println("The integration will terminate now.")

        ir, iphi, iz = cylindrical_cocos_indices(cocos(M))
        r = getindex.(sol.u,ir)
        phi = getindex.(sol.u,iphi)
        z = getindex.(sol.u,iz)
        ppara = getindex.(sol.u,4)
        mu = getindex.(sol.u,5)
        pitch = get_pitch.(M, gcp, ppara, mu, r, z)
        n = length(sol)
        dt = eltype(r)[(sol.t[min(i+1,n)] - sol.t[i]) for i=1:n]
        energy = get_kinetic_energy.(M, gcp, ppara, mu, r, z)

        path = OrbitPath(vacuum,drift,energy,pitch,r,z,phi,dt)
        return path, stat
    end

    if !success && !toa #Try non-adaptive
        verbose && println(" ")
        verbose && println("Adaptive did not find guiding-centre orbit. Trying non-adaptive... ")
        try
            reset!(stat,one_transit,r_callback)
            sol = solve(ode_prob, integrator, dt=dts, reltol=1e-8, abstol=1e-12, verbose=false, force_dtmin=true,
                        callback=cb,adaptive=false)
            if sol.t[end] >= tspan[2] && one_transit
                reset!(stat,one_transit,r_callback)
                ode_prob = remake(ode_prob; tspan=(tspan[1],100*tspan[2]))
                sol = solve(ode_prob, integrator, dt=dts, reltol=1e-8, abstol=1e-12, verbose=false, force_dtmin=true,
                            callback=cb,adaptive=false)
            end
            success = SciMLBase.successful_retcode(sol) &&
                      (stat.poloidal_complete || !one_transit || limit_phi) || stat.hits_boundary
            retcode = sol.retcode
        catch err
            verbose && (println("Non-Adaptive: dt = $(dts)"); println(err))
            if isa(err,InterruptException)
                rethrow(err)
            end
        end
    end

    for i=1:maxiter #Try progressivly smaller time step
        (success || toa) && break
        dts = dts/10
        verbose && println("Re-trying non-adaptive with smaller timestep... ")
        try
            reset!(stat,one_transit,r_callback)
            sol = solve(ode_prob, integrator, dt=dts, reltol=1e-8, abstol=1e-12, verbose=false,
                        callback=cb, force_dtmin=true, adaptive=false)
            success = SciMLBase.successful_retcode(sol) &&
                      (stat.poloidal_complete || !one_transit || limit_phi) || stat.hits_boundary
            retcode = sol.retcode
        catch err
            verbose && (println("Non-Adaptive: dt = $(dts)"); println(err))
            if isa(err,InterruptException)
                rethrow(err)
            end
        end
    end

    if !success
        verbose && @warn "Unable to find Guiding Center Orbit" gcp retcode
        stat.class = :incomplete
        return OrbitPath(;vacuum=vacuum,drift=drift), stat
    end
    stat.errcode=0

    if one_transit && stat.class != :lost && !stat.poloidal_complete && !toa #Try one more time
        verbose && println("Trying one last time with even smaller timestep... ")
        try
            reset!(stat,one_transit,r_callback)
            sol = solve(ode_prob, integrator, dt=dts/10, reltol=1e-8, abstol=1e-12, verbose=false,
                        callback=cb, force_dtmin=true, adaptive=false)
            success = SciMLBase.successful_retcode(sol) &&
                      (stat.poloidal_complete || !one_transit) || stat.hits_boundary
            retcode = sol.retcode
        catch err
            verbose && println(err)
            if isa(err,InterruptException)
                rethrow(err)
            end
        end
    end

    if one_transit && stat.class != :lost && !stat.poloidal_complete
        verbose && println(" ")
        verbose && println("maxiters: $(maxiters)")
        verbose && println("Solution array lengths: $(length(sol.t))")
        verbose && @warn "Orbit did not complete one transit in allotted time" gcp tmax retcode
    end

    if interp_dt > 0.0
        n = floor(Int,abs(sol.t[end]/interp_dt))
        if n > max_length
            @warn "Requested interp_dt conflicts with max_length"
        end
        n = min(max_length,n)
        sol = sol(range(tmin,sol.t[end],length=n)) # Interpolate solution onto n points evently spaced in time
    else
        n = min(max_length, length(sol))
        sol = sol(range(tmin,sol.t[end],length=n)) # Interpolate solution onto n points evently spaced in time
    end
    n = length(sol)

    #Needed for classification
    ir, iphi, iz = cylindrical_cocos_indices(cocos(M))
    r = getindex.(sol.u,ir)
    phi = getindex.(sol.u,iphi)
    z = getindex.(sol.u,iz)
    ppara = getindex.(sol.u,4)
    mu = getindex.(sol.u,5)
    pitch = get_pitch.(M, gcp, ppara, mu, r, z)

    if !r_callback # Check rm, if not automatically checked in r_callback
        ind = argmax(r)
        if r[ind] > stat.rm
            stat.rm = r[ind]
            stat.zm = z[ind]
            stat.pm = pitch[ind]
        end
    end

    # Classify the orbit
    if stat.class == :unknown && classify_orbit && one_transit
        stat.class = classify(r,z,pitch,magnetic_axis(M))
    end

    if !store_path
        return OrbitPath(;vacuum=vacuum,drift=drift), stat
    end

    # Get everything else in path
    dt = eltype(r)[(sol.t[min(i+1,n)] - sol.t[i]) for i=1:n] # Get the incremental time steps
    energy = get_kinetic_energy.(M, gcp, ppara, mu, r, z) # Get the energy of the guiding-centre particle

    # P_rz
#    prz = zeros(n)
#    dl = zeros(n)
#    @inbounds for i=1:n
#        v_gc = gc_ode([r[i],phi[i],z[i]],false,0.0)
#        v = norm(v_gc[1:2:3])
##        prz[i] = 1/(T_p*v)
#        dl[i] = v*dt[i]
#    end

    path = OrbitPath(vacuum,drift,energy,pitch,r,z,phi,dt) # Create an OrbitPath object

    if stat.class == :lost # Take care of lost orbits
        stat.tau_p = zero(stat.tau_p)
        stat.tau_t = zero(stat.tau_t)
        return path, stat
    end

    return path, stat # Return the orbit path and the integration status object
end

"""
    integrate(M, gcp; phi0=0.0, ..., debug = false)

All default values for the arguments are set in this function. They are passed into the actual integrate() function
written above.
"""
function integrate(M::AbstractEquilibrium, gcp::GCParticle; phi0=0.0, dt=cyclotron_period(M,gcp)*1e-1,
                   tmin=0.0,tmax=1e6*dt, integrator=Tsit5(), wall=nothing, interp_dt = 0.0,
                   classify_orbit=true, one_transit=false, store_path=true, max_length=500,
                   maxiter=3, toa=false, maxiters=Int(1e6), autodiff=true, r_callback=false, verbose=false,
                   vacuum=false, drift=true, limit_phi=false, maxphi=10*2*pi, debug=false)

    path, stat = integrate(M, gcp, phi0, dt, tmin, tmax, integrator, wall, interp_dt,
                           classify_orbit, one_transit, store_path, max_length, maxiter,
                           toa, maxiters, autodiff, r_callback, verbose, vacuum, drift, limit_phi, maxphi, debug)
    return path, stat
end

function integrate(M::AbstractEquilibrium, c::EPRCoordinate; kwargs...)
    gcp = GCParticle(c)
    return integrate(M, gcp; kwargs...)
end

function get_orbit(M::AbstractEquilibrium, gcp::GCParticle; kwargs...)
    path, stat = integrate(M, gcp; one_transit=true, r_callback=true, kwargs...)
    gcvalid = gcde_check(M, gcp, path)

    if stat.class == :incomplete || stat.class == :lost
        return Orbit(EPRCoordinate(typeof(gcp.r);amu=(gcp.m/mass_u),q=gcp.q),stat.class,stat.tau_p,stat.tau_t,path,gcvalid)
    end
    hc = HamiltonianCoordinate(M,gcp)
    KE = get_kinetic_energy(M,hc,stat.rm,stat.zm)
    tm = (stat.tau_p - stat.tm)/stat.tau_p
    c = EPRCoordinate(KE,stat.pm,stat.rm,stat.zm;t=tm,amu=(gcp.m/mass_u),q=gcp.q)

    return Orbit(c,stat.class,stat.tau_p,stat.tau_t,path,gcvalid)
end

function get_orbit(M::AbstractEquilibrium, c::EPRCoordinate; hard = false, rtol=1e-4, kwargs...)
    gcp = GCParticle(c)
    path, stat = integrate(M, gcp; one_transit=true, r_callback=false, kwargs...)
    gcvalid = gcde_check(M, gcp, path)

    rmax = stat.rm
    if stat.class != :incomplete && stat.class != :lost
        if rmax > c.r && (hard || !isapprox(rmax,c.r,rtol=rtol))
            stat.class = :invalid
        end
    else
        stat.tau_p=zero(stat.tau_p)
        stat.tau_t=zero(stat.tau_t)
    end
    return Orbit(c,stat.class,stat.tau_p,stat.tau_t,path,gcvalid)
end

function Base.show(io::IO, orbit::Orbit)
    classes = Dict(:trapped=>"Trapped ",:co_passing=>"Co-passing ",:ctr_passing=>"Counter-passing ",
                   :stagnation=>"Stagnation ",:potato=>"Potato ",:incomplete=>"Incomplete ",
                   :Invalid=>"Invalid ",:meta=>"Meta ",:lost=>"Lost ")
    class_str = orbit.class in keys(classes) ? classes[orbit.class] : string(orbit.class)

    println(io, class_str*"Orbit:")
    @printf(io, " τ_p= %.3f μs\n", orbit.tau_p*1e6)
    @printf(io, " τ_t= %.3f μs\n", orbit.tau_t*1e6)
    print(io,   " gc-valid: $(orbit.gcvalid)")
end
