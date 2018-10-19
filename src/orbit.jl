struct OrbitPath{T}
    r::Vector{T}
    z::Vector{T}
    phi::Vector{T}
    pitch::Vector{T}
    energy::Vector{T}
    dt::Vector{T}
    dl::Vector{T}
end

function OrbitPath(T::DataType=Float64)
    OrbitPath(T[],T[],T[],T[],T[],T[],T[])
end

Base.length(op::OrbitPath) = length(op.r)

struct Orbit{T,S<:AbstractOrbitCoordinate{Float64}}
    coordinate::S
    class::Symbol
    tau_p::T
    tau_t::T
end

function Orbit(c::AbstractOrbitCoordinate{T},class=:incomplete) where {T}
    return Orbit(c, class, zero(T), zero(T))
end

mutable struct OrbitStatus
    errcode::Int
    tau_c::Float64
    ri::SArray{Tuple{3},Float64,1,3}
    vi::SArray{Tuple{3},Float64,1,3}
    initial_dir::Int
    nr::Int
    nphi::Int
    naxis::Int
    rm::Float64
    zm::Float64
    pm::Float64
    tm::Float64
    tau_p::Float64
    tau_t::Float64
    poloidal_complete::Bool
    hits_boundary::Bool
    class::Symbol
end

OrbitStatus() = OrbitStatus(1,0.0,SVector{3}(0.0,0.0,0.0),SVector{3}(0.0,0.0,0.0),0,0,0,0,0.0,0.0,0.0,0.0,0.0,0.0,false,false,:incomplete)

function make_gc_ode(M::AxisymmetricEquilibrium, c::T, os::OrbitStatus) where {T<:AbstractOrbitCoordinate}
    oc = HamiltonianCoordinate(M, c)
    ode = function f(y,p::Bool,t)
        r = y[1]
        z = y[3]

        psi = M.psi_rz(r,z)
        g = M.g(psi)

        F = fields(M,r,z)
        B = F.B
        E = F.E

        babs = M.b(r,z)
        gradB = Interpolations.gradient(M.b,r,z)
        gradB = SVector{3}(gradB[1],0.0,gradB[2])
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
        v_gc = (vpara*B/babs + v_exb + v_grad + v_curv)*os.tau_c
        return SVector{3}(v_gc[1], v_gc[2]/r, v_gc[3])
    end
    return ode
end

function get_orbit(M::AxisymmetricEquilibrium, gcp::GCParticle,
                   dt, tmax, integrator, one_transit::Bool, store_path::Bool, maxiter::Int)

    r0 = @SVector [gcp.r,0.0,gcp.z]
    if !((M.r[1] < gcp.r < M.r[end]) && (M.z[1] < gcp.z < M.z[end]))
        @error "Starting point outside boundary: " r0
    end

    hc = HamiltonianCoordinate(M, gcp)

    os = OrbitStatus()
    #os.tau_c = ion_cyclotron_period(M,gcp)
    os.tau_c = 1.0
    os.ri = r0
    gc_ode = make_gc_ode(M,hc,os)
    os.vi = gc_ode(r0,false,0.0)
    os.initial_dir = abs(os.vi[1]) > abs(os.vi[3]) ? 1 : 3

    tspan = (0.0,(tmax*1e-6)/os.tau_c)
    ode_prob = ODEProblem(gc_ode,r0,tspan,one_transit)

    try
        sol = solve(ode_prob, integrator, dt=(dt*1e-6)/os.tau_c, reltol=1e-8, abstol=1e-12, verbose=false,
                    callback=standard_callback, save_everystep=store_path,adaptive=false)
        dtt = (dt/os.tau_c)*1e-6
        for i=1:maxiter
            (sol.retcode == :Success && os.class != :incomplete) && break
            dtt = dtt/10
            sol = solve(ode_prob, integrator, dt=dtt, reltol=1e-8, abstol=1e-12, verbose=false,
                        callback=standard_callback, save_everystep=store_path,adaptive=false)
        end
    catch err
        if !isa(err,InterruptException)
            println(err)
            @show gcp
            return OrbitPath(), os
        else
            throw(err)
        end
    end

    if sol.retcode != :Success
        @warn "Unable to find Guiding Center Orbit" gcp
        os.class = :incomplete
        return OrbitPath(), os
    end
    os.errcode=0

    if !store_path
        return OrbitPath(), os
    end

    n = floor(Int,sol.t[end]/((dt*1e-6)/os.tau_c))
    sol = sol(range(0.0,stop=sol.t[end],length=n))

    r = getindex.(sol.u,1)
    phi = getindex.(sol.u,2)
    z = getindex.(sol.u,3)
    dt = [(sol.t[min(i+1,n)] - sol.t[i])*os.tau_c for i=1:n]

    # P_rz
#    prz = zeros(n)
    dl = zeros(n)
    @inbounds for i=1:n
        v_gc = gc_ode([r[i],phi[i],z[i]],false,0.0)/os.tau_c
        v = norm(v_gc[1:2:3])
#        prz[i] = 1/(T_p*v)
        dl[i] = v*dt[i]
    end

    pitch = get_pitch(M, hc, r, z)
    energy = get_kinetic_energy(M, hc, r, z)

    path = OrbitPath(r,z,phi,pitch,energy,dt,dl)

    if os.class == :lost
        return path, os
    end

    #class = classify(path,M.axis)
    #if os.class == :unknown
    #    os.class = class
    #else
    #    if os.class != class
    #        @warn "Classification disagrees" os.class class gcp
    #    end
    #    os.class = class
    #end
    if os.class == :unknown
        os.class = classify(path,M.axis)
    end

    return path, os
end

function get_orbit(M::AxisymmetricEquilibrium, gcp::GCParticle; dt=0.1, tmax=1000.0, integrator=ImplicitMidpoint(),
                   one_transit=true, store_path=true,maxiter=3)
    path, os = get_orbit(M, gcp, dt, tmax, integrator, one_transit,store_path,maxiter)
    return path, os
end

function get_orbit(M::AxisymmetricEquilibrium, c::EPRCoordinate; dt=0.1, tmax=1000.0, integrator=ImplicitMidpoint(),
                   one_transit=true, store_path=true,maxiter=3)
    gcp = GCParticle(c.energy,c.pitch,c.r,c.z,c.m,c.q)
    path, os = get_orbit(M, gcp, dt, tmax, integrator, one_transit, store_path,maxiter)
    if os.class != :incomplete || os.class != :lost
        os.rm > c.r && (os.class = :degenerate)
    end
    return path, os
end

function down_sample(p::OrbitPath{T}; mean_dl=2.5, nmin=30) where {T}
    L = sum(p.dl)*100 # cm
    npart = length(p)
    npart == 0 && error("Orbit path has zero length")
    np = max(round(Int,L/float(mean_dl)),nmin)
    grps = partition(1:npart, ceil(Int,npart/float(np)))
    ngrps = length(grps)
    r = zeros(T,ngrps)
    z = zeros(T,ngrps)
    phi = zeros(T,ngrps)
    pitch = zeros(T,ngrps)
    energy = zeros(T,ngrps)
    dt = zeros(T,ngrps)
    dl = zeros(T,ngrps)
    for (i, inds) in enumerate(grps)
        fi = inds[1]
        r[i] = p.r[fi]
        z[i] = p.z[fi]
        phi[i] = p.phi[fi]
        pitch[i] = p.pitch[fi]
        energy[i] = p.energy[fi]
        dt[i] = sum(p.dt[inds])
        dl[i] = sum(p.dl[inds])
    end
    return OrbitPath(r,z,phi,pitch,energy,dt,dl)
end

function Base.show(io::IO, orbit::Orbit)
    classes = Dict(:trapped=>"Trapped ",:co_passing=>"Co-passing ",:ctr_passing=>"Counter-passing ",
                   :stagnation=>"Stagnation ",:potato=>"Potato ",:incomplete=>"Incomplete ",
                   :degenerate=>"Degenerate ",:meta=>"Meta ",:lost=>"Lost ")
    class_str = orbit.class in keys(classes) ? classes[orbit.class] : string(orbit.class)

    println(io, class_str*"Orbit Type:")
    @printf(io, " τₚ = %.3f μs\n", orbit.tau_p*1e6)
    @printf(io, " τₜ = %.3f μs", orbit.tau_t*1e6)
end
