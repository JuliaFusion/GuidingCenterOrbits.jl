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
    path::OrbitPath{T}
end

function Orbit(c::AbstractOrbitCoordinate{T},class=:incomplete) where {T}
    return Orbit(c, class, zero(T), zero(T), OrbitPath(T))
end

function Orbit(c::AbstractOrbitCoordinate{T}, op::OrbitPath{T}; class=:incomplete) where {T}
    return Orbit(c, class, zero(T), zero(T), op)
end

function make_gc_ode(M::AxisymmetricEquilibrium, c::T) where {T<:AbstractOrbitCoordinate}
    oc = HamiltonianCoordinate(M, c)
    function vgc(t, y, ydot)
        ri = [y[1],y[3]]
        r = y[1]
        z = y[3]
        psi = M.psi_rz(r,z)
        g = M.g(psi)

        F = fields(M,r,z)
        B = F.B
        E = F.E
        #J = Jfield(M,r,z)

        babs = M.b(r,z)
        gradB = Interpolations.gradient(M.b,r,z)
        gradB = SVector{3}(gradB[1],0.0,gradB[2])

        Wperp = oc.mu*babs
        vpara = -babs*(oc.p_phi - oc.q*e0*psi)/(oc.amu*mass_u*g)
        Wpara = 0.5*oc.amu*mass_u*vpara^2

        # ExB Drift
        v_exb = cross(E, B)/babs^2

        # GradB Drift
        b1 = cross(B,gradB)/(babs^3)
        v_grad = Wperp*b1/(oc.q*e0)

        # Curvature Drift
        #bb = (dot([gradB[1],0.0,gradB[2]],B) - cross(B,mu0*J))/(babs^2)
        #v_curv = 2*Wpara*cross(B,bb)/(babs^2)/(oc.q*e0)
        v_curv = 2*Wpara*b1/(oc.q*e0)

        # Guiding Center Velocity
        v_gc = vpara*B/babs .+ v_exb .+ v_grad .+ v_curv

        ydot[1] = v_gc[1]
        ydot[2] = v_gc[2]/ri[1]
        ydot[3] = v_gc[3]

        return Sundials.CV_SUCCESS
    end
end

function get_orbit(M::AxisymmetricEquilibrium, E, pitch_i, ri, zi, amu, q::Int, nstep::Int, tmax, one_transit)

    hc = HamiltonianCoordinate(M, E, pitch_i, ri, zi, amu=amu, q=q)

    #Function to check if in function
    hits_boundary = [false]
    in_boundary = function ib(x)
        if (M.r[1] < x[1] < M.r[end]) && (M.z[1] < x[3] < M.z[end])
            return true
        else
            hits_boundary[1] = true
            return false
        end
    end

    #Function to check if orbit is complete
    r0 = [ri, 0.0, zi]
    r1 = copy(r0)
    initial_dir_sgn = [0]
    initial_dir = [0]
    npol = [0]
    npol_steps = [0]
    poloidal_complete = [false]
    one_poloidal = function icp(r2)
        poloidal_complete[1] && return false
        npol_steps[1] = npol_steps[1] + 1
        if initial_dir[1] == 0
            r = r2 - r0
            initial_dir[1] = abs(r[1]) > abs(r[3]) ? 1 : 3
            initial_dir_sgn[1] = sign(r[initial_dir[1]])
            r1 .= convert(Vector, r2)
            return false
        end

        dr10 = r1 - r0
        dr20 = r2 - r0
        dr21 = r2 - r1
        ind = initial_dir[1]
        if sign(dr10[ind]*dr20[ind]) < 0
            npol[1] = npol[1] + 1
            if norm(dr10[1:2:3]) < norm(dr21[1:2:3]) && sign(dr20[ind]) == initial_dir_sgn[1]
                poloidal_complete[1] = npol[1] == 2 || npol[1] == 4
                return true
            end
        end
        r1 .= convert(Vector, r2)
        return false
    end

    #Function to call between timesteps
    cb = function callback(mem, t, x)
        !(one_poloidal(x) && one_transit) && in_boundary(x)
    end

    if !in_boundary(r0)
        error("Starting point outside boundary: ", r0)
    end

    f = make_gc_ode(M, hc)
    t = 1e-6*collect(range(0.0,stop=tmax,length=nstep))

    res = Sundials.cvode(f, r0, t, reltol=1e-8,abstol=1e-12, callback = cb)

    r = res[:,1]
    phi = res[:,2]
    z = res[:,3]

    n = length(r)
    dt = fill(t[2]-t[1],n)

    #calculate transits
    ydot = zeros(3)
    nlast = npol_steps[1]
    flag = f(0.0,[r[nlast],phi[nlast],z[nlast]], ydot)
    dtlast = sqrt((r[nlast] - r[1])^2 + (z[nlast] - z[1])^2)/norm(ydot[1:2:3])
    if one_transit
        dt[end] = dtlast
        tau_p = sum(dt)
    else
        tau_p = sum(dt[1:nlast-1]) + dtlast
    end
    phiend = phi[nlast] + ydot[2]*dtlast
    tau_t = 2pi*tau_p/abs(phiend - phi[1])

    # P_rz
#    prz = zeros(n)
    dl = zeros(n)
    @inbounds for i=1:n
        ydot .= 0.0
        flag = f(0.0, res[i,:], ydot)
        v = norm(ydot[1:2:3])
#        prz[i] = 1/(T_p*v)
        dl[i] = v*dt[i]
    end

    pitch = get_pitch(M, hc, r, z)
    energy = get_kinetic_energy(M, hc, r, z)
    rmax, ind = findmax(r)
    pitch_rmax = pitch[ind]
    energy_rmax = energy[ind]
    path = OrbitPath(r,z,phi,pitch,energy,dt,dl)

    c = EPRCoordinate(energy_rmax, pitch_rmax, rmax, z[ind], hc.amu, hc.q)

    if hits_boundary[1]
        return Orbit(c, path, class=:lost)
    end

    if !poloidal_complete[1]
        return Orbit(c, path, class=:incomplete)
    end

    class = classify(path, pitch, M.axis, n=nlast)

    return Orbit(c, class, tau_p, tau_t, path)
end

function get_orbit(M::AxisymmetricEquilibrium, E, p, r, z; amu=H2_amu, q=1, nstep=3000, tmax=500.0, store_path=true, one_transit=true)
    o = get_orbit(M, E, p, r, z, amu, q, nstep, tmax, one_transit)
    if !store_path
        o = Orbit(o.coordinate, o.class, o.tau_p, o.tau_t, OrbitPath(typeof(o.tau_p)))
    end
    return o
end

function get_orbit(M::AxisymmetricEquilibrium, c::EPRCoordinate; nstep=3000, tmax=500.0, store_path=true, one_transit=true)
    o = get_orbit(M, c.energy, c.pitch, c.r, c.z, c.amu, c.q, nstep, tmax, one_transit)
    if o.class != :incomplete || o.class != :lost
        maximum(o.path.r) > c.r && return Orbit(c,:degenerate)
    end
    if !store_path
        o = Orbit(o.coordinate, o.class, o.tau_p, o.tau_t, OrbitPath(typeof(o.tau_p)))
    end
    return o
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
