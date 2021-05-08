function _get_shifted_jacobian(M, o::Orbit; transform = x -> x, spline=true, tol=0.1,verbose=false)
    tm = o.coordinate.t
    tp = range(0.0,1.0,length=length(o.path))
    if spline
        # Create spline of orbit path
        A = hcat(o.path.energy,o.path.pitch,o.path.r,o.path.z)
        itp = extrapolate(scale(interpolate(A,(BSpline(Cubic(Periodic(OnGrid()))), NoInterp())), tp, 1:4), (Periodic(), Throw()))

        # Create new orbit path
        energy = itp.(tp .- tm, 1)
        pitch = itp.(tp .- tm, 2)
        r = itp.(tp .- tm, 3)
        z = itp.(tp .- tm, 4)

        dt = fill(step(tp)*o.tau_p,length(o.path))
        dt[end] = 0.0
        opath = OrbitPath(o.path.vacuum, o.path.drift, energy,pitch,r,z,r*0.0,dt)
    else
        os = get_orbit(M,o.coordinate,drift=o.path.drift,vacuum=o.path.vacuum)
        opath = os.path
    end

    # Calculate jacobian
    J = _get_jacobian(M, o.coordinate, opath, o.tau_p, transform, tol, verbose)

    # Shift jacobian
    tj = range(0.0,1.0,length=length(J))
    Jitp = extrapolate(scale(interpolate(J,BSpline(Linear())),tj), Periodic())
    Jshifted = max.(Jitp.(tp .+ tm),0.0)

    return Jshifted
end

function get_jacobian(M::AbstractEquilibrium, c::EPRCoordinate;
                      transform = x -> x, spline=true, tol=0.1,
                      verbose=false, kwargs...)
    o = get_orbit(M,c; kwargs...)
    if o.class == :invalid
        return _get_shifted_jacobian(M, o; transform=transform, spline=spline,tol=tol,verbose=verbose)
    end
    return _get_jacobian(M, o.coordinate, o.path, o.tau_p, transform, tol,verbose)
end

function get_jacobian(M::AbstractEquilibrium, o::Orbit;
                      transform = x -> x, spline=true, tol=0.1,
                      verbose=false)
    length(o.path) == 0 && return zeros(length(o.path))
    r0 = o.path.r[1]
    if r0 < o.coordinate.r && !isapprox(r0,o.coordinate.r,rtol=1e-4)
        return _get_shifted_jacobian(M, o; transform=transform, spline=spline,tol=tol,verbose=verbose)
    end
    return _get_jacobian(M, o.coordinate, o.path, o.tau_p, transform, tol,verbose)
end

function check_jacobian(J::Array; tol = 0.1)
    !(any(abs.(J .- circshift(J,-1))./median(J) .> tol))
end

function _get_jacobian(M::AbstractEquilibrium, c::Union{GCParticle,AbstractOrbitCoordinate},
                       o::OrbitPath, tau_p, transform, tol, verbose)
    np = length(o)
    detJ = zeros(np)
    if np == 0 || tau_p == 0.0
        return detJ
    end

    ## Do first half of orbit
    ed = Dual(o.energy[1],(1.0,0.0,0.0,0.0))
    pd = Dual(o.pitch[1], (0.0,1.0,0.0,0.0))
    rd = Dual(o.r[1],     (0.0,0.0,1.0,0.0))
    td = Dual(1e-30,      (0.0,0.0,0.0,1.0))
    zd = one(td)*o.z[1]

    gcp = GCParticle(ed,pd,rd,zd,c.m,c.q)
    path, stat = integrate(M, gcp; tmax=td*tau_p, interp_dt=0.0,
                           one_transit=false, store_path=true,
                           classify_orbit=false, drift=o.drift, vacuum=o.vacuum)

    x = transform([path.energy[end],path.pitch[end],path.r[end],path.z[end]])
    detJ[1] = max(abs(det(hcat((partials(xx) for xx in x)...))),0.0)

    ep = partials(path.energy[end])
    pp = partials(path.pitch[end])
    rp = partials(path.r[end])
    zp = partials(path.z[end])
    np1 = floor(Int,np/2)
    for i=2:np1
        ed = Dual(o.energy[i-1],ep)
        pd = Dual(o.pitch[i-1], pp)
        rd = Dual(o.r[i-1],     rp)
        zd = Dual(o.z[i-1],     zp)
        gcp = GCParticle(ed,pd,rd,zd,c.m,c.q)

        dt = o.dt[i-1]
        dt = dt + 1e-30*iszero(dt)
        path, stat = integrate(M, gcp; tmax=dt, interp_dt=0.0,
                               one_transit=false, store_path=true,
                               classify_orbit=false, drift=o.drift,vacuum=o.vacuum)

        x = transform([path.energy[end],path.pitch[end],path.r[end],path.z[end]])
        ep = partials(path.energy[end])
        pp = partials(path.pitch[end])
        rp = partials(path.r[end])
        zp = partials(path.z[end])
        detJ[i] = max(abs(det(hcat((partials(xx) for xx in x)...))),0.0)
    end

    ## Do second half of orbit
    ed = Dual(o.energy[1],(1.0,0.0,0.0,0.0))
    pd = Dual(o.pitch[1], (0.0,1.0,0.0,0.0))
    rd = Dual(o.r[1],     (0.0,0.0,1.0,0.0))
    td = Dual(-1e-30,     (0.0,0.0,0.0,1.0))
    zd = one(td)*o.z[1]
    gcp = GCParticle(ed,pd,rd,zd,c.m,c.q)
    path, stat = integrate(M, gcp; tmax=td*tau_p, interp_dt=0.0,
                           one_transit=false, store_path=true,
                           classify_orbit=false, drift=o.drift,vacuum=o.vacuum)

    x = transform([path.energy[end],path.pitch[end],path.r[end],path.z[end]])
    detJ[np] = max(abs(det(hcat((partials(xx) for xx in x)...))),0.0)

    ep = partials(path.energy[end])
    pp = partials(path.pitch[end])
    rp = partials(path.r[end])
    zp = partials(path.z[end])
    for i=(np-1):-1:(np1+1)
        ed = Dual(o.energy[i+1],ep)
        pd = Dual(o.pitch[i+1], pp)
        rd = Dual(o.r[i+1],     rp)
        zd = Dual(o.z[i+1],     zp)
        gcp = GCParticle(ed,pd,rd,zd,c.m,c.q)

        dt = o.dt[i]
        dt = -(dt + 1e-30*iszero(dt))
        path, stat = integrate(M, gcp; tmax=dt, interp_dt=0.0,
                               one_transit=false, store_path=true,
                               classify_orbit=false, drift=o.drift,vacuum=o.vacuum)

        x = transform([path.energy[end],path.pitch[end],path.r[end],path.z[end]])
        detJ[i] = max(abs(det(hcat((partials(xx) for xx in x)...))),0.0)

        ep = partials(path.energy[end])
        pp = partials(path.pitch[end])
        rp = partials(path.r[end])
        zp = partials(path.z[end])
    end

    if tol > 0.0 && !check_jacobian(detJ; tol=tol)
        verbose && @warn "Jacobian calculation failed. Setting to zero."
        detJ = zero(detJ)
    end

    return detJ
end

function transform(M::AbstractEquilibrium, o::Orbit, ::Type{EPRParticle}; kwargs...)
    J = _get_jacobian(M, o.coordinate, o.path, o.tau_p; kwargs...)[1]
    return EPRParticle(o.coordinate), J
end

function transform(M::AbstractEquilibrium, gcp::GCParticle, ::Type{EPRParticle}; kwargs...)
    o = get_orbit(M, gcp; kwargs...)
    return transform(M, o; kwargs...)
end

function transform(M::AbstractEquilibrium, c::EPRCoordinate, ::Type{HamiltonianCoordinate})
    KE_d = Dual(c.energy, (1.0, 0.0, 0.0))
    p_d =  Dual(c.pitch,  (0.0, 1.0, 0.0))
    r_d =  Dual(c.r,      (0.0, 0.0, 1.0))
    c_d = EPRCoordinate(KE_d, p_d, r_d, one(r_d)*c.z, one(r_d)*c.t, c.m, c.q)

    hc_d = HamiltonianCoordinate(M, c_d)

    hc = HamiltonianCoordinate(value(hc_d.energy), value(hc_d.mu), value(hc_d.p_phi), c.m, c.q)
    J = abs(det(hcat(partials(hc_d.energy), partials(hc_d.mu), partials(hc_d.p_phi))))

    return hc, J
end
