function eprz_to_eprt(M, energy, pitch, r, z; m=H2_amu, q=1, auto_diff = true, kwargs...)

    islost=[false]
    erred = [true]
    f = function c(x; kwargs...)
        erred[1] = false
        islost[1] = false
        gcp = GCParticle(x[1], x[2], x[3], x[4], H2_amu*mass_u, q)
        path, stat = integrate(M, gcp; autodiff=false,store_path=false, classify_orbit=false, one_transit=true, kwargs...)
        if stat.class == :lost
            islost[1] = true
            return x
        elseif stat.errcode == 1
            erred[1] = true
            return x
        end
        hc = HamiltonianCoordinate(M, gcp)
        KE = get_kinetic_energy(M, hc, stat.rm, stat.zm)
        return [KE, stat.pm, stat.rm, (stat.tau_p - stat.tm)/stat.tau_p]
    end

    x = [energy,pitch,r,z]

    if auto_diff
        jr = DiffResults.JacobianResult(x)
        jcfg = ForwardDiff.JacobianConfig(nothing, x)
        ForwardDiff.jacobian!(jr, x->f(x; kwargs...), x, jcfg, Val{false}())

        v = DiffResults.value(jr)
        if !erred[1]
            if islost[1]
                detJ = 0.0
            else
                detJ = abs(det(DiffResults.jacobian(jr)))
            end
        else
            detJ = 0.0
        end
    end
    if erred[1] && !auto_diff
        v = f(x; kwargs...)
        if islost[1] || erred[1]
            detJ = 0.0
        else
            J = DiffEqDiffTools.finite_difference_jacobian(x->f(x; adaptive=false, kwargs...), x,
                                                           Val{:central}, eltype(x), Val{false})
            detJ = abs(det(J))
        end
    end
    return v, detJ
end

function eprt_to_eprz(M, energy, pitch, r, t; tau_p, m=H2_amu, q=1, auto_diff = true, kwargs...)
    islost=[false]
    erred = [true]
    f = function c(x; kwargs...)
        erred[1] = false
        islost[1] = false
        gcp = EPRCoordinate(M,x[1], x[2], x[3], amu=m, q=q)
        t = (x[4] + 1e-30*iszero(x[4]))
        if t > 0.5
            t = -(1.0 - t)
        end
        path, stat = integrate(M, gcp; tmax = t*tau_p, interp_dt=0, store_path=true,
                               classify_orbit=false, one_transit=false, integrator=BS3(), kwargs...)
        if stat.class == :lost
            islost[1] = true
            return x
        elseif stat.errcode == 1
            erred[1] = true
            return x
        end
        return [path.energy[end], path.pitch[end], path.r[end], path.z[end]]
    end

    x = [energy,pitch,r,t]

    if auto_diff
        jr = DiffResults.JacobianResult(x)
        jcfg = ForwardDiff.JacobianConfig(nothing, x)
        ForwardDiff.jacobian!(jr, x->f(x; kwargs...), x, jcfg, Val{false}())

        v = DiffResults.value(jr)
        if !erred[1]
            if islost[1]
                detJ = 0.0
            else
                detJ = abs(det(DiffResults.jacobian(jr)))
            end
        else
            detJ = 0.0
        end
    end
    if erred[1] && !auto_diff
        v = f(x; kwargs...)
        if islost[1] || erred[1]
            detJ = 0.0
        else
            J = DiffEqDiffTools.finite_difference_jacobian(x->f(x;kwargs...), x,
                                                           Val{:central}, eltype(x), Val{false})
            detJ = abs(det(J))
        end
    end
    return v, detJ
end

function _get_shifted_jacobian(M, o::Orbit; kwargs...)
    # Create spline of orbit path
    t = range(0.0,stop=1.0,length=length(o.path))
    A = hcat(o.path.energy,o.path.pitch,o.path.r,o.path.z)
    itp = extrapolate(scale(interpolate(A,(BSpline(Cubic(Periodic(OnGrid()))), NoInterp())), t, 1:4), (Periodic(), Throw()))

    # Find maximum r
    res = optimize(x->-itp(x,3), 0.0, 1.0)
    tm = Optim.minimizer(res)

    # Create new orbit path
    energy = itp.(t .+ tm, 1)
    pitch = itp.(t .+ tm, 2)
    r = itp.(t .+ tm, 3)
    z = itp.(t .+ tm, 4)

    dt = fill(step(t)*o.tau_p,length(o.path))
    dt[end] = 0.0
    opath = OrbitPath(r,z,r*0.0,pitch,energy,dt)

    # Calculate jacobian
    J = _get_jacobian(M, opath, o.tau_p; kwargs...)

    # Shift jacobian
    Jitp = extrapolate(scale(interpolate(J,BSpline(Cubic(Periodic(OnGrid())))),t), Periodic())
    Jshifted = Jitp.(t .- tm)

    return Jshifted
end

function get_jacobian(M::AxisymmetricEquilibrium, c::EPRCoordinate; kwargs...)
    o = get_orbit(M,c; kwargs...)
    if o.class == :degenerate
        return _get_shifted_jacobian(M, o; kwargs...)
    end
    return _get_jacobian(M, o.path, o.tau_p; kwargs...)
end

function get_jacobian(M::AxisymmetricEquilibrium, o::Orbit; kwargs...)
    r0 = o.path.r[1]
    if r0 < o.coordinate.r && !isapprox(r0,o.coordinate.r,rtol=1e-4)
        return _get_shifted_jacobian(M, o; kwargs...)
    end
    return _get_jacobian(M, o.path, o.tau_p; kwargs...)
end

function _get_jacobian(M::AxisymmetricEquilibrium, o::OrbitPath, tau_p; kwargs...)
    np = length(o)
    detJ = zeros(np)
    if np == 0 || tau_p == 0.0
        return detJ
    end

    ## Do first half of orbit
    ed = ForwardDiff.Dual(o.energy[1],(1.0,0.0,0.0,0.0))
    pd = ForwardDiff.Dual(o.pitch[1], (0.0,1.0,0.0,0.0))
    rd = ForwardDiff.Dual(o.r[1],     (0.0,0.0,1.0,0.0))
    td = ForwardDiff.Dual(1e-30,      (0.0,0.0,0.0,1.0))
    zd = one(td)*o.z[1]
    gcp = GCParticle(ed,pd,rd,zd)
    path, stat = integrate(M, gcp; tmax=td*tau_p*1e6, interp_dt=0.0,
                           one_transit=false, store_path=true,
                           classify_orbit=false, kwargs...)

    ep = ForwardDiff.partials(path.energy[end])
    pp = ForwardDiff.partials(path.pitch[end])
    rp = ForwardDiff.partials(path.r[end])
    zp = ForwardDiff.partials(path.z[end])
    detJ[1] = abs(det(hcat(ep,pp,rp,zp)))
    np1 = floor(Int,np/2)
    for i=2:np1
        ed = ForwardDiff.Dual(o.energy[i-1],ep)
        pd = ForwardDiff.Dual(o.pitch[i-1], pp)
        rd = ForwardDiff.Dual(o.r[i-1],     rp)
        zd = ForwardDiff.Dual(o.z[i-1],     zp)
        gcp = GCParticle(ed,pd,rd,zd)

        dt = o.dt[i-1]
        dt = dt + 1e-30*iszero(dt)
        path, stat = integrate(M, gcp; tmax=dt*1e6, interp_dt=0.0,
                               one_transit=false, store_path=true,
                               classify_orbit=false, kwargs...)

        ep = ForwardDiff.partials(path.energy[end])
        pp = ForwardDiff.partials(path.pitch[end])
        rp = ForwardDiff.partials(path.r[end])
        zp = ForwardDiff.partials(path.z[end])
        detJ[i] = abs(det(hcat(ep,pp,rp,zp)))
    end

    ## Do second half of orbit
    ed = ForwardDiff.Dual(o.energy[1],(1.0,0.0,0.0,0.0))
    pd = ForwardDiff.Dual(o.pitch[1], (0.0,1.0,0.0,0.0))
    rd = ForwardDiff.Dual(o.r[1],     (0.0,0.0,1.0,0.0))
    td = ForwardDiff.Dual(-1e-30,     (0.0,0.0,0.0,1.0))
    zd = one(td)*o.z[1]
    gcp = GCParticle(ed,pd,rd,zd)
    path, stat = integrate(M, gcp; tmax=td*tau_p*1e6, interp_dt=0.0,
                           one_transit=false, store_path=true,
                           classify_orbit=false, kwargs...)

    ep = ForwardDiff.partials(path.energy[end])
    pp = ForwardDiff.partials(path.pitch[end])
    rp = ForwardDiff.partials(path.r[end])
    zp = ForwardDiff.partials(path.z[end])
    detJ[np] = abs(det(hcat(ep,pp,rp,zp)))
    for i=(np-1):-1:(np1+1)
        ed = ForwardDiff.Dual(o.energy[i+1],ep)
        pd = ForwardDiff.Dual(o.pitch[i+1], pp)
        rd = ForwardDiff.Dual(o.r[i+1],     rp)
        zd = ForwardDiff.Dual(o.z[i+1],     zp)
        gcp = GCParticle(ed,pd,rd,zd)

        dt = o.dt[i]
        dt = -(dt + 1e-30*iszero(dt))
        path, stat = integrate(M, gcp; tmax=dt*1e6, interp_dt=0.0,
                               one_transit=false, store_path=true,
                               classify_orbit=false, kwargs...)

        ep = ForwardDiff.partials(path.energy[end])
        pp = ForwardDiff.partials(path.pitch[end])
        rp = ForwardDiff.partials(path.r[end])
        zp = ForwardDiff.partials(path.z[end])
        detJ[i] = abs(det(hcat(ep,pp,rp,zp)))
    end
    return detJ
end
