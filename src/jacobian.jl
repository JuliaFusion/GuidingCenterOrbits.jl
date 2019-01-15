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

function get_orbit_jacobian(M::AxisymmetricEquilibrium, c::EPRCoordinate; kwargs...)
    ed = ForwardDiff.Dual(c.energy,(1.0,0.0,0.0,0.0))
    pd = ForwardDiff.Dual(c.pitch, (0.0,1.0,0.0,0.0))
    rd = ForwardDiff.Dual(c.r,     (0.0,0.0,1.0,0.0))
    td = ForwardDiff.Dual(0.0,     (0.0,0.0,0.0,1.0))
    zd = one(td)*c.z
    cdual = EPRCoordinate(ed,pd,rd,zd,c.m,c.q)
    path, stat = integrate(M, cdual; tmin = td, one_transit = true, kwargs...)

    np = length(path)
    detJ = zeros(np)
    for i=1:np
        detJ[i] = abs(det(hcat((ForwardDiff.partials(p) for p in (path.energy[i],path.pitch[i],path.r[i],path.z[i]))...)))
    end

    r = ForwardDiff.value.(path.r)
    z = ForwardDiff.value.(path.z)
    phi = ForwardDiff.value.(path.phi)
    energy = ForwardDiff.value.(path.energy)
    pitch = ForwardDiff.value.(path.pitch)
    dt = ForwardDiff.value.(path.dt)
    tau_p = ForwardDiff.value(stat.tau_p)
    tau_t = ForwardDiff.value(stat.tau_t)

    rmax = maximum(r)
    if stat.class != :incomplete && stat.class != :lost
        rmax > c.r && !isapprox(rmax,c.r,rtol=1e-4) && (stat.class = :degenerate)
    else
        tau_p=zero(tau_p)
        tau_t=zero(tau_t)
    end
    path = OrbitPath(r,z,phi,pitch,energy,dt)
    return Orbit(c,stat.class,tau_p,tau_t,path), detJ
end
