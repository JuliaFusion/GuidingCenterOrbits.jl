function eprz_to_eprt(M, energy, pitch, r, z; m=H2_amu, q=1, kwargs...)

    f = function c(x; kwargs...)
        gcp = GCParticle(x[1], x[2], x[3], x[4], H2_amu*mass_u, q)
        path, stat = integrate(M, gcp; store_path=false, one_transit=true, kwargs...)
        if stat.class in (:lost,:incomplete)
            return x
        end
        hc = HamiltonianCoordinate(M, gcp)
        KE = get_kinetic_energy(M, hc, stat.rm, stat.zm)
        return [KE, stat.pm, stat.rm, (stat.tau_p - stat.tm)/stat.tau_p]
    end

    x = [energy,pitch,r,z]

    jr = DiffResults.JacobianResult(x)
    jcfg = ForwardDiff.JacobianConfig(nothing, x)
    ForwardDiff.jacobian!(jr, x->f(x; kwargs...), x, jcfg, Val{false}())

    v = DiffResults.value(jr)
    detJ = abs(det(DiffResults.jacobian(jr)))

    return v, detJ
end
