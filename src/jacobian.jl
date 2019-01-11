function eprz_to_eprt(M, energy, pitch, r, z; m=H2_amu, q=1, auto_diff = true, kwargs...)

    islost=[false]
    erred = [true]
    f = function c(x; kwargs...)
        erred[1] = false
        islost[1] = false
        gcp = GCParticle(x[1], x[2], x[3], x[4], H2_amu*mass_u, q)
        path, stat = integrate(M, gcp; store_path=false, classify_orbit=false, one_transit=true, kwargs...)
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
