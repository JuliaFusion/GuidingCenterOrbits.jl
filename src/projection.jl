const OrbitProjection = Equilibrium.Boundary

const in_orbit = in_boundary

function orbit_projection(M::AbstractEquilibrium, gcp; dx = 0.01, dy =0.01)
    hc = HamiltonianCoordinate(M,gcp)
    rlims, zlims = limits(M)

    x = range(rlims...,step=dx)
    y = range(zlims...,step=dy)

    mu = zeros(length(x), length(y))
    for ix in eachindex(x), iy in eachindex(y)
        B = norm(Bfield(M, x[ix], y[iy]))
        ir, iphi, iz = cylindrical_cocos_indices(cocos(M))
        g = poloidal_current(M, M(x[ix], y[iy])) # RBphi
        energy = 1e3*hc.energy*e0
        mu[ix,iy] = max(energy/B - (B/(2*hc.m))*((hc.p_phi - hc.q*e0*M(x[ix],y[iy]))/g)^2,0.0)
    end

    cl = Contour.contour(x,y,mu,hc.mu)
    ls = Contour.lines(cl)
    oproj = OrbitProjection([(gcp.r, gcp.z)])
    if length(ls) == 0
        return oproj
    end
    omin = Inf
    for (i,l) in enumerate(ls)
        r, z = Contour.coordinates(l)
        omin_tmp = minimum((r .- gcp.r).^2 .+ (z .- gcp.z).^2)
        if omin_tmp < omin
            omin = omin_tmp
            oproj = OrbitProjection(collect(zip(r,z)))
        end
    end

    return oproj
end
