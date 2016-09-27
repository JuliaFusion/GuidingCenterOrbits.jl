abstract AbstractOrbitCoordinate{T}

immutable HamiltonianCoordinate{T} <: AbstractOrbitCoordinate{T}
    energy::T
    mu::T
    p_phi::T
end

immutable RomeCoordinate{T} <: AbstractOrbitCoordinate{T}
    v::T
    pitch::T
    psi::T
end

immutable CQL3DCoordinate{T} <: AbstractOrbitCoordinate{T}
    energy::T
    pitch::T
    R::T
end

type Orbit{T,S<:AbstractOrbitCoordinate{Float64},R<:AbstractOrbitCoordinate{Float64}}
    hcoord::S
    ccoord::R
    r::Vector{T}
    z::Vector{T}
    phi::Vector{T}
    t::Vector{T}
    hits_boundary::Bool
    complete::Bool
end

function get_pitch(c::HamiltonianCoordinate, M::AxisymmetricEquilibrium, r, z; amu = 2.0141, Z=1.0)
    psi = M.psi([r,z])
    g = M.g([psi])
    babs = M.b([r,z])

    pitch = -babs*(c.p_phi - Z*e0*psi)/(sqrt(2e3*e0*c.energy*mass_u*amu)*g*M.sigma)
    return pitch
end

function cql3d_to_hamiltonian(c::CQL3DCoordinate, M::AxisymmetricEquilibrium, amu=2.0141, Z=1.0)
    psi = M.psi([c.R,M.axis[2]])
    babs = M.b([c.R,M.axis[2]])
    g = M.g([psi])

    E = c.energy
    mu = e0*1e3*E*(1-c.pitch^2)/babs
    Pphi = -M.sigma*sqrt(2e3*e0*E*mass_u*amu)*g*c.pitch/babs + Z*e0*psi
    hc = HamiltonianCoordinate(E,mu,Pphi)
    return hc
end

function create_energy(M::AxisymmetricEquilibrium, c::CQL3DCoordinate, amu=2.0141, Z=1.0)
    oc = cql3d_to_hamiltonian(c,M,amu,Z)
    function energy_function(x)
        psi = M.psi(x)
        g = M.g(psi)
        babs = M.b(x)
        E = ((((oc.p_phi - Z*e0*psi)^2)/(2*mass_u*amu))*(babs/g)^2 + oc.mu*babs)/(e0*1e3)
        return E
    end
end

function create_mu(M::AxisymmetricEquilibrium, c::CQL3DCoordinate, amu=2.0141, Z=1.0)
    oc = cql3d_to_hamiltonian(c,M,amu,Z)
    function mu_function(x)
        psi = M.psi(x)
        g = M.g(psi)
        babs = M.b(x)
        mu = (1e3*e0*oc.energy/babs) - (babs/(2*mass_u*amu))*((oc.p_phi - Z*e0*psi)/g)^2
        return mu
    end
end

function gc_velocity(M::AxisymmetricEquilibrium,c::HamiltonianCoordinate,ri; amu=2.0141, Z=1.0)
    E = c.energy
    pphi = c.p_phi
    mu = c.mu

    psi = M.psi(ri)
    g = M.g(psi)
    B = M.B(ri)
    gradB = ForwardDiff.gradient(M.b,ri)
    babs = norm(B)

    Wperp = mu*babs
    Wpara = 1e3*e0*E - Wperp
    vpara = -babs*(pphi - Z*e0*psi)/(amu*mass_u*g)

    vd = (1/(Z*e0))*(Wperp + 2*Wpara)*cross(B,[gradB[1],0.0,gradB[2]])/(babs^3)
    vg = vpara*B/babs + vd
    return vg
end

function make_gc_ode(M::AxisymmetricEquilibrium, c::CQL3DCoordinate, amu=2.0141, Z=1.0)
    oc = cql3d_to_hamiltonian(c, M, amu, Z)
    res = ForwardDiff.GradientResult(rand(2))
    function vgc(t, y, ydot)
        ri = [y[1],y[3]]
        psi = M.psi(ri)
        g = M.g(psi)
        B = M.B(ri)
        ForwardDiff.gradient!(res, M.b, ri)
        gradB = ForwardDiff.gradient(res)
        babs = ForwardDiff.value(res)

        Wperp = oc.mu*babs
        Wpara = 1e3*e0*oc.energy - Wperp
        vpara = -babs*(oc.p_phi - Z*e0*psi)/(amu*mass_u*g)

        vd = (1/(Z*e0))*(Wperp + 2*Wpara)*cross(B,[gradB[1],0.0,gradB[2]])/(babs^3)
        vg = vpara*B/babs + vd
        ydot[1] = vg[1]
        ydot[2] = vg[2]/ri[1]
        ydot[3] = vg[3]
        return Sundials.CV_SUCCESS
    end
end

function calc_orbit(M::AxisymmetricEquilibrium, wall::Polygon, c::CQL3DCoordinate; amu=2.0141, Z=1.0, nstep=3000, tmax=500.0)

    hamilc = cql3d_to_hamiltonian(c, M, amu, Z)

    bdry = function ib(x)
        x[1] - c.R <= 1e-3 &&
        in_polygon(x,wall) #&&
        #(M.r_domain[1] < x[1] < M.r_domain[2]) &&
        #(M.z_domain[1] < x[2] < M.z_domain[2])
    end

    if !bdry([c.R,M.axis[2]])
        error("Starting point outside wall")
    end

    f = make_gc_ode(M, c, amu, Z)
    y0 = [c.R, 0.0, M.axis[2]]
    t = 1e-6*collect(linspace(0.0,tmax,nstep))

    res = Sundials.cvode(f, y0, t)

    r = res[:,1]
    phi = res[:,2]
    z = res[:,3]

    n = length(r)

    r0 = [c.R,M.axis[2]]
    cnt = 1
    ncross = 0
    hits_boundary = false
    complete = false
    @inbounds for i=2:n-1
        cnt = cnt+1
        r1 = [r[i],z[i]]
        r2 = [r[i+1],z[i+1]]

        if !bdry(r2)
            hits_boundary = true
            break
        end

        dr10 = r1 - r0
        dr20 = r2 - r0
        dr21 = r2 - r1

        if sign(dr10[2]*dr20[2]) < 0
            ncross = ncross + 1
            if abs(dr10[1]) < 0.01 && sign(dr20[2]) == sign((z[2]-z[1]))
                complete = true
                break
            end
        end
    end
    if ncross > 2
        complete = false
    end

    if !complete && !hits_boundary
        warn("Orbit does not make one polodial transit: ",c)
        cnt = cnt + 1
    end

    rp = r[1:cnt]
    zp = z[1:cnt]
    phip = phi[1:cnt]
    tp = t[1:cnt]

    return Orbit(hamilc, c, rp, zp, phip, tp, hits_boundary, complete)
end

function calc_orbit_contour(M::AxisymmetricEquilibrium, wall::Polygon, c::CQL3DCoordinate;
                    amu=2.0141, Z=1.0, s_range=(0.01,0.01), tol=1e-8,
                    step_range = (1e-6,0.01), max_step=3000,verbose=true)

    raxis, zaxis = M.axis
    c.R < raxis && error("r < raxis")

    energy = create_energy(M, c, amu, Z)

    if !isapprox(energy([c.R,zaxis]), energy.oc.energy,atol=1e-16)
        warn(energy([c.R,zaxis]), " !≈ ",energy.oc.energy)
        return Nullable{Orbit}()
    end

    bdry = function ib(x)
        x[1] - c.R <= 1e-3 &&
        in_polygon(x,wall) #&&
        #(M.r_domain[1] < x[1] < M.r_domain[2]) &&
        #(M.z_domain[1] < x[2] < M.z_domain[2])
    end

    p = follow_contour(energy, [c.R, zaxis], tol=tol,
        step_range=step_range, boundary = bdry,
        maxiter=max_step,verbose=verbose)

    if isnull(p)
        r = Float64[]
        z = Float64[]
        phi = Float64[]
        t = Float64[]
        orb = Orbit(c,c,r,z,phi,t,true,false)
        return orb
    end

    path = getfield(get(p),:vertices)
    if length(path) == 1
        println(c)
    end

    dr = [path[2][1] - path[1][1],path[2][2] - path[1][2]]
    vgc = gc_velocity(M, energy.oc, collect(path[1]))
    if (dr[1]*vgc[1] + dr[2]*vgc[3]) < 0.0
        reverse!(path)
    end

    np = length(path)
    r = [path[i][1] for i=1:np]
    z = [path[i][2] for i=1:np]

    if (std(r) < s_range[1] && std(z) < s_range[2])
        # Stagnation orbit
        r = [path[1][1],path[1][1]]
        z = [path[1][2],path[1][2]]
        vt = gc_velocity(M, energy.oc, collect(path[1]))[2]
        t = [2pi*r[1]/vt, 0.0]
        phi = [0.0, 2pi]
        orb = Orbit(c, c, r, z, phi, t, false,true)
        return orb
    end

    phi = zeros(np)
    dt = zeros(np)
    vgc = zeros(3,np)
    vgc[:,1] = gc_velocity(M, energy.oc, collect(path[1]))
    @inbounds for i=2:np
        ri = collect(path[i])
        vgc[:,i] = gc_velocity(M, energy.oc, ri)
        d = sqrt((r[i] - r[i-1])^2 + (z[i] - z[i-1])^2)

        vp1 = sqrt(vgc[1,i-1]^2 + vgc[3,i-1]^2)
        vp2 = sqrt(vgc[1,i]^2 + vgc[3,i]^2)
        vt1 = vgc[2,i-1]
        vt2 = vgc[2,i]
        dt[i-1] = d/(0.5*(vp1+vp2))

        phi[i] = phi[i-1] + atan2((dt[i-1])*(vt1 + vt2)/2, r[i-1])
    end

    orb = Orbit(c, c, r, z, phi, cumsum(dt), false, true)
    return orb
end

function plot_orbit(o::Orbit;rlim=(1.0,2.5),zlim=(-1.25,1.25),xlim = (-2.3,2.3),ylim=(-2.3,2.3))
#    pyplot()
    t = o.t*1e6
    n = length(t)
    r = o.r
    z = o.z
    x = r.*cos.(o.phi)
    y = r.*sin.(o.phi)
    l = @layout [a{0.38w} b]
    p1 = plot(layout=l)

    rr = vec(hcat(r[1:end-1],r[2:end],fill(NaN,n-1))')
    zz = vec(hcat(z[1:end-1],z[2:end],fill(NaN,n-1))')
    plot!(p1,rr,zz,line_z=t,c=:bluesreds,label="",
          colorbar=false,subplot=1,xlim=rlim,
          ylim=zlim,xlabel="R [m]",ylabel="Z [m]",aspect_ratio=1.0)

    xx = vec(hcat(x[1:end-1],x[2:end],fill(NaN,n-1))')
    yy = vec(hcat(y[1:end-1],y[2:end],fill(NaN,n-1))')
    plot!(p1,xx,yy,line_z=t,c=:bluesreds,label="",subplot=2,
          xlim=xlim,ylim=ylim,xlabel="X [m]",
          ylabel="Y [m]",aspect_ratio=1.0)

    return p1
end
