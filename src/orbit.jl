abstract AbstractOrbitCoordinate{T}

immutable EPRZCoordinate{T} <: AbstractOrbitCoordinate{T}
    energy::T
    pitch::T
    r::T
    z::T
end

immutable CQL3DCoordinate{T} <: AbstractOrbitCoordinate{T}
    energy::T
    pitch_w::T
    psi_w::T
    r_w::T
    z_w::T
end

function CQL3DCoordinate(energy, pitch, R, M::AxisymmetricEquilibrium)
    psi = M.psi([R,M.axis[2]])
    rmin = R-0.01
    rmax = R+0.01
    r = linspace(rmin,rmax,1000)
    zmin = M.axis[2]-0.01
    zmax = M.axis[2]+0.01
    z = linspace(zmin,zmax,1000)
    psirz = [M.psi([rr,zz]) for rr in r, zz in z]
    l = Contour.contour(r,z,psirz,psi)
    rc, zc = coordinates(lines(l)[1])
    i = indmax(rc)
    r_w = rc[i]
    z_w = zc[i]
    return CQL3DCoordinate(energy, pitch, psi, r_w, z_w)
end

immutable HamiltonianCoordinate{T} <: AbstractOrbitCoordinate{T}
    energy::T
    mu::T
    p_phi::T
end

function HamiltonianCoordinate(c::CQL3DCoordinate, M::AxisymmetricEquilibrium; amu=H2_amu, Z=1.0)
    psi = c.psi_w
    babs = M.b([c.r_w,c.z_w])
    g = M.g([psi])

    E = c.energy
    mu = e0*1e3*E*(1-c.pitch_w^2)/babs
    Pphi = -M.sigma*sqrt(2e3*e0*E*mass_u*amu)*g*c.pitch_w/babs + Z*e0*psi
    hc = HamiltonianCoordinate(E,mu,Pphi)
    return hc
end

function HamiltonianCoordinate(c::HamiltonianCoordinate, M::AxisymmetricEquilibrium; amu=H2_amu, Z=1.0)
    return c
end

function HamiltonianCoordinate(c::EPRZCoordinate, M::AxisymmetricEquilibrium; amu=H2_amu, Z=1.0)
    psi = M.psi([c.r,c.z])
    babs = M.b([c.r,c.z])
    g = M.g([psi])

    E = c.energy
    mu = e0*1e3*E*(1-c.pitch^2)/babs
    Pphi = -M.sigma*sqrt(2e3*e0*E*mass_u*amu)*g*c.pitch/babs + Z*e0*psi
    hc = HamiltonianCoordinate(E,mu,Pphi)
    return hc
end

function CQL3DCoordinate(c::EPRZCoordinate, M::AxisymmetricEquilibrium; amu=H2_amu, Z=1.0)

    hamilc = HamiltonianCoordinate(c, M, amu=amu, Z=Z)

    #Function to check if in function
    hits_boundary = [false]
    in_boundary = function ib(x)
        if (M.r_domain[1] < x[1] < M.r_domain[2]) && (M.z_domain[1] < x[3] < M.z_domain[2])
            return true
        else
            hits_boundary[1] = true
            return false
        end
    end

    #Function to check if orbit is complete
    r0 = [c.r, 0.0, c.z]
    r1 = copy(r0)
    initial_dir = [0]
    npol = [0]
    orbit_complete = [false]
    is_complete = function ic(r2)
        if initial_dir[1] == 0
            initial_dir[1] = sign(r2[3]-r0[3])
            r1[:] = convert(Vector, r2)
            return false
        end

        dr10 = r1 - r0
        dr20 = r2 - r0
        dr21 = r2 - r1

        if sign(dr10[3]*dr20[3]) < 0
            npol[1] = npol[1] + 1
            if norm(dr10[1:2:3]) < norm(dr21[1:2:3]) && sign(dr20[3]) == initial_dir[1]
                orbit_complete[1] = npol[1] == 2 || npol[1] == 4
                return true
            end
        end
        r1[:] = convert(Vector, r2)
        return false
    end

    #Function to call between timesteps
    cb = function callback(x)
        !is_complete(x) && in_boundary(x)
    end

    if !in_boundary(r0)
        error("Starting point outside boundary")
    end

    f = make_gc_ode(M, c, amu=amu, Z=Z)
    t = 1e-6*collect(linspace(0.0,600.0,24000))

    res = Sundials.cvode(f, r0, t, reltol=1e-8,abstol=1e-12, callback = cb)
    if !orbit_complete[1] && !hits_boundary[1]
        warn("EPRZCoordinate cannot be expressed as a CQL3DCoordinate: ",npol[1]," ",hits_boundary)
        return (Nullable{CQL3DCoordinate}(),:loss)
    end

    r = res[:,1]
    phi = res[:,2]
    z = res[:,3]

    pitch = [get_pitch(c, M, rr, zz, amu=amu, Z=Z) for (rr,zz) in zip(r,z)]
    class = classify_orbit(r,z,pitch,M.axis)

    rind = indmax(r)
    r_w = r[rind]
    z_w = z[rind]
    pitch_w = pitch[rind]
    psi_w = M.psi([r_w,z_w])
    return (Nullable{CQL3DCoordinate}(CQL3DCoordinate(c.energy, pitch_w, psi_w, r_w, z_w)), class)
end

type Orbit{T,S<:AbstractOrbitCoordinate{Float64},R<:AbstractOrbitCoordinate{Float64}}
    hcoord::S
    ccoord::R
    r::Vector{T}
    z::Vector{T}
    phi::Vector{T}
    dt::Vector{T}
    pitch::Vector{T}
    class::Symbol
end

function get_pitch(c::HamiltonianCoordinate, M::AxisymmetricEquilibrium, r, z; amu = H2_amu, Z=1.0)
    psi = M.psi([r,z])
    g = M.g([psi])
    babs = M.b([r,z])

    f = -babs/(sqrt(2e3*e0*c.energy*mass_u*amu)*g*M.sigma)
    pitch = f*(c.p_phi - Z*e0*psi)
    pitchabs = sqrt(max(1.0-(c.mu*babs/(1e3*e0*c.energy)), 0.0))
    if !isapprox(abs(pitch), pitchabs, atol=1.e-1)
        warn("abs(pitch) != abspitch: ",pitchabs," ",pitch)
    end
    return clamp(pitch,-1.0,1.0)
end

function get_pitch{T<:AbstractOrbitCoordinate}(c::T, M::AxisymmetricEquilibrium, r, z; amu=H2_amu, Z=1.0)
    hcoord = HamiltonianCoordinate(c, M, amu=amu, Z=Z)
    return get_pitch(hcoord, M, r, z, amu=amu, Z=Z)
end

function classify_orbit{T}(r::Vector{T},z::Vector{T},pitch::Vector{T},axis)
    op = Polygon()
    for p in zip(r,z)
        push!(op.vertices,p)
    end
    push!(op.vertices,(r[1],z[1]))

    if in_polygon(axis, op)
        if all(sign.(pitch) .== sign(pitch[1]))
            if sign(pitch[1]) > 0.0
                class = Symbol("co_passing")
            else
                class = Symbol("ctr_passing")
            end
        else
            class = Symbol("potato")
        end
    else
        if all(sign.(pitch) .== sign(pitch[1]))
            class = Symbol("stagnation")
        else
            class = Symbol("trapped")
        end
    end

    return class
end

function make_gc_ode{T<:AbstractOrbitCoordinate}(M::AxisymmetricEquilibrium, c::T; amu=H2_amu, Z=1.0)
    oc = HamiltonianCoordinate(c, M, amu = amu, Z=Z)
    res = DiffBase.GradientResult(rand(2))
    function vgc(t, y, ydot)
        ri = [y[1],y[3]]
        psi = M.psi(ri)
        g = M.g(psi)
        B = M.B(ri)
        ForwardDiff.gradient!(res, M.b, ri)
        gradB = DiffBase.gradient(res)
        babs = DiffBase.value(res)

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

function calc_orbit(M::AxisymmetricEquilibrium, wall::Polygon, c::CQL3DCoordinate; amu=H2_amu, Z=1.0, nstep=3000, tmax=500.0)

    hamilc = HamiltonianCoordinate(c, M, amu=amu, Z=Z)

    #Function to check if in function
    hits_boundary = [false]
    in_boundary = function ib(x)
        if x[1] - c.r_w <= 0.0 && in_polygon(x[1:2:3],wall)
            return true
        else
            hits_boundary[1] = true
            return false
        end
    end

    #Function to check if orbit is complete
    r0 = [c.r_w, 0.0, c.z_w]
    r1 = copy(r0)
    initial_dir = [0]
    npol = [0]
    orbit_complete = [false]
    is_complete = function ic(r2)
        if initial_dir[1] == 0
            initial_dir[1] = sign(r2[3]-r0[3])
            r1[:] = convert(Vector, r2)
            return false
        end

        dr10 = r1 - r0
        dr20 = r2 - r0
        dr21 = r2 - r1

        if sign(dr10[3]*dr20[3]) < 0
            npol[1] = npol[1] + 1
            if norm(dr10[1:2:3]) < norm(dr21[1:2:3]) && sign(dr20[3]) == initial_dir[1]
                orbit_complete[1] = npol[1] == 2 || npol[1] == 4
                return true
            end
        end
        r1[:] = convert(Vector, r2)
        return false
    end

    #Function to call between timesteps
    cb = function callback(mem, t, x)
        !is_complete(x) && in_boundary(x)
    end

    if !in_boundary(r0)
        error("Starting point outside wall")
    end

    f = make_gc_ode(M, c, amu=amu, Z=Z)
    t = 1e-6*collect(linspace(0.0,tmax,nstep))

    res = Sundials.cvode(f, r0, t, reltol=1e-8,abstol=1e-12, callback = cb)

    r = res[:,1]
    phi = res[:,2]
    z = res[:,3]
    n = length(r)

    pitch = zeros(length(r))
    for i=1:n
        pitch[i] = get_pitch(hamilc,M,r[i],z[i])
    end

    if hits_boundary[1] || !orbit_complete[1]
        class = Symbol("loss")
    else
        class = classify_orbit(r,z,pitch,M.axis)
    end

    dt = fill(t[2]-t[1],n)

    #correct for last timestep
    ydot = zeros(3)
    flag = f(0.0, r0, ydot)
    dt[end] = sqrt((r[end] - r[1])^2 + (z[end] - z[1])^2)/norm(ydot[1:2:3])

    return Orbit(hamilc, c, r, z, phi, dt, pitch, class)
end

function create_energy(M::AxisymmetricEquilibrium, c::CQL3DCoordinate, amu=H2_amu, Z=1.0)
    oc = HamiltonianCoordinate(c,M,amu,Z)
    function energy_function(x)
        psi = M.psi(x)
        g = M.g(psi)
        babs = M.b(x)
        E = ((((oc.p_phi - Z*e0*psi)^2)/(2*mass_u*amu))*(babs/g)^2 + oc.mu*babs)/(e0*1e3)
        return E
    end
end

function create_mu(M::AxisymmetricEquilibrium, c::CQL3DCoordinate, amu=H2_amu, Z=1.0)
    oc = HamiltonianCoordinate(c,M,amu=amu,Z=Z)
    function mu_function(x)
        psi = M.psi(x)
        g = M.g(psi)
        babs = M.b(x)
        mu = (1e3*e0*oc.energy/babs) - (babs/(2*mass_u*amu))*((oc.p_phi - Z*e0*psi)/g)^2
        return mu
    end
end

function gc_velocity(M::AxisymmetricEquilibrium,c::HamiltonianCoordinate,ri; amu=H2_amu, Z=1.0)
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

function plot_orbit(o::Orbit;rlim=(1.0,2.5),zlim=(-1.25,1.25),xlim = (-2.3,2.3),ylim=(-2.3,2.3))
#    pyplot()
    t = cumsum(o.dt)*1e6
    n = length(t)
    r = o.r
    z = o.z
    x = r.*cos.(o.phi)
    y = r.*sin.(o.phi)
    l = @layout [a{0.38w} b]
    p1 = plot(layout=l)

    rr = vec(hcat(r[1:end-1],r[2:end],fill(NaN,n-1))')
    zz = vec(hcat(z[1:end-1],z[2:end],fill(NaN,n-1))')
    plot!(p1,rr,zz,line_z=t,c=:bluesreds,label=o.class,
          colorbar=false,subplot=1,xlim=rlim,
          ylim=zlim,xlabel="R [m]",ylabel="Z [m]",aspect_ratio=1.0)

    xx = vec(hcat(x[1:end-1],x[2:end],fill(NaN,n-1))')
    yy = vec(hcat(y[1:end-1],y[2:end],fill(NaN,n-1))')
    plot!(p1,xx,yy,line_z=t,c=:bluesreds,label=o.class,subplot=2,
          xlim=xlim,ylim=ylim,xlabel="X [m]",
          ylabel="Y [m]",aspect_ratio=1.0)

    return p1
end
