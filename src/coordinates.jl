abstract type AbstractOrbitCoordinate{T} end

struct EPRCoordinate{T} <: AbstractOrbitCoordinate{T}
    energy::T
    pitch::T
    r::T
    z::T
    amu::T
    q::Int
end

function EPRCoordinate(energy, pitch, r, z; amu=H2_amu, q = 1)
    return EPRCoordinate(energy, pitch, r, z, amu, q)
end

function EPRCoordinate(M::AxisymmetricEquilibrium, energy, pitch, R ; amu=H2_amu, q=1)
    psi = M.psi_rz[R,M.axis[2]]
    rmin = R-0.01
    rmax = R+0.01
    r = linspace(rmin,rmax,1000)
    zmin = M.axis[2]-0.01
    zmax = M.axis[2]+0.01
    z = linspace(zmin,zmax,1000)
    psirz = [M.psi_rz[rr,zz] for rr in r, zz in z]
    l = Contour.contour(r,z,psirz,psi)
    rc, zc = coordinates(lines(l)[1])
    i = indmax(rc)
    return EPRCoordinate(energy, pitch, rc[i], zc[i], amu, q)
end

function Base.show(io::IO, c::EPRCoordinate)
    println(io,typeof(c))
    @printf(io," E = %.3f keV\n",c.energy)
    @printf(io," pitch = %.3f\n",c.pitch)
    @printf(io," Rmax = %.3f m\n",c.r)
end

struct HamiltonianCoordinate{T} <: AbstractOrbitCoordinate{T}
    energy::T
    mu::T
    p_phi::T
    amu::T
    q::Int
end

function HamiltonianCoordinate(energy, mu, p_phi; amu=H2_amu, q=1)
    return HamiltonianCoordinate(energy, mu, p_phi, amu, q)
end

function HamiltonianCoordinate(M::AxisymmetricEquilibrium, c::EPRCoordinate)
    psi = M.psi_rz[c.r,c.z]
    babs = M.b[c.r,c.z]
    g = M.g[psi]

    E = c.energy
    mu = e0*1e3*E*(1-c.pitch^2)/babs
    Pphi = -M.sigma*sqrt(2e3*e0*E*mass_u*c.amu)*g*c.pitch/babs + c.q*e0*psi
    return HamiltonianCoordinate(E,mu,Pphi,c.amu,c.q)
end

function HamiltonianCoordinate(M::AxisymmetricEquilibrium, c::HamiltonianCoordinate)
    return c
end

function HamiltonianCoordinate(M::AxisymmetricEquilibrium, E, p, r, z; amu=H2_amu, q=1)
    psi = M.psi_rz[r,z]
    babs = M.b[r,z]
    g = M.g[psi]

    mu = e0*1e3*E*(1-p^2)/babs
    Pphi = -M.sigma*sqrt(2e3*e0*E*mass_u*amu)*g*p/babs + q*e0*psi

    return HamiltonianCoordinate(E,mu,Pphi,amu,q)
end

function Base.show(io::IO, c::HamiltonianCoordinate)
    println(io,typeof(c))
    @printf(io," E = %.3f keV\n",c.energy)
    @printf(io," μ₀ = %.3f\n",c.mu)
    @printf(io," Pᵩ = %.3f\n",c.p_phi)
end
