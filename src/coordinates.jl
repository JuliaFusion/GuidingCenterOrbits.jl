abstract type AbstractOrbitCoordinate{T} end

immutable EPRCoordinate{T} <: AbstractOrbitCoordinate{T}
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

function EPRCoordinate(M::AxisymmetricEquilibrium, energy, pitch, R ; amu=H2_amu, q=1, dz=0.2)
    zaxis = M.axis[2]
    zmax = zaxis + dz
    zmin = zaxis - dz
    res = optimize(x->M.psi_rz[R,x], zmin, zmax)
    Z = Optim.minimizer(res)
    (Z == zmax || Z == zmin) && error(@sprintf("Unable to find starting Z value with dz = %.2f. Increase dz",dz))
    return EPRCoordinate(energy, pitch, R, Z, amu, q)
end

function Base.show(io::IO, c::EPRCoordinate)
    println(io,typeof(c))
    @printf(io," E = %.3f keV\n",c.energy)
    @printf(io," pitch = %.3f\n",c.pitch)
    @printf(io," Rmax = %.3f m\n",c.r)
end

immutable HamiltonianCoordinate{T} <: AbstractOrbitCoordinate{T}
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

function normalized_hamiltonian(M::AxisymmetricEquilibrium, E, p, r, z; amu=H2_amu, q=1)
    psi = M.psi_rz[r,z]
    babs = M.b[r,z]
    g = M.g[psi]

    mu = abs(M.B0)*(1-p^2)/babs
    Pphi = M.sigma*(-M.sigma*sqrt(2e3*e0*E*mass_u*amu)*g*p/babs + q*e0*psi)/(e0*M.psi[end])

    return E, Pphi, mu
end

function normalized_hamiltonian(M::AxisymmetricEquilibrium, hc::HamiltonianCoordinate)
    E = hc.energy
    mu = (abs(M.B0)/(E*e0*1e3))*hc.mu
    Pphi = (M.sigma/(e0*M.psi[end]))*hc.p_phi

    return E, Pphi, mu
end

function Base.show(io::IO, c::HamiltonianCoordinate)
    println(io,typeof(c))
    @printf(io," E = %.3f keV\n",c.energy)
    @printf(io," μ₀ = %.3e\n",c.mu)
    @printf(io," Pᵩ = %.3e\n",c.p_phi)
end
