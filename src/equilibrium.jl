struct AxisymmetricEquilibrium{T<:Real, S<:Range{T},
                               R<:AbstractArray{T,2},
                               Q<:AbstractArray{T,1}}
    r::S               # Radius/R range
    z::S               # Elevation/Z range
    psi::S             # "Ribbon" Polodial Flux range (polodial flux from magnetic axis)
    psi_rz::R          # "Ribbon" Polodial Flux on RZ grid (polodial flux from magnetic axis)
    g::Q               # Polodial Current
    p::Q               # Plasma pressure
    b::R               # Magnetic field magnitude
    j::R               # Plasma Current magnitude
    axis::NTuple{2, T} # Magnetic Axis (raxis,zaxis)
    sigma::Int         # sign(dot(J,B))
end

function AxisymmetricEquilibrium(r::Range{T}, z::Range{T}, psi::Range{T}, psi_rz, g, p, axis::NTuple{2,T}) where {T <: Real}

    psi_max = maximum(psi_rz)
    dpsi = step(psi)
    psi_ext = psi[1]:dpsi:psi_max
    g_ext = [i <= length(g) ? g[i] : g[end] for i=1:length(psi_ext)]
    p_ext = [i <= length(p) ? p[i] : 0.0 for i=1:length(psi_ext)]

    psi_rz_itp = scale(interpolate(psi_rz, BSpline(Cubic(Line())), OnGrid()), r, z)
    g_itp = scale(interpolate(g_ext, BSpline(Cubic(Line())), OnGrid()), psi_ext)
    p_itp = scale(interpolate(p_ext, BSpline(Cubic(Line())), OnGrid()), psi_ext)

    b = [norm(Bfield(psi_rz_itp,g_itp,rr,zz)) for rr in r, zz in z]
    b_itp = scale(interpolate(b, BSpline(Cubic(Line())), OnGrid()), r, z)

    j = [norm(Jfield(psi_rz_itp,g_itp,p_itp,rr,zz)) for rr in r, zz in z]
    j_itp = scale(interpolate(j, BSpline(Cubic(Line())), OnGrid()), r, z)

    rr = axis[1] + (r[end] - axis[1])/10
    zz = axis[2] + (z[end] - axis[2])/10

    sigma = Int(sign(dot(Jfield(psi_rz_itp,g_itp,p_itp,rr,zz), Bfield(psi_rz_itp,g_itp,rr,zz))))

    AxisymmetricEquilibrium(r, z, psi_ext, psi_rz_itp, g_itp, p_itp, b_itp, j_itp, axis, sigma)
end

function Bfield(psi_rz, g, r, z)
    psi = psi_rz[r,z]
    gval = g[psi]
    grad_psi = -gradient(psi_rz, r, z) #negative because "ribbon" polodial flux definintion

    br = -grad_psi[2]/r
    bz = grad_psi[1]/r
    bt = gval/r

    return [br,bt,bz]
end

function Bfield(M::AxisymmetricEquilibrium, r, z)
    Bfield(M.psi_rz, M.g, r, z)
end

function Jfield(psi_rz, g, p, r, z)
    psi = psi_rz[r,z]
    gval = g[psi]
    grad_psi = -gradient(psi_rz, r, z) #negative because "ribbon" polodial flux definintion

    gp = -gradient(g, psi)[1]
    pp = -gradient(p, psi)[1]

    jr = -gp*grad_psi[2]/(r*mu0)
    jz = gp*grad_psi[1]/(r*mu0)
    jt = r*pp + gval*gp/(r*mu0)

    return [jr,jt,jz]
end

function Jfield(M::AxisymmetricEquilibrium, r, z)
    Jfield(M.psi_rz, M.g, M.p, r, z)
end
