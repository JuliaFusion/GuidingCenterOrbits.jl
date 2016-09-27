function B_function(f_psi,f_g)
    res = ForwardDiff.GradientResult(rand(2))
    function B(ri)
        ForwardDiff.gradient!(res,f_psi, ri)
        psi = ForwardDiff.value(res)
        grad_psi = -ForwardDiff.gradient(res) #negative because "ribbon" polodial flux definintion

        br = -grad_psi[2]/ri[1]
        bz = grad_psi[1]/ri[1]
        bt = f_g([psi])/ri[1]

        return [br,bt,bz]
    end
end

function J_function(f_psi,f_g, f_p)
    res = ForwardDiff.GradientResult(rand(2))
    resd = ForwardDiff.DerivativeResult(1.0)
    function J(ri)
        ForwardDiff.gradient!(res,f_psi, ri)
        psi = ForwardDiff.value(res)
        grad_psi = -ForwardDiff.gradient(res) #negative because "ribbon" polodial flux definintion

        g = f_g([psi])
        resd = ForwardDiff.derivative!(resd, f_g, psi)
        gp = -ForwardDiff.derivative(resd)
        resd = ForwardDiff.derivative!(resd, f_p, psi)
        pp = -ForwardDiff.derivative(resd)

        jr = -gp*grad_psi[2]/(ri[1]*mu0)
        jz = gp*grad_psi[1]/(ri[1]*mu0)
        jt = ri[1]*pp + g*gp/(ri[1]*mu0)

        return [jr,jt,jz]
    end
end

type AxisymmetricEquilibrium{T}
    r_domain::NTuple{2,T}   #(rmin, rmax)
    z_domain::NTuple{2,T}   #(zmin, zmax)
    psi_domain::NTuple{2,T} #(psimin, psimax)

    psi::Function  # "Ribbon" Polodial Flux (polodial flux from magnetic axis)
    g::Function    # Polodial Current
    b::Function    # Magnetic field magnitude
    p::Function    # Plasma pressure
    B::Function    # Magnetic Field Vector
    J::Function    # Current Density Vector
    axis::NTuple{2, T} #Magnetic Axis (raxis,zaxis)
    sigma::Int     #sign(dot(J,B))

end

function AxisymmetricEquilibrium{S}(r_domain::NTuple{2,S},
                                    z_domain::NTuple{2,S},
                                    psi_domain::NTuple{2,S},
                                    psi, g, p, axis::NTuple{2,S})
    B = B_function(psi, g)
    J = J_function(psi, g, p)

    r = linspace(r_domain...,100)
    z = linspace(z_domain...,100)
    babs = [norm(B([rr,zz])) for rr in r, zz in z]
    b = interpolate(babs, (r,z), BicubicSpline(), Irregular{2,Val{size(babs)}}, Value(0.0))

    ri = collect(axis)  + [r_domain[2] - axis[1], z_domain[2] - axis[2]]/10

    sigma = sign(dot(J(ri), B(ri)))

    AxisymmetricEquilibrium{S}(r_domain,z_domain,psi_domain, psi, g, b, p, B, J, axis, sigma)
end
