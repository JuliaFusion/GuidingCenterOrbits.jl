function r_condition(u,t,integ)
    v = integ.f(u,integ.p,integ.t)
    v[1]/sqrt(v[1]^2 +v[2]^2)
end
function r_affect!(integ)
    os = integ.f.f.os
    #println("r callback ",integ.t*1e6," ",integ.u)
    if !os.poloidal_complete
        os.nr += 1
        if (integ.u[1] > os.rm)
            os.rm = integ.u[1]
            os.zm = integ.u[3]
            os.tm = integ.t
        end
    else
        integ.p && terminate!(integ)
    end
end
r_cb = ContinuousCallback(r_condition,r_affect!,abstol=1e-8)

function phi_condition(u,t,integ)
    #dz_cur = get_du(integ)[3]
    #dz_prev = integ.f(integ.uprev,integ.p,integ.tprev)[3]
    #sign(dz_cur) != sign(dz_prev)
    v = integ.f(u,integ.p,integ.t)
    v[2]/norm(v)
end
function phi_affect!(integ)
    os = integ.f.f.os
    #println("phi callback ",integ.t*1e6," ",integ.u)
    if !os.poloidal_complete
        os.nphi += 1
    else
        integ.p && terminate!(integ)
    end
end
#z_cb = DiscreteCallback(z_condition,z_affect!,save_positions=(false,false))
phi_cb = ContinuousCallback(phi_condition,phi_affect!,rootfind=false)

function poloidal_condition(u,t,integ)
    os = integ.f.f.os
    i = os.initial_dir
    u[i] - os.ri[i]
end
function poloidal_affect!(integ)
    os = integ.f.f.os
    vi = os.vi
    vc = integ.f(integ.u,integ.p,integ.t)
    dp = dot(vi,vc)/(norm(vi)*norm(vc))
    vi_rz = SVector(vi[1],vi[3])
    vc_rz = SVector(vc[1],vc[3])
    dprz = dot(vi_rz,vc_rz)/(norm(vi_rz)*norm(vc_rz))
    #println("poloidal callback ",integ.t*1e6," ",integ.u," ",lr," ",lz)
    if !os.poloidal_complete && (os.nr >= 2 && (os.naxis == 0 || os.naxis >= 2) &&
                                 dp > 0.99999 && dprz > 0.99999)
    #if !os.poloidal_complete && ((lr != 0 && lz != 0 && iseven(lr) && iseven(lz)) || lz > 10 || lr > 10)
        os.poloidal_complete=true
        os.tau_p = integ.t
        os.tau_t = 2pi*os.tau_p/abs(integ.u[2] - integ.sol.u[1][2])
        M = integ.f.f.M
        oc = integ.f.f.oc
        os.pm = get_pitch(M,oc,os.rm,os.zm)
        if os.nphi > 0
            if os.naxis == 0 || os.naxis == 4
                os.class = :trapped
            elseif os.naxis == 2
                os.class = :potato
            else
                os.class = :unknown
            end
        else
            if os.naxis == 2
                if os.pm > 0
                    os.class = :co_passing
                else
                    os.class = :ctr_passing
                end
            elseif os.naxis == 0
                os.class = :stagnation
            else
                os.class = :unknown
            end
        end
        integ.p && terminate!(integ)
    end
end
pol_cb = ContinuousCallback(poloidal_condition,poloidal_affect!,abstol=1e-8)

function maxis_condition(u,t,integ)
    raxis = integ.f.f.M.axis[1]
#    sign(u[1] - raxis) != sign(integ.uprev[1] - raxis)
    u[1] - raxis
end
function maxis_affect!(integ)
    os = integ.f.f.os
    #println("maxis callback ",integ.t*1e6," ",integ.u)
    if !os.poloidal_complete
        os.naxis += 1
    else
        integ.p && terminate!(integ)
    end
end
#maxis_cb = DiscreteCallback(maxis_condition,maxis_affect!,save_positions=(false,false))
maxis_cb = ContinuousCallback(maxis_condition,maxis_affect!,rootfind=true,abstol=1e-8)

function out_of_bounds_condition(u,t,integ)
    M = integ.f.f.M
    !((M.r[1] < u[1] < M.r[end]) && (M.z[1] < u[3] < M.z[end]))
end
function out_of_bounds_affect!(integ)
    integ.f.f.os.hits_boundary=true
    integ.f.f.os.class = :lost
    terminate!(integ)
end
oob_cb = DiscreteCallback(out_of_bounds_condition, out_of_bounds_affect!,save_positions=(false,false))

standard_callback = CallbackSet(r_cb, phi_cb, maxis_cb, pol_cb, oob_cb)
