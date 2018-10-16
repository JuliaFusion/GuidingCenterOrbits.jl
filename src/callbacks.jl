function r_condition(u,t,integ)
    integ.f(u,integ.p,t)[1]
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
r_cb = ContinuousCallback(r_condition,r_affect!,abstol=1e-10,save_positions=(false,false))

function z_condition(u,t,integ)
    #dz_cur = get_du(integ)[3]
    #dz_prev = integ.f(integ.uprev,integ.p,integ.tprev)[3]
    #sign(dz_cur) != sign(dz_prev)
    integ.f(u,integ.p,integ.t)[3]
end
function z_affect!(integ)
    os = integ.f.f.os
    #println("z callback ",integ.t*1e6," ",integ.u)
    if !os.poloidal_complete
        os.nz += 1
    else
        integ.p && terminate!(integ)
    end
end
#z_cb = DiscreteCallback(z_condition,z_affect!,save_positions=(false,false))
z_cb = ContinuousCallback(z_condition,z_affect!,abstol=1e-10,save_positions=(false,false))

function poloidal_condition(u,t,integ)
    #i = integ.f.f.os.initial_dir
    u[3] - integ.sol.u[1][3]
end
function poloidal_affect!(integ)
    os = integ.f.f.os
    lr = os.nr
    lz = os.nz
    #println("poloidal callback ",integ.t*1e6," ",integ.u," ",lr," ",lz)
    if !os.poloidal_complete && ((lr != 0 && lz != 0 && iseven(lr) && iseven(lz)) || lz > 10 || lr > 10)
        os.poloidal_complete=true
        os.tau_p = integ.t
        os.tau_t = 2pi*os.tau_p/abs(integ.u[2] - integ.sol.u[1][2])
        M = integ.f.f.M
        oc = integ.f.f.oc
        os.pm = get_pitch(M,oc,os.rm,os.zm)
        if os.pm < 0.0
            os.class = :ctr_passing
        else
            if os.naxis == 4 || os.naxis == 0
                if os.nr >= 4
                    os.class = :trapped
                else
                    os.class = :stagnation
                end
            else
                if os.nr >= 4
                    os.class = :potato
                else
                    os.class = :co_passing
                end
            end
        end
        integ.p && terminate!(integ)
    end
end
pol_cb = ContinuousCallback(poloidal_condition,poloidal_affect!,abstol=1e-12,save_positions=(false,false))

function maxis_condition(u,t,integ)
    raxis = integ.f.f.M.axis[1]
    sign(u[1] - raxis) != sign(integ.uprev[1] - raxis)
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
maxis_cb = DiscreteCallback(maxis_condition,maxis_affect!,save_positions=(false,false))
#maxis_cb = ContinuousCallback(maxis_condition,maxis_affect!,save_positions=(false,false),rootfind=true)

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

standard_callback = CallbackSet(r_cb, z_cb, maxis_cb, pol_cb, oob_cb)
