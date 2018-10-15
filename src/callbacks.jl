function r_condition(u,t,integ)
    integ.f(u,integ.p,t)[1]
end
function r_affect!(integ)
    os = integ.f.f.os
    if !os.poloidal_complete
        os.nr += 1
        if (integ.u[1] > os.rm)
            os.rm = integ.u[1]
            os.zm = integ.u[3]
            os.tm = integ.t
        end
    end
end
r_cb = ContinuousCallback(r_condition,r_affect!,abstol=1e-8,save_positions=(false,false))

function z_condition(u,t,integ)
    integ.f(u,integ.p,t)[3]
end
function z_affect!(integ)
    os = integ.f.f.os
    if !os.poloidal_complete
        os.nz += 1
    end
end
z_cb = ContinuousCallback(z_condition,z_affect!,abstol=1e-8.save_positions=(false,false),rootfind=false)

function poloidal_condition(u,t,integ)
    u[3] - integ.sol.u[1][3]
end
function poloidal_affect!(integ)
    os = integ.f.f.os
    lr = os.nr
    lz = os.nz
    if !os.poloidal_complete && lr != 0 && lz != 0 && iseven(lr) && iseven(lz) || lz > 10
        os.poloidal_complete=true
        os.tau_p = integ.t
        os.tau_t = 2pi*os.tau_p/abs(integ.u[2] - integ.sol.u[1][2])
        M = integ.f.f.M
        oc = integ.f.f.oc
        os.pm = get_pitch(M,oc,os.rm,os.zm)
        if os.naxis == 0
            os.class = :stagnation
            if os.nr >= 4
                os.class = :trapped
            end
        elseif os.naxis == 2
            if os.nr >= 4
                os.class = :potato
            else
                if os.pm > 0
                    os.class = :co_passing
                else
                    os.class = :ctr_passing
                end
            end
        elseif os.naxis == 4
            os.class = :trapped
        end
        integ.p && terminate!(integ)
    end
end
pol_cb = ContinuousCallback(poloidal_condition,poloidal_affect!,save_positions=(false,false))

function maxis_condition(u,t,integ)
    M = integ.f.f.M
    u[1] - M.axis[1]
end
function maxis_affect!(integ)
    os = integ.f.f.os
    if !os.poloidal_complete
        os.naxis += 1
    end
end
maxis_cb = ContinuousCallback(maxis_condition,maxis_affect!,save_positions=(false,false),rootfind=false)

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

standard_callback = CallbackSet(r_cb, z_cb, pol_cb, maxis_cb, oob_cb)
