function r_condition(u,t,integ)
    v = integ.f(u,integ.p,integ.t)
    v[1]/sqrt(v[1]^2 +v[2]^2)
end
function r_affect!(integ)
    stat = integ.f.f.stat
    #println("r callback ",integ.t*1e6," ",integ.u)
    if !stat.poloidal_complete
        stat.nr += 1
        if (integ.u[1] > stat.rm)
            stat.rm = integ.u[1]
            stat.zm = integ.u[3]
            stat.tm = integ.t
        end
    else
        integ.p && terminate!(integ)
    end
end
r_cb = ContinuousCallback(r_condition,r_affect!)#,abstol=1e-8)

function phi_condition(u,t,integ)
    #dz_cur = get_du(integ)[3]
    #dz_prev = integ.f(integ.uprev,integ.p,integ.tprev)[3]
    #sign(dz_cur) != sign(dz_prev)
    v = integ.f(u,integ.p,integ.t)
    v[2]/norm(v)
end
function phi_affect!(integ)
    stat = integ.f.f.stat
    #println("phi callback ",integ.t*1e6," ",integ.u)
    if !stat.poloidal_complete
        stat.nphi += 1
    else
        integ.p && terminate!(integ)
    end
end
#z_cb = DiscreteCallback(z_condition,z_affect!,save_positions=(false,false))
phi_cb = ContinuousCallback(phi_condition,phi_affect!,rootfind=false)

function poloidal_condition(u,t,integ)
    stat = integ.f.f.stat
    i = stat.initial_dir
    u[i] - stat.ri[i]
end
function poloidal_affect!(integ)
    stat = integ.f.f.stat
    vi = stat.vi
    vc = integ.f(integ.u,integ.p,integ.t)
    dp = dot(vi,vc)/(norm(vi)*norm(vc))
    vi_rz = SVector(vi[1],vi[3])
    vc_rz = SVector(vc[1],vc[3])
    dprz = dot(vi_rz,vc_rz)/(norm(vi_rz)*norm(vc_rz))
    #println("poloidal callback ",integ.t*1e6," ",integ.u," ",lr," ",lz)
    #println(dp," ",dprz)
    if !stat.poloidal_complete && (stat.nr >= 2 && (stat.naxis == 0 || stat.naxis >= 2) &&
                                   (dp > 0.9999 && dprz > 0.999))
    #if !stat.poloidal_complete && ((lr != 0 && lz != 0 && iseven(lr) && iseven(lz)) || lz > 10 || lr > 10)
        stat.poloidal_complete=true
        stat.tau_p = integ.t
        stat.tau_t = 2pi*stat.tau_p/abs(integ.u[2] - integ.sol.u[1][2])
        M = integ.f.f.M
        oc = integ.f.f.oc
        stat.pm = get_pitch(M,oc,stat.rm,stat.zm)
        if stat.nphi > 0
            if stat.naxis == 0 || stat.naxis == 4
                stat.class = :trapped
            elseif stat.naxis == 2
                stat.class = :potato
            else
                stat.class = :unknown
            end
        else
            if stat.naxis == 2
                if stat.pm > 0
                    stat.class = :co_passing
                else
                    stat.class = :ctr_passing
                end
            elseif stat.naxis == 0
                stat.class = :stagnation
            else
                stat.class = :unknown
            end
        end
        integ.p && terminate!(integ)
    end
end
pol_cb = ContinuousCallback(poloidal_condition,poloidal_affect!)#,abstol=1e-8)

function maxis_condition(u,t,integ)
    raxis = integ.f.f.M.axis[1]
#    sign(u[1] - raxis) != sign(integ.uprev[1] - raxis)
    u[1] - raxis
end
function maxis_affect!(integ)
    stat = integ.f.f.stat
    #println("maxis callback ",integ.t*1e6," ",integ.u)
    if !stat.poloidal_complete
        stat.naxis += 1
    else
        integ.p && terminate!(integ)
    end
end
#maxis_cb = DiscreteCallback(maxis_condition,maxis_affect!,save_positions=(false,false))
maxis_cb = ContinuousCallback(maxis_condition,maxis_affect!,rootfind=true)#,abstol=1e-8)

function out_of_bounds_condition(u,t,integ)
    M = integ.f.f.M
    !((M.r[1] < u[1] < M.r[end]) && (M.z[1] < u[3] < M.z[end]))
end
function out_of_bounds_affect!(integ)
    integ.f.f.stat.hits_boundary=true
    integ.f.f.stat.class = :lost
    terminate!(integ)
end
oob_cb = DiscreteCallback(out_of_bounds_condition, out_of_bounds_affect!,save_positions=(false,false))

function wall_condition(wall, u, t, integ)
    !in_vessel(wall,(u[1],u[3]))
end
function wall_affect!(integ)
    integ.f.f.stat.hits_boundary=true
    integ.f.f.stat.class = :lost
    terminate!(integ)
end
wall_callback(wall) = DiscreteCallback((u,t,integ)->wall_condition(wall,u,t,integ),wall_affect!,save_positions=(false,false))

transit_callback = CallbackSet(r_cb, phi_cb, maxis_cb, pol_cb, oob_cb)
