function r_condition(u,t,integ)
    v = integ.f(u,integ.p,t)
    v[1]/sqrt(v[1]^2 +v[3]^2)
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
            M = integ.f.f.M
            gcp = integ.f.f.gcp
            stat.pm = get_pitch(M,gcp,integ.u[4],integ.u[5],integ.u[1],integ.u[3])
        end
    else
        integ.p && terminate!(integ)
    end
end
r_cb = ContinuousCallback(r_condition,r_affect!,abstol=1e-6)

function poloidal_condition(u,t,integ)
    stat = integ.f.f.stat
    i = stat.initial_dir
    (t != 0.0)*(u[i] - stat.ri[i])
end
function poloidal_affect!(integ)
    stat = integ.f.f.stat
    vi = stat.vi
    ri = stat.ri
    vc = integ.f(integ.u,integ.p,integ.t)
    dp = dot(vi,vc)/(norm(vi)*norm(vc))
    vi_rz = SVector(vi[1],vi[3])
    vc_rz = SVector(vc[1],vc[3])
    dprz = dot(vi_rz,vc_rz)/(norm(vi_rz)*norm(vc_rz))
    #println("poloidal callback ",integ.t*1e6," ",integ.u," ",lr," ",lz)
    #println(dp," ",dprz)
    if !stat.poloidal_complete && (stat.nr >= 2) &&
        (dp > 0.99 && dprz > 0.99) &&
        (abs(integ.u[1]-ri[1]) < 0.01)

        stat.poloidal_complete=true
        stat.tau_p = integ.t
        stat.tau_t = 2pi*stat.tau_p/abs(integ.u[2] - integ.sol.u[1][2])
        stat.class = :unknown
        integ.p && terminate!(integ)
    end
end
pol_cb = ContinuousCallback(poloidal_condition,poloidal_affect!,abstol=1e-6)

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

transit_callback = CallbackSet(r_cb, pol_cb, oob_cb)
