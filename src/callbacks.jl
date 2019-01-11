function angle_condition(u,t,integ)
    ri = integ.sol.u[1]
    p0 = ri[1] + ri[3]*im
    p0_s = ri[1] + (ri[3] - 1.0)*im
    p1 = u[1] + u[3]*im
    p1_s = u[1] + (u[3] - 1.0)*im
    d = p1/p0
    theta = angle(d)
    d_s = p1_s/p0_s
    theta_s = angle(d_s)
    (t != 0.0)*(theta*theta_s + theta + theta_s)
end
function angle_affect!(integ)
    stat = integ.f.f.stat
    vi = stat.vi
    ri = stat.ri
    vc = integ.f(integ.u,integ.p,integ.t)
    vi_rz = SVector(vi[1],vi[3])
    vc_rz = SVector(vc[1],vc[3])
    dprz = dot(vi_rz,vc_rz)/(norm(vi_rz)*norm(vc_rz))
    if !stat.poloidal_complete && dprz > 0.95 && abs(integ.u[1]-ri[1]) < 0.01
        stat.poloidal_complete=true
        stat.tau_p = integ.t
        stat.tau_t = 2pi*stat.tau_p/abs(integ.u[2] - integ.sol.u[1][2])
        M = integ.f.f.M
        oc = integ.f.f.oc
        stat.pm = get_pitch(M,oc,stat.rm,stat.zm)
        stat.class = :unknown
        integ.p && terminate!(integ)
    end
end
angle_cb = ContinuousCallback(angle_condition,angle_affect!,abstol=1e-6)

function r_condition(u,t,integ)
    v = integ.f(u,integ.p,integ.t)
    v[1]/sqrt(v[1]^2 + v[3]^2)
end
function r_affect!(integ)
    stat = integ.f.f.stat
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
r_cb = ContinuousCallback(r_condition,r_affect!,abstol=1e-6)

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

transit_callback = CallbackSet(angle_cb, r_cb, oob_cb)
