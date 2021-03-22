"""
    freeze_condition(u,t,integ)

Take the current step u, the current time step t, the previous step integ.uprev
and the previous time step integ.tprev, and compute the spatial and temporal difference.
If the spatial difference is less than a thousandth of 1Å, or if the temporal difference
is less than a thousandth of a femtosecond, return true. Otherwise return false.

Inputs:
    - u = the next step suggested by the integrator. u is an array containing the elements
    specified in GuidingCenterOrbits.jl/orbit.jl/integrate(). The elements of u fully characterize
    the guiding-centre particle. The elements of u are [r, phi, z, p_para, mu] where r,phi,z are the
    cylindrical coordinates, p_para is the momentum parallel to the magnetic field and mu is the
    magnetic moment.
    - t = the current time step of the integrator
    - integ = the integrator. Read more at the 'Integrator interface' page of the DifferentialEquations.jl package website.
"""
function freeze_condition(u,t,integ)
    Δu = abs.(u .- integ.uprev) # The difference between the current and previous time-step (R,phi,Z,p_para,mu)
    Δr = abs(Δu[1])
    Δphi = abs(Δu[2])
    Δz = abs(Δu[3])
    Δspace = sqrt(Δr^2 + Δphi^2 + Δz^2)
    #println(Δspace)
    Δtime = abs(t - integ.tprev)
    #println(Δtime)
    (Δspace < 1.0e-13) ||  (Δtime < 1.0e-18) # If the integration has stagnated in space or time...
end
function freeze_affect!(integ)
    #println("I did it!")
    terminate!(integ) # If the freeze_condition returs true, terminate the integration!
end
# Read more about Callbacks at the DifferentialEquations.jl 'Event handling and Callback Functions' page
brr_cb = DiscreteCallback(freeze_condition,freeze_affect!,save_positions=(false,false)) # save_positions = (false,false). If freeze_condition, do not save position before or after suggested time step.

"""
    phi_condition(u,t,integ)

Check whether the guiding-centre particle has gone around the tokamak toroidally more than
maxphi in phi angle, without completing even one poloidal orbit. If so, return true.

Inputs:
    Same as for freeze_condition.
"""
function phi_condition(maxphi, u, t, integ)
    stat = integ.f.f.stat # The GCStatus mutable struct defined in orbit.jl
    !stat.poloidal_complete && (abs(u[2]-stat.ri[2]) > maxphi) # If it has gone around the tokamak more than (maxphi/(2*pi)) times toroidally, but not once poloidally...
end
function phi_affect!(integ)
    integ.p && terminate!(integ) # If only one poloidal transit desired, terminate the integration
end
phi_callback(maxphi) = DiscreteCallback((u,t,integ)->phi_condition(maxphi,u,t,integ),phi_affect!,save_positions=(false,false))
# Check callback syntax at https://diffeq.sciml.ai/stable/features/callback_functions/ for further info

"""
    r_condition(u,t,integ)

Return the ratio of the r component of the velocity, to the total speed of the particle.

Inputs:
    Same as freeze_condition.
"""
function r_condition(u,t,integ)
    v = integ.f(u,integ.p,t)
    v[1]/sqrt(v[1]^2 +v[3]^2)
end
function r_affect!(integ)
    stat = integ.f.f.stat # Get the stat struct. Please see GuidingCenterOrbits.jl/orbit.jl/GCStatus{} for fields.
    if !stat.poloidal_complete && stat.nr < 20 # prevent infinite loop
        stat.nr += 1
        if (integ.u[1] > stat.rm) # If the guiding-centre particles has gone outside its previous rm coordinate
            stat.rm = integ.u[1] # Update the rm coordinate with the current r coordinate
            stat.zm = integ.u[3] # Update the zm coordinate with the current z coordinate
            stat.tm = integ.t # Update the poloidal transit time with the current time
            M = integ.f.f.M # Get the axisymmetric equilibrium
            gcp = integ.f.f.gcp # Get the guiding-centre particle itself
            stat.pm = get_pitch(M,gcp,integ.u[4],integ.u[5],integ.u[1],integ.u[3]) # Update the pm coordinate with the current pitch
        end
    else
        integ.p && terminate!(integ)
    end
end
r_cb = ContinuousCallback(r_condition,r_affect!,abstol=1e-6)

"""
    poloidal_condition(u,t,integ)

If it's not the first time step of the integration, return the difference between
the current speed and the initial speed, both evaluated in the initial velocity
direction of the guiding-centre particle.

Inputs:
    Same as freeze_condition.
"""
function poloidal_condition(u,t,integ)
    stat = integ.f.f.stat
    i = stat.initial_dir
    (t != 0.0)*(u[i] - stat.ri[i])
end
function poloidal_affect!(integ)
    stat = integ.f.f.stat # Get the stat struct. Please see GuidingCenterOrbits.jl/orbit.jl/GCStatus{} for fields.
    vi = stat.vi # Get the initial velocity vector
    ri = stat.ri # Get the initial guiding-centre elements vector
    vc = integ.f(integ.u,integ.p,integ.t) # Get the current velocity vector
    dp = dot(vi,vc)/(norm(vi)*norm(vc)) # Compute the normalized dot product between the current velocity vector and the initial velocity vector
    vi_rz = SVector(vi[1],vi[3]) # Create a 2D vector which is the initial velocity vector projected onto the cross-section of the tokamak
    vc_rz = SVector(vc[1],vc[3]) # Create a 2D vector which is the current velocity vector projected onto the cross-section of the tokamak
    dprz = dot(vi_rz,vc_rz)/(norm(vi_rz)*norm(vc_rz)) # Compute the normalized dot product between the 2D projected initial and current velocities
    if !stat.poloidal_complete && (stat.nr >= 2) &&
        (dp > 0.99 && dprz > 0.99) && # If the dot products are almost 1...
        (abs(integ.u[1]-ri[1]) < 0.01) # ... and the particle is within 1 cm of the initial r coordinate

        stat.poloidal_complete=true # Deem the particle as poloidally complete
        stat.tau_p = integ.t # Set the current time to be the poloidal transit time
        stat.tau_t = 2pi*stat.tau_p/abs(integ.u[2] - integ.sol.u[1][2]) # Compute the toroidal transit time
        stat.class = :unknown # We don't yet know what class the orbit is (potato, banana... etc)
        integ.p && terminate!(integ) # Terminate the integration, if only one poloidal transit is desired
    end
end
pol_cb = ContinuousCallback(poloidal_condition,poloidal_affect!,abstol=1e-6)

"""
    out_of_bounds_condition(u,t,integ)

Check whether the guiding-centre particle is still inside the boundaries provided
by the axisymmetric equilibrium M. If that is not the case, return true.

Inputs:
    Same as freeze_condition.
"""
function out_of_bounds_condition(u,t,integ)
    M = integ.f.f.M
    !((M.r[1] < u[1] < M.r[end]) && (M.z[1] < u[3] < M.z[end]))
end
function out_of_bounds_affect!(integ)
    integ.f.f.stat.hits_boundary=true # The particle hit the boundary
    integ.f.f.stat.class = :lost # Classify it as lost
    terminate!(integ) # Terminate the integration
end
oob_cb = DiscreteCallback(out_of_bounds_condition, out_of_bounds_affect!,save_positions=(false,false))

"""
    wall_condition(u,t,integ)

Check whether the guiding-centre particle has hit the tokamak wall. If that is the case,
return true.

Inputs:
    Same as freeze_condition.
"""
function wall_condition(wall, u, t, integ)
    !in_vessel(wall,(u[1],u[3]))
end
function wall_affect!(integ)
    integ.f.f.stat.hits_boundary=true
    integ.f.f.stat.class = :lost
    terminate!(integ)
end
wall_callback(wall) = DiscreteCallback((u,t,integ)->wall_condition(wall,u,t,integ),wall_affect!,save_positions=(false,false))

transit_callback = CallbackSet(r_cb, pol_cb, oob_cb, brr_cb) # A collection of callbacks
