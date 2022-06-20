__precompile__()

module GuidingCenterOrbits

#import standard libraries
using Printf
using LinearAlgebra
using Statistics

using Equilibrium
using OrdinaryDiffEq
using ForwardDiff
using Optim
using Base.Iterators
using Interpolations
using StaticArrays
using Contour
using HDF5
using ProgressMeter
#using Plots

const e0 = 1.60217733e-19 # Coulombs / Joules
const mu0 = 4*pi*1e-7 # N/A^2
const c0 = 2.99792458e8 # m/s

const mass_u = 1.6605402e-27 # kg
const e_amu = 5.48579909070e-4 # amu
const H1_amu = 1.007276466879 # amu
const H2_amu = 2.0141017778 # amu
const H3_amu = 3.01550071632 # amu
const He3_amu = 3.01602931914 # amu
const He4_amu = 4.00150617913 # amu
const B5_amu = 10.81 # amu
const C6_amu = 12.011 # amu

include("particles.jl")
export AbstractParticle, Particle
export Electron, Proton, Deuteron, Triton, Alpha
export GCParticle, EPRParticle
export GCElectron, GCProton, GCDeuteron, GCTriton, GCAlpha

include("coordinates.jl")
export AbstractOrbitCoordinate, EPRCoordinate, HamiltonianCoordinate, GCEPRCoordinate

include("callbacks.jl")
export standard_callback, wall_callback, phi_callback 

include("orbit.jl")
export Orbit, OrbitPath, integrate, get_orbit, gc_velocity, GCStatus, write_Orbs, read_Orbs

include("fullorbit.jl")
export FullOrbitPath, get_full_orbit

include("projection.jl")
export OrbitProjection, orbit_projection, in_orbit

include("utils.jl")
export get_pitch, hits_wall_path, hits_wall, get_kinetic_energy
export lorentz_factor, cyclotron_frequency, cyclotron_period, normalize
export perpendicular_vectors, gyro_step, larmor_radius, gcde_check

using ForwardDiff
using ForwardDiff: Dual, partials, value

include("jacobian.jl")
export get_jacobian, transform

end
