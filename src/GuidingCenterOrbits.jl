__precompile__()

module GuidingCenterOrbits

#import standard libraries
using Printf
using LinearAlgebra

using Equilibrium
using OrdinaryDiffEq
using Optim
using Base.Iterators
using Interpolations
using StaticArrays
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
const B5_amu = 10.81 # amu
const C6_amu = 12.011 # amu

include("particles.jl")
export AbstractParticle, Particle
export Electron, Proton, Deuteron, Triton, Alpha
export GCParticle, EPRParticle
export GCElectron, GCProton, GCDeuteron, GCTriton, GCAlpha

include("coordinates.jl")
export AbstractOrbitCoordinate, EPRCoordinate, HamiltonianCoordinate

include("callbacks.jl")
export standard_callback

include("orbit.jl")
export Orbit, OrbitPath, integrate, get_orbit, gc_velocity

include("fullorbit.jl")
export FullOrbitPath

include("utils.jl")
export get_pitch, hits_wall_path, hits_wall, get_kinetic_energy
export cyclotron_frequency, cyclotron_period, normalize

using ForwardDiff
using ForwardDiff: Dual, partials, value
using DiffResults, DiffEqDiffTools

include("jacobian.jl")
export eprz_to_eprt, eprt_to_eprz, get_jacobian, transform

end
