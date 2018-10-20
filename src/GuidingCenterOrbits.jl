__precompile__()

module GuidingCenterOrbits

#import standard libraries
using Printf
using LinearAlgebra

using Equilibrium
using DifferentialEquations
using Optim
using Base.Iterators
using Interpolations
using StaticArrays
#using Plots

const mass_u = 1.6605402e-27
const e0 = 1.60217733e-19
const mu0 = 4*pi*1e-7
const H2_amu = 2.0141017778

include("coordinates.jl")
export EPRCoordinate, HamiltonianCoordinate, normalized_hamiltonian
export Particle, GCParticle, EPRParticle

include("callbacks.jl")
export standard_callback

include("orbit.jl")
export Orbit, OrbitPath, integrate, get_orbit

include("fullorbit.jl")
export FullOrbitPath

include("utils.jl")
export get_pitch, hits_wall_path, hits_wall, get_kinetic_energy

end
