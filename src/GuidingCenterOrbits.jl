__precompile__()

module GuidingCenterOrbits

using Equilibrium
using Interpolations
using Sundials
using Optim
using Base.Iterators
#using Plots

const mass_u = 1.6605402e-27
const e0 = 1.60217733e-19
const mu0 = 4*pi*1e-7
const H2_amu = 2.0141017778

include("coordinates.jl")
include("orbit.jl")
include("utils.jl")

export EPRCoordinate, HamiltonianCoordinate
export normalized_hamiltonian
export get_pitch, hits_wall, get_kinetic_energy
export Orbit, OrbitPath, get_orbit, down_sample, plot_orbit

end
