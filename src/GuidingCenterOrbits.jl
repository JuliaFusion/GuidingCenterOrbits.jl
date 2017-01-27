module GuidingCenterOrbits

using DiffBase
using ForwardDiff
using InterpolatingFunctions
using Contour
using Sundials
using Plots

const mass_u = 1.6605402e-27
const e0 = 1.60217733e-19
const mu0 = 4*pi*1e-7
const H2_amu = 2.0141017778

include("polygons.jl")
include("equilibrium.jl")
include("coordinates.jl")
include("orbit.jl")
include("utils.jl")

export Polygon, in_polygon
export AxisymmetricEquilibrium
export EPRCoordinate, HamiltonianCoordinate
export get_pitch, hits_wall
export Orbit, OrbitPath, get_orbit, plot_orbit

end
