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
include("contour.jl")
include("equilibrium.jl")
include("orbit.jl")

export Polygon, in_polygon, orientation, area
export AxisymmetricEquilibrium
export follow_contour
export calc_orbit, plot_orbit, get_pitch, Orbit
export CQL3DCoordinate, HamiltonianCoordinate, EPRZCoordinate

end
