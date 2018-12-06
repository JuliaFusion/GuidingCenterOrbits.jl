# GuidingCenterOrbits

[![Build Status](https://travis-ci.org/lstagner/GuidingCenterOrbits.jl.svg?branch=master)](https://travis-ci.org/lstagner/GuidingCenterOrbits.jl)

[![Coverage Status](https://coveralls.io/repos/lstagner/GuidingCenterOrbits.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/lstagner/GuidingCenterOrbits.jl?branch=master)

[![codecov.io](http://codecov.io/github/lstagner/GuidingCenterOrbits.jl/coverage.svg?branch=master)](http://codecov.io/github/lstagner/GuidingCenterOrbits.jl?branch=master)

#Basic Usage
```
using EFIT, Equilibrium, GuidingCenterOrbits

#Read in equilibrium
M, wall = read_geqdsk(EFIT.test_gfile)

#Define initial conditions
gcp = GCParticle(80.0,0.2,1.9,0.0)

#Calculate trajectory
path, stat = integrate(M, gcp, tmax=1000.0)

#Calculate orbit trajectory for single poloidal transit
o = get_orbit(M,gcp)
```
