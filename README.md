# DiffEqBridge.jl

[![Build Status](https://travis-ci.org/JuliaDiffEq/DiffEqBridge.jl.svg?branch=master)](https://travis-ci.org/JuliaDiffEq/DiffEqBridge.jl)
[![Coverage Status](https://coveralls.io/repos/JuliaDiffEq/DiffEqBridge.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/JuliaDiffEq/DiffEqBridge.jl?branch=master)
[![codecov.io](http://codecov.io/github/JuliaDiffEq/DiffEqBridge.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaDiffEq/DiffEqBridge.jl?branch=master)

This package contains bindings for Bridge.jl to allow it to be used with the
JuliaDiffEq common interface. For more information on using the solvers from this
package, see the [DifferentialEquations.jl documentation](https://juliadiffeq.github.io/DiffEqDocs.jl/latest/).

## Common API Usage

This library adds the common interface to Bridge.jl's solvers. [See the DifferentialEquations.jl documentation for details on the interface](http://docs.juliadiffeq.org/latest/index.html). Following the Black-Scholes example from [the SDE tutorial](http://docs.juliadiffeq.org/latest/tutorials/ode_example.html), we can solve this using `BridgeEuler` via the following:

```julia
α=1
β=1
u0=1/2
f(t,u) = α*u
g(t,u) = β*u
dt = 1//2^(4)
tspan = (0.0,1.0)
prob = SDEProblem(f,g,u0,(0.0,1.0))
sol = solve(prob,BridgeEuler(),dt=dt)
using Plots; plot(sol,vars=(1,2,3))
```

The options available in `solve` are documented [at the common solver options page](http://docs.juliadiffeq.org/latest/basics/common_solver_opts.html). The available methods are documented [at the ODE solvers page](http://docs.juliadiffeq.org/latest/solvers/ode_solve.html#BridgeDiffEq.jl-1)
and [at the SDE solvers page](http://docs.juliadiffeq.org/latest/solvers/sde_solve.html#BridgeDiffEq.jl-1).
