using BridgeDiffEq, StaticArrays
using Base.Test

α=1
β=1
u0=1/2
f(t,u) = α*u
g(t,u) = β*u
dt = 1//2^(4)
tspan = (0.0,1.0)

prob = ODEProblem(f,u0,(0.0,1.0))
sol = solve(prob,BridgeR3(),dt=dt)
sol = solve(prob,BridgeBS3(),dt=dt)

prob = SDEProblem(f,g,u0,(0.0,1.0))
@time sol = solve(prob,BridgeEuler(),dt=dt)
sol = solve(prob,BridgeHeun(),dt=dt)
sol = solve(prob,BridgeSRK(),dt=dt)



u0 = @SVector [2.0,3.0]
prob = ODEProblem(f,u0,(0.0,1.0))
sol = solve(prob,BridgeR3(),dt=dt)
sol = solve(prob,BridgeBS3(),dt=dt)

prob = SDEProblem(f,g,u0,(0.0,1.0))
sol = solve(prob,BridgeEuler(),dt=dt)
sol = solve(prob,BridgeHeun(),dt=dt)
