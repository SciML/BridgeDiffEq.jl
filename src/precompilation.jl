using PrecompileTools: PrecompileTools, @compile_workload, @setup_workload
using DiffEqBase: ODEProblem, SDEProblem

@setup_workload begin
    f(u, p, t) = u
    g(u, p, t) = u
    u0 = 0.5
    dt = 1 // 16
    tspan = (0.0, 1.0)

    @compile_workload begin
        # Scalar problems
        prob_ode = ODEProblem(f, u0, tspan)
        solve(prob_ode, BridgeR3(), dt = dt)
        solve(prob_ode, BridgeBS3(), dt = dt)

        prob_sde = SDEProblem(f, g, u0, tspan)
        solve(prob_sde, BridgeEuler(), dt = dt)
        solve(prob_sde, BridgeHeun(), dt = dt)
        solve(prob_sde, BridgeSRK(), dt = dt)

        # SVector problems (common vector use case)
        u0_vec = StaticArrays.@SVector [1.0, 1.0]
        prob_ode_vec = ODEProblem(f, u0_vec, tspan)
        solve(prob_ode_vec, BridgeR3(), dt = dt)
        solve(prob_ode_vec, BridgeBS3(), dt = dt)

        prob_sde_vec = SDEProblem(f, g, u0_vec, tspan)
        solve(prob_sde_vec, BridgeEuler(), dt = dt)
        solve(prob_sde_vec, BridgeHeun(), dt = dt)
    end
end
