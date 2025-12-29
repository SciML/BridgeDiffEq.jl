using PrecompileTools

@setup_workload begin
    f(u, p, t) = u
    g(u, p, t) = u
    u0 = 0.5
    dt = 1 // 16
    tspan = (0.0, 1.0)

    @compile_workload begin
        prob_ode = ODEProblem(f, u0, tspan)
        solve(prob_ode, BridgeR3(), dt = dt)
        solve(prob_ode, BridgeBS3(), dt = dt)

        prob_sde = SDEProblem(f, g, u0, tspan)
        solve(prob_sde, BridgeEuler(), dt = dt)
        solve(prob_sde, BridgeHeun(), dt = dt)
        solve(prob_sde, BridgeSRK(), dt = dt)
    end
end
