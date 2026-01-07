using BridgeDiffEq
using JET
using StaticArrays
using Test

@testset "JET static analysis" begin
    # Setup test problems
    f(u, p, t) = u
    g(u, p, t) = u
    u0 = 0.5
    dt = 1 // 16
    tspan = (0.0, 1.0)

    prob_ode = ODEProblem(f, u0, tspan)
    prob_sde = SDEProblem(f, g, u0, tspan)

    @testset "Scalar ODE type stability" begin
        rep = @report_opt target_modules = (BridgeDiffEq,) solve(
            prob_ode, BridgeR3(),
            dt = dt
        )
        @test length(JET.get_reports(rep)) == 0

        rep = @report_opt target_modules = (BridgeDiffEq,) solve(
            prob_ode, BridgeBS3(),
            dt = dt
        )
        @test length(JET.get_reports(rep)) == 0
    end

    @testset "Scalar SDE type stability" begin
        rep = @report_opt target_modules = (BridgeDiffEq,) solve(
            prob_sde, BridgeEuler(),
            dt = dt
        )
        @test length(JET.get_reports(rep)) == 0

        rep = @report_opt target_modules = (BridgeDiffEq,) solve(
            prob_sde, BridgeHeun(),
            dt = dt
        )
        @test length(JET.get_reports(rep)) == 0

        rep = @report_opt target_modules = (BridgeDiffEq,) solve(
            prob_sde, BridgeSRK(),
            dt = dt
        )
        @test length(JET.get_reports(rep)) == 0
    end

    @testset "Vector type stability" begin
        u0_vec = @SVector [2.0, 3.0]
        prob_ode_vec = ODEProblem(f, u0_vec, tspan)
        prob_sde_vec = SDEProblem(f, g, u0_vec, tspan)

        rep = @report_opt target_modules = (BridgeDiffEq,) solve(
            prob_ode_vec, BridgeR3(),
            dt = dt
        )
        @test length(JET.get_reports(rep)) == 0

        rep = @report_opt target_modules = (BridgeDiffEq,) solve(
            prob_ode_vec, BridgeBS3(),
            dt = dt
        )
        @test length(JET.get_reports(rep)) == 0

        rep = @report_opt target_modules = (BridgeDiffEq,) solve(
            prob_sde_vec,
            BridgeEuler(), dt = dt
        )
        @test length(JET.get_reports(rep)) == 0

        rep = @report_opt target_modules = (BridgeDiffEq,) solve(
            prob_sde_vec,
            BridgeHeun(), dt = dt
        )
        @test length(JET.get_reports(rep)) == 0
    end
end
