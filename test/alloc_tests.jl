using BridgeDiffEq
using StaticArrays
using Test

# Allocation regression tests
# These tests ensure allocations don't regress from known baseline values
# The solve function inherently allocates (creating solutions), but we want
# to ensure allocations stay reasonable and don't regress

@testset "Allocation Regression Tests" begin
    # Setup
    f(u, p, t) = u
    g(u, p, t) = u
    dt = 1 // 16
    tspan = (0.0, 1.0)

    @testset "Scalar ODE Allocations" begin
        u0 = 0.5
        prob = ODEProblem(f, u0, tspan)

        # Warm up
        solve(prob, BridgeR3(), dt = dt)
        solve(prob, BridgeBS3(), dt = dt)

        # Measure allocations
        allocs_r3 = @allocated solve(prob, BridgeR3(), dt = dt)
        allocs_bs3 = @allocated solve(prob, BridgeBS3(), dt = dt)

        # These should not exceed reasonable bounds (in bytes)
        # Baseline: ~22KB for BridgeR3, ~27KB for BridgeBS3
        @test allocs_r3 < 50_000  # 50KB upper bound
        @test allocs_bs3 < 50_000  # 50KB upper bound
    end

    @testset "Scalar SDE Allocations" begin
        u0 = 0.5
        prob = SDEProblem(f, g, u0, tspan)

        # Warm up
        solve(prob, BridgeEuler(), dt = dt)
        solve(prob, BridgeHeun(), dt = dt)
        solve(prob, BridgeSRK(), dt = dt)

        # Measure allocations
        allocs_euler = @allocated solve(prob, BridgeEuler(), dt = dt)
        allocs_heun = @allocated solve(prob, BridgeHeun(), dt = dt)
        allocs_srk = @allocated solve(prob, BridgeSRK(), dt = dt)

        # These should not exceed reasonable bounds (in bytes)
        # Baseline: ~17KB for Euler, ~24KB for Heun, ~26KB for SRK
        @test allocs_euler < 50_000  # 50KB upper bound
        @test allocs_heun < 50_000  # 50KB upper bound
        @test allocs_srk < 50_000  # 50KB upper bound
    end

    @testset "Vector ODE Allocations" begin
        u0 = @SVector [2.0, 3.0]
        prob = ODEProblem(f, u0, tspan)

        # Warm up
        solve(prob, BridgeR3(), dt = dt)
        solve(prob, BridgeBS3(), dt = dt)

        # Measure allocations
        allocs_r3 = @allocated solve(prob, BridgeR3(), dt = dt)
        allocs_bs3 = @allocated solve(prob, BridgeBS3(), dt = dt)

        # These should not exceed reasonable bounds (in bytes)
        @test allocs_r3 < 60_000  # 60KB upper bound (slightly higher for vectors)
        @test allocs_bs3 < 60_000  # 60KB upper bound
    end

    @testset "Vector SDE Allocations" begin
        u0 = @SVector [2.0, 3.0]
        prob = SDEProblem(f, g, u0, tspan)

        # Warm up
        solve(prob, BridgeEuler(), dt = dt)
        solve(prob, BridgeHeun(), dt = dt)

        # Measure allocations
        allocs_euler = @allocated solve(prob, BridgeEuler(), dt = dt)
        allocs_heun = @allocated solve(prob, BridgeHeun(), dt = dt)

        # These should not exceed reasonable bounds (in bytes)
        @test allocs_euler < 60_000  # 60KB upper bound
        @test allocs_heun < 60_000  # 60KB upper bound
    end

    @testset "No regression in allocation per timestep" begin
        # Test with larger problem to ensure allocations scale linearly
        f2(u, p, t) = u
        g2(u, p, t) = u
        u0 = 0.5
        dt_fine = 1 // 64
        tspan_long = (0.0, 10.0)

        prob_ode = ODEProblem(f2, u0, tspan_long)
        prob_sde = SDEProblem(f2, g2, u0, tspan_long)

        # Warm up
        solve(prob_ode, BridgeR3(), dt = dt_fine)
        solve(prob_sde, BridgeEuler(), dt = dt_fine)

        # Measure allocations
        allocs_ode = @allocated solve(prob_ode, BridgeR3(), dt = dt_fine)
        allocs_sde = @allocated solve(prob_sde, BridgeEuler(), dt = dt_fine)

        nsteps = Int((tspan_long[2] - tspan_long[1]) / Float64(dt_fine))

        # Check bytes per step is reasonable (should be O(1) per step for fixed-size state)
        bytes_per_step_ode = allocs_ode / nsteps
        bytes_per_step_sde = allocs_sde / nsteps

        # These values should stay small - mostly just the solution storage
        @test bytes_per_step_ode < 100  # Less than 100 bytes per step
        @test bytes_per_step_sde < 150  # Less than 150 bytes per step (SDE has more state)
    end
end
