using BridgeDiffEq, StaticArrays
using Test

const GROUP = get(ENV, "GROUP", "all")

@testset "BridgeDiffEq.jl" begin
    if GROUP == "all" || GROUP == "core"
        @testset "Core Functionality" begin
            α = 1
            β = 1
            u0 = 1 / 2
            f(u, p, t) = α * u
            g(u, p, t) = β * u
            dt = 1 // 2^(4)
            tspan = (0.0, 1.0)

            @testset "Scalar ODE" begin
                prob = ODEProblem(f, u0, (0.0, 1.0))
                sol = solve(prob, BridgeR3(), dt = dt)
                @test length(sol.t) == 17
                sol = solve(prob, BridgeBS3(), dt = dt)
                @test length(sol.t) == 17
            end

            @testset "Scalar SDE" begin
                prob = SDEProblem(f, g, u0, (0.0, 1.0))
                sol = solve(prob, BridgeEuler(), dt = dt)
                @test length(sol.t) == 17
                sol = solve(prob, BridgeHeun(), dt = dt)
                @test length(sol.t) == 17
                sol = solve(prob, BridgeSRK(), dt = dt)
                @test length(sol.t) == 17
            end

            @testset "Vector ODE" begin
                u0_vec = @SVector [2.0, 3.0]
                prob = ODEProblem(f, u0_vec, (0.0, 1.0))
                sol = solve(prob, BridgeR3(), dt = dt)
                @test length(sol.t) == 17
                sol = solve(prob, BridgeBS3(), dt = dt)
                @test length(sol.t) == 17
            end

            @testset "Vector SDE" begin
                u0_vec = @SVector [2.0, 3.0]
                prob = SDEProblem(f, g, u0_vec, (0.0, 1.0))
                sol = solve(prob, BridgeEuler(), dt = dt)
                @test length(sol.t) == 17
                sol = solve(prob, BridgeHeun(), dt = dt)
                @test length(sol.t) == 17
            end
        end
    end

    if GROUP == "all" || GROUP == "nopre"
        @testset "Allocation Tests" begin
            include("alloc_tests.jl")
        end
    end

    if GROUP == "all"
        @testset "Explicit Imports" begin
            include("explicit_imports.jl")
        end
    end

    if GROUP == "all" || GROUP == "jet"
        @testset "JET Static Analysis" begin
            include("jet_tests.jl")
        end
    end
end
