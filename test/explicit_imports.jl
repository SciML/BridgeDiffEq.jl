using ExplicitImports
using BridgeDiffEq
using Test

@testset "ExplicitImports" begin
    @test check_no_implicit_imports(BridgeDiffEq) === nothing
    @test check_no_stale_explicit_imports(BridgeDiffEq) === nothing
end
