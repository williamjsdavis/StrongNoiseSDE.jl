using StrongNoiseSDE
using Test

@test 1==1

# Probability tests
@testset "Probability internals" begin
    @test StrongNoiseSDE.a(1,2) == -0.75
    @test StrongNoiseSDE.a(10,0.5) == 9.0
end
