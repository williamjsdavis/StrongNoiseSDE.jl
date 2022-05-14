using StrongNoiseSDE
using Test

@test 1==1

# Probability tests
@testset "Probability internals" begin
    @test StrongNoiseSDE.a(1,2) == -0.75
    @test StrongNoiseSDE.a(10,0.5) == 9.0
    @test StrongNoiseSDE.Δ(1.,2.,3.) == 8.0
    @test StrongNoiseSDE.Δ(2.,3.,-4.) == -41.0
    @test StrongNoiseSDE.Δ(2.,4.,2.) == 0.0
    @test StrongNoiseSDE.b(1.,2.,3.,4.) == 0.25
    @test StrongNoiseSDE.b(2.,-3.,4.,5.) == 3.2
    @test StrongNoiseSDE.b(3.,4.,-5.,8.) == 4.25
    @test StrongNoiseSDE.b(-5.,7.,8.,-10.) == -2.2
    @test StrongNoiseSDE.h₀(1.,2.,3.,4.) == 1.1386783618209275
    @test StrongNoiseSDE.h₀(2.,3.,-4.,5.) == 0.828728793408412
    @test StrongNoiseSDE.h₀(-3.,4.,5.,6.) == 0.06299419857612897
    @test StrongNoiseSDE.h₀(0.5,0.1,0.2,0.3) == 18.574186872187475
end
@testset "Probability internals exceptions" begin
    @test_throws DomainError StrongNoiseSDE.h₀(2.,-3.,4.,5.)
    @test_throws DomainError StrongNoiseSDE.h₀(3.,4.,5.,-6.)
end
