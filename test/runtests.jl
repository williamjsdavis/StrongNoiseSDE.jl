using StrongNoiseSDE
using JLD2
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

# End to end test
@testset "End to end tests" begin
    N = 10
    γ̂₁ = zeros(N)
    intercepts = zeros(N)
    m̂₁ = zeros(N)
    m̂₂ = zeros(N)
    σ̂²ᵧ₁ = ones(N)
    σ̂²ᵧ₂ = ones(N)
    σ̂²m₁ = ones(N)
    σ̂²m₂ = ones(N)
    x_domain = LinRange(-10,10,100)
    y_centers = LinRange(-6,6,N)
    F_pac = StrongNoiseSDE.get_F(γ̂₁,intercepts,m̂₁,m̂₂,σ̂²ᵧ₁,σ̂²ᵧ₂,σ̂²m₁,σ̂²m₂,x_domain,y_centers)
    
    p_test = [1.,-1.,1.,-1.,1.,0.2]
    F_test = F_pac(p_test...)
    @test isreal(F_test)
    @test isfinite(F_test)
end

@testset "Precomputed predictions: multiplicative + noise" begin
    # Precomputed
    @load "test-predictions-1.jld2" γ₁ γ₂ m₁ m₂ y_centers p_true
    @test size(p_true) == (6,)
    @test all(p_true .== [0,-1,1,0,1,1])
    N = size(y_centers)
    @test size(γ₁) == N
    @test size(γ₂) == N
    @test size(m₁) == N
    @test size(m₂) == N

    # Calculated
    x_domain = LinRange(-10,10,500)
    preditions = StrongNoiseSDE.get_predictions(x_domain,y_centers)
    γ₁_calc, γ₂_calc, m₁_calc, m₂_calc = preditions(p_true...)

    # Compare
    #NOTE: set some approximate equality in the future?
    @test all(γ₁ .== γ₁_calc)
    @test all(γ₂ .== γ₂_calc)
    @test all(m₁ .== m₁_calc)
    @test all(m₂ .== m₂_calc)
end

