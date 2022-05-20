module StrongNoiseSDE

include("utils.jl")

# Probabilities

export get_P, f̄_σ_posterior

include("probabilities.jl")


# Observations

export calculate_m₁, calculate_m₂, calculate_γ₁, calculate_γ₂

include("observations.jl")

# Predictions

export get_predictions

include("predictions.jl")

# Objective function
export get_F
function get_F(γ̂₁,intercepts,m̂₁,m̂₂,σ̂²ᵧ₁,σ̂²ᵧ₂,σ̂²m₁,σ̂²m₂,x_domain,y_centers)
    
    M = length(y_centers)
    P_set(d10,d11,d20,d21,d22) = get_P(d10,d11,d20,d21,d22,x_domain)
    f̄_σ_set(d10,d11,d20,d21,d22) = get_f̄_σ_posterior(P_set(d10,d11,d20,d21,d22),x_domain)
    
    γ̂₂(σ) = intercepts .- σ^2
    γ₁_f(d10,d11,d20,d21,d22,σ) = calculate_γ₁(y_centers,f̄_σ_set(d10,d11,d20,d21,d22),
        d10,d11,d20,d21,d22,σ,x_domain)
    γ₂_f(d10,d11,d20,d21,d22,σ) = calculate_γ₂(y_centers,f̄_σ_set(d10,d11,d20,d21,d22),
        d10,d11,d20,d21,d22,σ,x_domain)
    m₁_f(d10,d11,d20,d21,d22,γ₁,γ₂) = calculate_m₁(y_centers,d10,d11,d20,d21,d22,γ₁,γ₂)
    m₂_f(d10,d11,d20,d21,d22,γ₁,γ₂) = calculate_m₂(y_centers,d10,d11,d20,d21,d22,γ₁,γ₂)
    
    function F(d10,d11,d20,d21,d22,σ)
        # Data (+globals)
        γ̂₂_setσ = γ̂₂(σ)
        
        # Predictions
        γ₁ = γ₁_f(d10,d11,d20,d21,d22,σ)
        γ₂ = γ₂_f(d10,d11,d20,d21,d22,σ)
        m₁ = m₁_f(d10,d11,d20,d21,d22,γ₁,γ₂)
        m₂ = m₂_f(d10,d11,d20,d21,d22,γ₁,γ₂)
        
        
        # Loss
        loss_sum = 0.
        for (i,y) in enumerate(y_centers)
            loss_sum += (γ̂₁[i] - γ₁[i])^2/σ̂²ᵧ₁[i] + 
                         (γ̂₂_setσ[i] - γ₂[i])^2/σ̂²ᵧ₂[i] + 
                         (m̂₁[i] - m₁[i])^2/σ̂²m₁[i] + 
                         (m̂₂[i] - m₂[i])^2/σ̂²m₂[i]
        end
        return loss_sum/M
    end
    return F
end

export get_F_visualise
#import Plots: plot, plot!, scatter, scatter!, grid
function plot_functions(x_points,m̂₁,m̂₂,γ̂₁,intercepts,σ̂²m₁,σ̂²m₂,σ̂²ᵧ₁,σ̂²ᵧ₂)
    p1 = scatter(x_points, m̂₁, yerror=sqrt.(σ̂²m₁),
        label="Data", title="m̂₁")
    p2 = scatter(x_points, m̂₂, yerror=sqrt.(σ̂²m₂),
        label="Data", title="m̂₂")
    p3 = scatter(x_points, γ̂₁, yerror=sqrt.(σ̂²ᵧ₁),
        label="Data", title="γ̂₁")
    p4 = scatter(x_points, intercepts, yerror=sqrt.(σ̂²ᵧ₂),
        label="Data", title="γ̂₂+σ²")
    return p1, p2, p3, p4
end
function get_F_visualise(γ̂₁,intercepts,m̂₁,m̂₂,σ̂²ᵧ₁,σ̂²ᵧ₂,σ̂²m₁,σ̂²m₂,x_domain,y_centers)
    
    M = length(y_centers)
    P_set(d10,d11,d20,d21,d22) = get_P(d10,d11,d20,d21,d22,x_domain)
    f̄_σ_set(d10,d11,d20,d21,d22) = get_f̄_σ_posterior(P_set(d10,d11,d20,d21,d22),x_domain)
    
    γ̂₂(σ) = intercepts .- σ^2
    γ₁_f(d10,d11,d20,d21,d22,σ) = calculate_γ₁(y_centers,f̄_σ_set(d10,d11,d20,d21,d22),
        d10,d11,d20,d21,d22,σ,x_domain)
    γ₂_f(d10,d11,d20,d21,d22,σ) = calculate_γ₂(y_centers,f̄_σ_set(d10,d11,d20,d21,d22),
        d10,d11,d20,d21,d22,σ,x_domain)
    m₁_f(d10,d11,d20,d21,d22,γ₁,γ₂) = calculate_m₁(y_centers,d10,d11,d20,d21,d22,γ₁,γ₂)
    m₂_f(d10,d11,d20,d21,d22,γ₁,γ₂) = calculate_m₂(y_centers,d10,d11,d20,d21,d22,γ₁,γ₂)
    m₁_f_eq10(d10,d11,d20,d21,d22,σ) = calculate_m₁_eq10(y_centers,f̄_σ_set(d10,d11,d20,d21,d22),
        d10,d11,d20,d21,d22,σ,x_domain)
    m₂_f_eq10(d10,d11,d20,d21,d22,σ) = calculate_m₂_eq10(y_centers,f̄_σ_set(d10,d11,d20,d21,d22),
        d10,d11,d20,d21,d22,σ,x_domain)
    
    
    function F(d10,d11,d20,d21,d22,σ)
        # Data (+globals)
        γ̂₂_setσ = γ̂₂(σ)
        
        # Predictions
        γ₁ = γ₁_f(d10,d11,d20,d21,d22,σ)
        γ₂ = γ₂_f(d10,d11,d20,d21,d22,σ)
        m₁ = m₁_f(d10,d11,d20,d21,d22,γ₁,γ₂)
        m₂ = m₂_f(d10,d11,d20,d21,d22,γ₁,γ₂)
        m₁_eq10 = m₁_f_eq10(d10,d11,d20,d21,d22,σ)
        m₂_eq10 = m₂_f_eq10(d10,d11,d20,d21,d22,σ)
        
        p1, p2, p3, p4 = plot_functions(y_centers,m̂₁,m̂₂,γ̂₁,intercepts,σ̂²m₁,σ̂²m₂,σ̂²ᵧ₁,σ̂²ᵧ₂)
        scatter!(p1, y_centers, m₁, label="Pred.")
        scatter!(p1, y_centers, m₁_eq10, label="Pred. eq10")
        scatter!(p2, y_centers, m₂, label="Pred.")
        scatter!(p2, y_centers, m₂_eq10, label="Pred. eq10")
        scatter!(p3, y_centers, γ₁, label="Pred.")
        scatter!(p4, y_centers, γ₂ .+ σ^2, label="Pred.")
        plot(p1, p2, p3, p4, layout=grid(2,2), size=(800,500)) |> display
        
        # Loss
        loss_sum = 0.
        for (i,y) in enumerate(y_centers)
            loss_sum += (γ̂₁[i] - γ₁[i])^2/σ̂²ᵧ₁[i] + 
                         (γ̂₂_setσ[i] - γ₂[i])^2/σ̂²ᵧ₂[i] + 
                         (m̂₁[i] - m₁[i])^2/σ̂²m₁[i] + 
                         (m̂₂[i] - m₂[i])^2/σ̂²m₂[i]
        end
        return loss_sum/M
    end
    return F
end
end
