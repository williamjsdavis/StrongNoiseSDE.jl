# Predictions
function get_predictions(x_domain,y_centers)
    
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
    
    function predictions(d10,d11,d20,d21,d22,σ)
        γ₁ = γ₁_f(d10,d11,d20,d21,d22,σ)
        γ₂ = γ₂_f(d10,d11,d20,d21,d22,σ)
        m₁ = m₁_f(d10,d11,d20,d21,d22,γ₁,γ₂)
        m₂ = m₂_f(d10,d11,d20,d21,d22,γ₁,γ₂)
        return γ₁, γ₂, m₁, m₂
    end
    return predictions
end
