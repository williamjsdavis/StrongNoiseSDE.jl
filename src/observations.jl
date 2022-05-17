# Observed data
# NOTE: Data structure is just arrays at the moment. Might be worthwhile to make a struct in the future

function calculate_m₁(y_centers,d10,d11,d20,d21,d22,γ₁,γ₂)
    m₁ = zeros(length(y_centers))
    for (i,y) in enumerate(y_centers)
        m₁[i] = d10 + d11*(y + γ₁[i])
    end
    return m₁
end
function calculate_m₂(y_centers,d10,d11,d20,d21,d22,γ₁,γ₂)
    m₂ = zeros(length(y_centers))
    for (i,y) in enumerate(y_centers)
        m₂[i] = 2*(
            γ₁[i]*d10 + 
            (γ₂[i] + y*γ₁[i])*d11 + 
            d20 + 
            (γ₁[i] + y)*d21 + 
            (2*y*γ₁[i] + γ₂[i] + y^2)*d22
        )
    end
    return m₂
end
function calculate_m₁_eq10(y_centers,f̄_σ_posterior,d10,d11,d20,d21,d22,σ,x_domain)
    m₁ = zeros(length(y_centers))
    for (i,y) in enumerate(y_centers)
        f₁_set(x) = (d10+d11*x)*f̄_σ_posterior(x,y,σ,x_domain)
        m₁[i] = integrate(x_domain, f₁_set.(x_domain))
    end
    return m₁
end
function calculate_m₂_eq10(y_centers,f̄_σ_posterior,d10,d11,d20,d21,d22,σ,x_domain)
    m₂ = zeros(length(y_centers))
    for (i,y) in enumerate(y_centers)
        f₂_set(x) = ((x-y)*(d10+d11*x)+(d20+d21*x+d22*x^2))*f̄_σ_posterior(x,y,σ,x_domain)
        m₂[i] = 2*integrate(x_domain, f₂_set.(x_domain))
    end
    return m₂
end
function calculate_γ₁(y_centers,f̄_σ_posterior,d10,d11,d20,d21,d22,σ,x_domain)
    γ₁ = zeros(length(y_centers))
    for (i,y) in enumerate(y_centers)
        f₁_set(x) = (x-y)*f̄_σ_posterior(x,y,σ,x_domain)
        γ₁[i] = integrate(x_domain, f₁_set.(x_domain))
    end
    return γ₁
end
function calculate_γ₂(y_centers,f̄_σ_posterior,d10,d11,d20,d21,d22,σ,x_domain)
    γ₂ = zeros(length(y_centers))
    for (i,y) in enumerate(y_centers)
        f₂_set(x) = (x-y)^2*f̄_σ_posterior(x,y,σ,x_domain)
        γ₂[i] = integrate(x_domain, f₂_set.(x_domain))
    end
    return γ₂
end


