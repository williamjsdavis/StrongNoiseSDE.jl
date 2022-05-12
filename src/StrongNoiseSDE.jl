module StrongNoiseSDE

export get_P, f̄_σ_posterior

Δ(d20,d21,d22) = 4*d20*d22 - d21^2
a(d11,d22) = d11/(2*d22) - 1
b(d10,d11,d21,d22) = d10 - d21*d11/(2*d22)
h₀(x,d20,d21,d22) = 2*(atan((2*d22*x + d21)/sqrt(Δ(d20,d21,d22))) + π/2)/sqrt(Δ(d20,d21,d22))

D2(x,d20,d21,d22) = d20 + d21*x + d22*x^2
P_g(x,d10,d11,d20,d21,d22) = D2(x,d20,d21,d22)^a(d11,d22)*exp(b(d10,d11,d21,d22)*h₀(x,d20,d21,d22))
#P_g(x,d10,d11,d20,d21,d22) = begin
#    println((d10,d11,d20,d21,d22))
#    return D2(x,d20,d21,d22)^a(d11,d22)*exp(b(d10,d11,d21,d22)*h₀(x,d20,d21,d22))
#end

function integrate(x, y)
    retval = (x[2] - x[1]) * (y[1] + y[2])
    for i in 2:(length(y) - 1)
        retval += (x[i+1] - x[i]) * (y[i] + y[i+1])
    end
    return 0.5 * retval
end

# Conditional probabilities
get_N(x_domain,P) = 1/integrate(x_domain,P.(x_domain))
function get_P(d10,d11,d20,d21,d22,x_domain)
    P_set(x) = P_g(x,d10,d11,d20,d21,d22)
    N = get_N(x_domain,P_set)
    P_norm(x) = N*P_set(x)
    return P_norm
end

# Conditional probabilities
function get_f̄_σ_posterior(P_norm,x_domain)
    f_σ_likelihood(y,x,σ) = exp(-(y-x)^2/(2*σ^2))/(σ*sqrt(2*π))
    integrand(y,x,σ) = f_σ_likelihood(y,x,σ)*P_norm(x)
    denom(y,σ,x_domain) = integrate(x_domain, integrand.(y,x_domain,σ))
    f̄_σ_posterior(x,y,σ,x_domain) = integrand(y,x,σ)/denom(y,σ,x_domain)
    return f̄_σ_posterior
end

# Variable data
export calculate_m₁, calculate_m₂, calculate_γ₁, calculate_γ₂

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
