# For calculating probabilities and conditional probabilities

# Probabilities
# NOTE: Be careful with unicode symbols
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


