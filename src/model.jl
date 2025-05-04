using LinearAlgebra, SparseArrays

mutable struct Parameters
    Δ::Float64
    Γ::Float64 
    Ω::Float64
    
    function Parameters(;Δ::Float64=0., Γ::Float64=1/6, Ω::Float64=1.0)
        new(Δ, Γ, Ω)
    end
end

σₓ = sparse(ComplexF64[0. 1.; 1. 0.])
σ₊ = sparse(ComplexF64[0. 0.; 1. 0.])
σ₋ = sparse(ComplexF64[0. 0.; 0. 1.])

function makeHamiltonian(p::Parameters)
    H = -p.Ω/2. * σₓ - p.Δ*σ₊*σ₋
    return H
end

function makeLindbladian(p::Parameters)
    L = √p.Γ * σ₋
    return L
end

function steadyPe(p::Parameters)
    Pe=abs(p.Ω)^2/4 / (p.Δ^2+p.Γ^2/4+abs(p.Ω)^2/2)
    return Pe
end

export Parameters, makeHamiltonian, makeLindbladian, steadyPe