Base.@kwdef mutable struct Parameters
    σ::Float64 = 0.19
    Λ::Float64 = 1.0
    m0::Float64 = 0.0
end

Base.@kwdef struct Model
    pot::Function
    param::Parameters
    n::Int64
end

Base.@kwdef struct ResultGap
    Tlist::Vector{Float64}
    massgap::Vector{Float64}
end