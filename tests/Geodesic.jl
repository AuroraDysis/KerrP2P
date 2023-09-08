# We set E = 1

using UnPack
using DifferentialEquations

mutable struct KerrCache{T<:AbstractFloat}
    M::T
    a::T
    λ::T
    rh::T
end

function calc_initial_data(M::T, a::T, pos::Array{T, 1}) where T <: AbstractFloat
    t, r, θ, ϕ = pos

    r, θ = pos[2], pos[3]
    sinθ, cosθ = sincos(θ)
    Σ = r^2 + a^2 * cosθ^2
    Δ = r^2 - 2 * M * r + a^2

    gtt = -(1 - 2 * M * r / Σ)
    grr = Σ / Δ
    gθθ = Σ
    gϕϕ = (r^2 + a^2 + 2 * M * r * a^2 * sinθ^2 / Σ) * sinθ^2
    gϕt = -2 * M * r * a / Σ

    λ = 0.0
    η = 0.0

    ν_r = -1
    ν_θ = 0

    R = (a^2 - a * λ + r^2)^2 - (λ - a)^2 + η^2
    Theta = η + (a * cosθ)^2 - (λ * cosθ / sinθ)^2

    tdot = ((r^2 + a^2) / Δ * (r^2 + a^2 - a * λ) + a * (λ - a * sinθ^2)) / Σ
    rdot = ν_r * sqrt(R) / Σ
    θdot = ν_θ * sqrt(Theta) / Σ
    ϕdot = (a / Δ * (r^2 + a^2 - a * λ) + λ / sinθ^2 - a) / Σ

    xdot = [tdot, rdot, θdot, ϕdot]

    gdn = zeros(4, 4)
    gdn[1, 1] = gtt
    gdn[2, 2] = grr
    gdn[3, 3] = gθθ
    gdn[4, 4] = gϕϕ
    gdn[4, 1] = gdn[1, 4] = gϕt
    
    println(transpose(xdot) * gdn * xdot)

    pr = grr * rdot
    pθ = gθθ * θdot

    rh = M + sqrt(M^2 - a^2)
    KerrCache{T}(T(M), T(a), T(λ), T(rh))

    [t, r, θ, ϕ, tdot, rdot, θdot, ϕdot]
end

function eom_back!(du::Array{T, 1}, u::Array{T, 1}, p::KerrCache{T}, t::T) where T <: AbstractFloat
    @unpack M, a, λ = p
    t, r, θ, ϕ, pr, pθ = u

    @inbounds begin
        du[1] = 1 + (4 * r * (a^2 - a * λ + r^2)) / ((a^2 + (-2 + r) * r) * (a^2 + 2 * r^2 + a^2 * cos(2 * θ)))
        du[2] = (pr * (a^2 + (-2 + r) * r)) / (r^2 + a^2 * cos(θ)^2)
        du[3] = pθ / (r^2 + a^2 * cos(θ)^2)
        du[4] = (2 * (a * (-a * λ + 2 * r) + λ * (a^2 + (-2 + r) * r) * csc(θ)^2)) / ((a^2 + (-2 + r) * r) * (a^2 + 2 * r^2 + a^2 * cos(2 * θ)))
        du[5] = ((pθ^2 + pr^2 * (a^2 - r)) * r - a^2 * pr^2 * (-1 + r) * cos(θ)^2) / (r^2 + a^2 * cos(θ)^2)^2
        du[6] = (-((2 * a^2 * (a^4 * pr^2 + 4 * a * λ * r + a^2 * (pθ^2 + 2 * (-1 + pr^2 * (-2 + r)) * r) + r * (pθ^2 * (-2 + r) + (pr^2 * (-2 + r)^2 - 2 * r) * r)) * sin(2 * θ))) / ((a^2 + (-2 + r) * r) * (a^2 + 2 * r^2 + a^2 * cos(2 * θ))^2)) + (λ^2 * (cot(θ) * csc(θ)^2 + (4 * a^4 * r * sin(2 * θ)) / ((a^2 + (-2 + r) * r) * (a^2 + 2 * r^2 + a^2 * cos(2 * θ))^2))) / (a^2 + r^2)
    end
   end

calc_initial_data(1.0, 0.8, [0.0, 10.0, pi / 2, 0.0])
