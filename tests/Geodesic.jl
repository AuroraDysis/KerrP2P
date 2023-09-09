# We set E = 1

using UnPack
using DifferentialEquations

mutable struct KerrCache{T<:AbstractFloat}
    M::T
    a::T
    λ::T
    rh::T # horizon radius
    ro::T # observer radius
end

function stop_condition_r(u, t, integrator)
    @unpack rh, ro = integrator.p
    u[2] <= 1.01 * rh || u[2] >= ro
end

function calc_ray(M::T, a::T, pos::Array{T, 1}, λ::T, η::T, ν_r::Int, ν_θ::Int, ro::T, atol::T, rtol::T) where T <: AbstractFloat
    t, r, θ, ϕ = pos

    sinθ, cosθ = sincos(θ)

    # metric components
    Σ = r^2 + a^2 * cosθ^2
    Δ = r^2 - 2 * M * r + a^2

    gtt = -(1 - 2 * M * r / Σ)
    grr = Σ / Δ
    gθθ = Σ
    gϕϕ = (r^2 + a^2 + 2 * M * r * a^2 * sinθ^2 / Σ) * sinθ^2
    gϕt = -2 * M * r * a / Σ

    R = (r^2 + a^2 - a * λ)^2 - Δ * (η + (λ - a)^2)
    Theta = η + (a * cosθ)^2 - (λ * cosθ / sinθ)^2

    tdot = ((r^2 + a^2) / Δ * (r^2 + a^2 - a * λ) + a * (λ - a * sinθ^2)) / Σ
    rdot = ν_r * sqrt(R) / Σ
    θdot = ν_θ * sqrt(Theta) / Σ
    ϕdot = (a / Δ * (r^2 + a^2 - a * λ) + λ / sinθ^2 - a) / Σ

    # null geodesic condition
    cond = rdot^2 * grr + tdot^2 * gtt + θdot^2 * gθθ + 2 * tdot * ϕdot * gϕt + ϕdot^2 * gϕϕ
    @assert isapprox(cond, 0.0, atol=10 * eps(T)) "null geodesic condition is not satisfied"

    pr = grr * rdot
    pθ = gθθ * θdot

    rh = M + sqrt(M^2 - a^2)
    cache = KerrCache{T}(M, a, λ, rh, ro)

    u0 = [t, r, θ, ϕ, pr, pθ]
    prob = ODEProblem(eom_back!, u0, (zero(T), typemax(T)), cache)
    cb = ContinuousCallback(stop_condition_r, igt -> terminate!(igt))
    # only save the last step
    sol = solve(prob, Vern8(), reltol=1e-10, abstol=1e-10, save_start = false, save_everystep = false, callback=cb)

    print(sol[end])
end

function eom_back!(du::Array{T, 1}, u::Array{T, 1}, p::KerrCache{T}, t::T) where T <: AbstractFloat
    @unpack M, a, λ = p
    t, r, θ, ϕ, pr, pθ = u

    @inbounds begin
        du[1] = 4 * r .* (-λ .* a + a .^ 2 + r .^ 2) ./ ((a .^ 2 + r .* (r - 2)) .* (a .^ 2 .* cos(2 * θ) + a .^ 2 + 2 * r .^ 2)) + 1
        du[2] = pr .* (a .^ 2 + r .* (r - 2)) ./ (a .^ 2 .* cos(θ) .^ 2 + r .^ 2)
        du[3] = pθ ./ (a .^ 2 .* cos(θ) .^ 2 + r .^ 2)
        du[4] = (2 * λ .* (a .^ 2 + r .* (r - 2)) .* csc(θ) .^ 2 + 2 * a .* (-λ .* a + 2 * r)) ./ ((a .^ 2 + r .* (r - 2)) .* (a .^ 2 .* cos(2 * θ) + a .^ 2 + 2 * r .^ 2))
        du[5] = (-a .^ 2 .* pr .^ 2 .* (r - 1) .* cos(θ) .^ 2 + r .* (pr .^ 2 .* (a .^ 2 - r) + pθ .^ 2)) ./ (a .^ 2 .* cos(θ) .^ 2 + r .^ 2) .^ 2 + (-6 * λ .^ 2 .* a .^ 4 .* r + 4 * λ .^ 2 .* r .* (a .^ 2 + r .* (r - 2)) .^ 2 .* csc(θ) .^ 2 + 12 * λ .* a .^ 2 .* r .^ 2 .* (λ + a) + 2 * a .^ 4 .* (-λ + a) .^ 2 + 2 * a .^ 2 .* (-λ .^ 2 .* a .^ 2 .* r + a .^ 2 .* (-λ + a) .^ 2 + 2 * a .* r .^ 2 .* (λ + a) + r .^ 4 - 4 * r .^ 3) .* cos(2 * θ) - 6 * a .* r .^ 4 .* (-4 * λ + a) + 8 * a .* r .^ 3 .* (-λ .^ 2 .* a - 4 * λ + a) - 4 * r .^ 6) ./ ((a .^ 2 + r .* (r - 2)) .^ 2 .* (a .^ 2 .* cos(2 * θ) + a .^ 2 + 2 * r .^ 2) .^ 2)
        du[6] = λ .^ 2 .* (4 * a .^ 4 .* r .* sin(2 * θ) ./ ((a .^ 2 + r .* (r - 2)) .* (a .^ 2 .* cos(2 * θ) + a .^ 2 + 2 * r .^ 2) .^ 2) + cot(θ) .* csc(θ) .^ 2) ./ (a .^ 2 + r .^ 2) - 2 * a .^ 2 .* (4 * λ .* a .* r + a .^ 4 .* pr .^ 2 + a .^ 2 .* (pθ .^ 2 + 2 * r .* (pr .^ 2 .* (r - 2) - 1)) + r .* (pθ .^ 2 .* (r - 2) + r .* (pr .^ 2 .* (r - 2) .^ 2 - 2 * r))) .* sin(2 * θ) ./ ((a .^ 2 + r .* (r - 2)) .* (a .^ 2 .* cos(2 * θ) + a .^ 2 + 2 * r .^ 2) .^ 2)
    end
end

function main()
    λ = -0.7511316141196980351821846354387491739005901369953212373679122913636
    η = 26.5724289697094692725198762793446996667830541429568972667550981404403
    ro = 1000.0
    atol = 1e-10
    rtol = 1e-10
    calc_ray(1.0, 0.8, [0.0, 10.0, pi / 2, 0.0], λ, η, -1, -1, ro, atol, rtol)
end

main()
