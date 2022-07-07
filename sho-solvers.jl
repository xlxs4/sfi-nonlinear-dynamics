# Explore different solvers in a Simple Harmonic Oscillator system.

##
Base.@kwdef mutable struct SHO
    Δt::Float64
    x::Float64
    v::Float64
end

function forwardstep!(s::SHO)
    dx = s.v
    dv = -4s.x

    s.x += s.Δt * dx
    s.v += s.Δt * dv
end

function forwardstep(s::SHO)
    x = s.x
    v = s.v

    dx = v
    dv = -4x

    x += s.Δt * dx
    v += s.Δt * dv

    return (x, v)
end

function backwardstep!(s::SHO, r::Tuple{Float64,Float64})
    x, v = r
    dx = v
    dv = -4x

    s.x += s.Δt * dx
    s.v += s.Δt * dv
end

function backwardstep(s::SHO, r::Tuple{Float64,Float64})
    x, v = r
    dx = v
    dv = -4x

    x = s.x
    v = s.v

    x += s.Δt * dx
    v += s.Δt * dv

    return (x, v)
end

function fw(Δt::Float64, x::Float64, v::Float64, n::Int64)
    s = SHO(Δt, x, v)
    w = Vector{Tuple{Float64,Float64}}(undef, n + 1)
    w[1] = (s.x, s.v)
    for i ∈ 2:n+1
        forwardstep!(s)
        w[i] = (s.x, s.v)
    end

    return w
end

function bw(Δt::Float64, x::Float64, v::Float64, n::Int64)
    s = SHO(Δt, x, v)
    w = Vector{Tuple{Float64,Float64}}(undef, n + 1)

    w[1] = (s.x, s.v)
    for i ∈ 2:n+1
        r = forwardstep(s)
        backwardstep!(s, r)

        w[i] = (s.x, s.v)
    end

    return w
end

function trapezoidal(Δt::Float64, x::Float64, v::Float64, n::Int64)
    s = SHO(Δt, x, v)
    w = Vector{Tuple{Float64,Float64}}(undef, n + 1)

    w[1] = (s.x, s.v)
    for i ∈ 2:n+1
        rₑ = forwardstep(s)
        rᵢ = backwardstep(s, rₑ)

        xₑ, vₑ = rₑ
        xᵢ, vᵢ = rᵢ

        s.x = ((xₑ + xᵢ) / 2)
        s.v = ((vₑ + vᵢ) / 2)

        w[i] = (s.x, s.v)
    end

    return w
end

function experiment(Δt::Float64, n::Int64)
    x = -1.0
    v = -2.0

    wₑ = fw(Δt, x, v, n)
    wᵢ = bw(Δt, x, v, n)
    wₜ = trapezoidal(Δt, x, v, n)

    return (wₑ, wᵢ, wₜ)
end
##

##
using CairoMakie

n = 10
wₑ, wᵢ, wₜ = experiment(0.05, n)

lines(wₑ, label = "Explicit Euler")
lines!(wᵢ, label = "Implicit Euler")
lines!(wₜ, label = "Trapezoidal")

axislegend()
current_figure()
##