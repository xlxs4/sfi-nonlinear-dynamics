# Play around with the Hénon map.

##
using Plots
using SIMD
##

##
function henon!(α, β, xs, ys, x₀, y₀, n)
    f(x, y) = 1 - α * x^2 + y
    g(y, x) = β * x

    xs[:, 1] .= x₀
    ys[:, 1] .= y₀

    @inbounds @simd for j in 1:n-1
        xs[j+1] = f(xs[j], ys[j])
        ys[j+1] = g(ys[j], xs[j])
    end
    xs, ys
end
##

##
α = 1.4
β = 0.3

x₀ = 0.2
y₀ = 0.2

n = 40

xs = zeros(n)
ys = zeros(n)

henon!(α, β, xs, ys, x₀, y₀, n)
##

##
plot(1:n, xs)
##

##
plot(1:n, ys)
##