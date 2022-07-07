# Partition the logistic map.
# Some symbolic dynamics fun.

##
using Plots

plot(-x * -x^3)

function iteratelogistic(x0::Float64, r::Float64, n::Int)
    f(x) = r * x * (1 - x)
    orbit = zeros(n)
    orbit[1] = x0
    for i ∈ 2:n
        orbit[i] = f(orbit[i-1])
    end
    orbit
end

function generatepartition(orbit::Vector{Float64})
    [x < 0.5 ? 'L' : 'R' for x ∈ orbit]
end

x0 = 0.2
y0 = 0.4
r = 4.0
n = 4000

orbitx = iteratelogistic(x0, r, n)
orbity = iteratelogistic(y0, r, n)
Δorbit = orbitx - orbity

t = 1:n
plot(t, orbitx, title = "Different initial condition orbits", label = "x", lw = 3)
plot!(t, orbity, label = "y", lw = 3)

plot(t, Δorbit, title = "Difference between orbits", label = "difference", lw = 3)

function countpartition(partition::Vector{Char})
    lcount = 0
    rcount = 0
    for l ∈ partition
        if l == 'L'
            lcount += 1
        else
            rcount += 1
        end
    end
    (lcount, rcount)
end

println(countpartition(generatepartition(orbitx)))
##