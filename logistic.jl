# A more basic introduction to the logistic map.

##
using Plots, Statistics
##

##
function logisticmap(n, x₀, r, k = 0)
    n = n + 1
    xₛ = zeros(n)

    for i in 1:k
        t = x₀
        x₀ = r * t * (1 - t)
    end

    xₛ[1] = x₀
    for i in 2:n
        xₛ[i] = r * xₛ[i-1] * (1 - xₛ[i-1])
    end
    return xₛ
end
##

##
function experiment()
    r = 3.72
    n = 5000

    x₀ = 0.2
    x̂₀ = 0.200001

    xₙ = logisticmap(n, x₀, r)
    x̂ₙ = logisticmap(n, x̂₀, r)

    return xₙ, x̂ₙ, n
end
##

##
result = experiment()
xₙ, x̂ₙ, n = result

plot(abs.(xₙ - x̂ₙ), seriestype = :scatter, markersize = 3.3, title = "|xₙ - x̂ₙ|")

mean(abs.(xₙ - x̂ₙ))
##

##
function bifurcation(x₀, rmin, rmax, rstep, n, k)
    range = rmin:rstep:rmax
    iter = length(range)
    res = zeros(iter, n)

    currentindex = 1
    for r = range
        t = logisticmap(n, x₀, r, k)
        for i in 1:n
            res[currentindex, i] = t[i]
        end
        currentindex += 1
    end

    return res
end
##

##
x₀ = 0.5
rmin = 3.84
rmax = 3.8571
rstep = 0.0001
k = 5000
n = 2000
res = bifurcation(x₀, rmin, rmax, rstep, k, n)
plot(rmin:rstep:rmax,
    res,
    seriestype = :scatter,
    markercolor = :blue,
    markerstrokecolor = :match,
    markersize = 1.1,
    markerstrokewidth = 0,
    markeralpha = 0.3,
    legend = false,
    title = "Logistic map bifurcation diagram")
# xlims!(2.9, 4)
##