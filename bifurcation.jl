# A fast way to plot the bifurcation diagram for the logistic map.
# Remember that some observations can be generalized
# for all 1D maps with a single quadratic maximum.

##
using Plots
using SIMD

function bifur!(rs, x, xs, maxiter)
    x[:, 1] .= xs
    for j in 1:maxiter-1 # The j loop needs to be sequential, so should be the outer loop.
        @inbounds @simd for k in 1:length(rs) # The k loop should be the innermost since it is first on x.
            x[k, j+1] = rs[k] * x[k, j] * (1 - x[k, j])
        end
    end
    return x
end

rs = LinRange(2, 3.9, 6001)
xs = 0.2
maxiter = 1000
x = zeros(length(rs), maxiter)

bifur!(rs, x, xs, maxiter)

plot(rs, x[:, end-100:end], # Avoid transients to get a better plot.
    seriestype = :scatter,
    markercolor = :blue,
    markerstrokecolor = :match,
    markersize = 1.1,
    markerstrokewidth = 0,
    legend = false,
    markeralpha = 0.3)
# xlims!(0.5, 2)
##
