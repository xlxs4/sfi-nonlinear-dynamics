# A simple implementation of the box counting algorithm,
# also known as Minkowski-Bouligand or Minkowski (fractal) dimension.

# Determines the Capacity Dimension of a 2D trajectory of a dynamical system.

##
using LinearAlgebra
using SparseArrays

function dcap(t::Tuple{Vector{Float64},Vector{Float64}}, ϵ::Float64)::Tuple{Vector{Int64},Vector{Int64}}
    x, z = t
    xmin = minimum(x)
    zmin = minimum(z)

    x_to_ind(x) = ceil(Int64, (x - xmin) / ϵ)
    z_to_ind(z) = ceil(Int64, (z - zmin) / ϵ)

    return (x_to_ind.(x), z_to_ind.(z))
end


function N(t::Tuple{Vector{Float64},Vector{Float64}}, ϵ::Float64)
    x, z = dcap(t, ϵ)
    w = spzeros(Float64, length(x), length(z))

    for (i, j) in (zip(x, z))
        if (i != 0 && j != 0)
            w[i, j] = 1
        end
    end

    return sum(w)
end
##

##
using CairoMakie
using DataFrames
using DelimitedFiles
using GLM

t = readdlm("CapDimData.dat", ',')

x = t[:, 1]
z = t[:, 3]
t = (x, z)

rlen = 100
r = zeros(Float64, rlen)

i = 1
# ϵrange = range(0.01, 0.65, length = rlen)
ϵrange = range(0.28, 0.5, length = rlen)

for ϵ ∈ ϵrange
    r[i] = N(t, ϵ)
    i += 1
end

ϵlog(ϵ) = 1 / log(ϵ) # helper: f(x) = exp(1/x)
ϵrange = ϵlog.(ϵrange)

r = log.(r)

data = DataFrame(X = ϵrange, Y = r)
ols = lm(@formula(Y ~ X), data)
fit = predict(ols)

lines(ϵrange, r, label = "Log/log Scaling Region of Box Counting")
lines!(ϵrange, fit, label = "Ordinary Least Squares Linear Fit")
axislegend(position = :rb)
current_figure()
##