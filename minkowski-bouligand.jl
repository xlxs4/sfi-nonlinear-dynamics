# A simple implementation of the box counting algorithm,
# also known as Minkowski-Bouligand or Minkowski (fractal) dimension.

# Determines the Capacity Dimension of a 2D trajectory of a dynamical system.

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