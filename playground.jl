# Play with time-series embedding, producing the corresponding trajectory in reconstruction space.
# (You can verify the Taken theorem condition)
# This is delay-coordinate embedding.

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
using DynamicalSystems
using DelimitedFiles
using LinearAlgebra

s = vec(readdlm("amplitude.dat"))
r = embed(s, 7, 8)

lines(r[:, 1], r[:, 3], linewidth = 0.6)
##

##
using CairoMakie
using DynamicalSystems

r = embed(x, 3, 18)
lines(r[:, 1], r[:, 3], linewidth = 0.6)
##

##
using DelimitedFiles
using DynamicalSystems

s = vec(readdlm("amplitude.dat"))

t = embed(s, 2, 8)

x, z = t[:,1], t[:,2]
r = N((x, z), 0.5)
##