# Play around with the Lotka-Volterra system.
# No spatial diffusion here.

##
using DifferentialEquations
using Plots
##

##
function lotkavolterra!(du, u, p, t)
    α, β, γ, δ = p

    du[1] = α * u[1] - β * u[1] * u[2]
    du[2] = -γ * u[2] + δ * β * u[1] * u[2]
end
##

##
p = (1.0, 2.0, 1.5, 1.25)
u0 = [0.4; 0.2]
tspan = (0.0, 10.0)

prob = ODEProblem(lotkavolterra!, u0, tspan, p)
sol = solve(prob)
##

##
plot(sol, vars = (1, 2), title = "Lotka-Volterra Phase Space")
##

##
plot(sol, title = "Lotka-Volterra ODE")
##