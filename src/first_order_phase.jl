include("discretization.jl")
include("functions.jl")

N= 10
Rmax=2
d=3
Λ=10
mass=-10
lam=4 
H=1
L=InfiniteVolume()
ϕ =SymmetricInterval(N,Rmax).x
U0 = @. H+ mass*ϕ + lam*ϕ^3
p=(SymmetricInterval(N,Rmax),FlowEquation(d, Λ, L))
tspan = (0,40)
prob = ODEProblem(not_symmetric,U0,tspan,p)

sol = solve(prob)

using Plots

plot(sol.u)