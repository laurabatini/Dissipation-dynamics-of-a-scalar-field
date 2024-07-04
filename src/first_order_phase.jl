include("discretization.jl")
include("functions.jl")
using QuadGK
using Integrals

N= 100
Rmax=4
d=3
Λ=10
mass=-10
lam=4 
H=1
L=InfiniteVolume()
grid=SymmetricInterval(N,Rmax)
ϕ =grid.x
U0 = @. H+ mass*ϕ + lam*ϕ^3
p=(grid,FlowEquation(d, Λ, L))
tspan = (0,200)
prob = ODEProblem(not_symmetric,U0,tspan,p)

sol = solve(prob)

using Plots

plot(sol.u)

function integrale(sol,disc,idx,t)
    
    y=sol(t)[1:idx]
    x= disc.x[1:idx]
    
    solve(SampledIntegralProblem(y, x), TrapezoidalRule()).u
end

integrale(sol,grid,2,1)

integrale.(Ref(sol),Ref(grid),2:100,Ref(1))

p1=plot(integrale.(Ref(sol),Ref(grid),2:100,Ref(1)))
for t  in 1:3:100
 plot!(p1,integrale.(Ref(sol),Ref(grid),2:100,Ref(t)))
 xlims!(25,75)
end     
p1