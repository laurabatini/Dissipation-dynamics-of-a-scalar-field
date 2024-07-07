include("discretization.jl")
include("functions_new.jl")
using QuadGK
using Integrals

N= 400
Rmax=2
d=3
Λ=10
mass=-1
lam=4 
H=-10
L=InfiniteVolume()
grid=SymmetricInterval(N,Rmax)
ϕ =grid.x
U0 = @. H+ mass*ϕ + lam*ϕ^3
plot(ϕ,U0)
u0=transpose(stack(((U0,zeros(length(U0))))))

p=(grid,FlowEquation(d, Λ, L))
tspan = (0,400)
prob = ODEProblem(not_symmetric,u0,tspan,p)

sol = solve(prob)

using Plots
plot(sol.u[1][1,:])
plot(sol.u[end][1,:])

p1=plot(sol.u[2][2,:])
for t  in 3:2:length(sol.t)
   plot!(p1,sol.u[t][2,:])
end
p1

p1=plot(sol.u[2][1,:])
for t  in 3:20:length(sol.t)
   plot!(p1,sol.u[t][1,:])
end
p1
savefig(p1,"fancy_Stuff.pdf")
plot(sol.u[end][1,:])

function integrale_potential(sol,disc,idx,t)
    
    y=sol(t)[1,1:idx]
    x= disc.x[1:idx]
    
    solve(SampledIntegralProblem(y, x), TrapezoidalRule()).u
end

integrale(sol,grid,2,1)

integrale.(Ref(sol),Ref(grid),2:100,Ref(1))

p1=plot(integrale.(Ref(sol),Ref(grid),2:400,Ref(1)))
for t  in 1:3:100
 plot!(p1,integrale.(Ref(sol),Ref(grid),2:400,Ref(t)))
 #xlims!(25,75)
end     
p1