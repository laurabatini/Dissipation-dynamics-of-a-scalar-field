include("functions.jl")
using DelimitedFiles


function Flow_finite_V(;N,d,Λ,Rmax,mass, tfin,abstol,reltol,method, L)
    
    line=OriginInterval(N,Rmax)

    prova=SemiDiscrete(OriginInterval(N,Rmax),FlowEquation(d,Λ,L)) # finite V

    x,flowequation=prova.spatial_discretization,prova.equation_of_motion
    dx=Δx(x)
    U=initial_conditions(mass,1.0,line) # U=initial_condition(mass,1.0,line)
    tspan = (0,tfin)
    prob = ODEProblem(prova,U,tspan,0)
    sol = solve(prob,method,save_everystep=false,save_start=false,save_end=true,abstol = abstol, reltol = reltol)
    sol=sol[1]
    return sol[1, :], sol[2,:] # Uprime and logX 

end

function Flow_finite_V(N, d, Λ, Rmax, mass, tfin, abstol, reltol, method, L)
    line = OriginInterval(N, Rmax)
    prova = SemiDiscrete(OriginInterval(N, Rmax), FlowEquation(d, Λ, L))
    dx = Δx(line)
    U = initial_conditions(mass, 1.0, line)
    tspan = (0, tfin)
    prob = ODEProblem(prova, U, tspan, 0)
    sol = solve(prob, method, save_everystep=false, save_start=false, save_end=true, abstol=abstol, reltol=reltol)
    return sol.u[1, :], sol.u[2, :]
end


function Flow_infinite_V(;N=200,d=3,Λ=10,Rmax=2.0,mass=.8, tfin=10,abstol=1e-9,reltol=1e-10,method=FBDF())
    
    line=OriginInterval(N,Rmax)
    prova=SemiDiscrete(line,FlowEquation(d,Λ)) 
 
    U=initial_conditions(mass,1.0,line) #U=ic_costerm(mass,1.0,line)
    tspan = (0,tfin)
    prob = ODEProblem(prova,U,tspan,0)
    sol = solve(prob,method,
    save_everystep=false,save_start=false,save_end=true,abstol = abstol, reltol = reltol)
    result=sol[1]
    return result[1,:], result[2,:] # Uprime and logX 
end
