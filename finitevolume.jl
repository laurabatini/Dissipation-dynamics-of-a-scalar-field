using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())



include("discretization.jl")
using SpecialFunctions
using ForwardDiff
using LinearAlgebra
using DifferentialEquations
using SparseArrays
using Plots


function ndimcount(dim, R) 
   
   r = floor(Int, R)
    if  dim == 1
       return 1 + 2*r; 
    end
    n = ndimcount( static(dim-1), R)
    for  i in 1:r 
      n += 2*ndimcount( static(dim-1), sqrt( R^2 - i^2))
    end
    return n
end 

using StaticNumbers


Ω(d)= 2* pi^(d/2) / gamma(d/2)  


struct FlowEquation{D,S,O,L}

end

struct InfiniteVolume
end

struct FiniteVolume
end



@inline dimension(::FlowEquation{D,S,O,L}) where {D,S,O,L}= D

@inline scale(::FlowEquation{D,S,O,L}) where {D,S,O,L}= S

@inline SolidAngle(::FlowEquation{D,S,O,InfiniteVolume}) where {D,S,O,L}= O


@inline volume(::FlowEquation{D,S,O,FiniteVolume}) where {D,S,O,L}= O^D

@inline lenght(::FlowEquation{D,S,O,FiniteVolume}) where {D,S,O,L}= O

@inline volume(::FlowEquation{D,S,O,InfiniteVolume} where {D,S,O,} )  = InfiniteVolume()



function FlowEquation(d,Lambda)
    FlowEquation{d,Lambda,Ω(d),InfiniteVolume}()
end

function FlowEquation(d,Lambda,L) 
    if isa(L,Real)
        FlowEquation{d,Lambda,L,FiniteVolume}()  
    elseif isa(L,InfiniteVolume)
     FlowEquation{d,Lambda,Ω(d),InfiniteVolume}()
    else 
        throw(ArgumentError("The lenght of the box is infinity or as specific lenght "))
    end

end

volume(FlowEquation(4,2))

 

function potential_flow(du,t,::FlowEquation{D,S,O,L}) where {D,S,O,L}
    k=S*exp(-t)
    - k^2 *1/ (du +k^2)
end


function Propagator(du,t,flow::FlowEquation{D,S,O,T}) where {D,S,O,T}
    k=S*exp(-t)
    1/ (du +k^2)
end


function Prefactor(t,flow::FlowEquation{D,S,O,InfiniteVolume}) where {D,S,O}
    k=S*exp(-t)
    O/(2 * pi)^D * k^(D)/D 
end

function Prefactor(t,flow::FlowEquation{D,S,L,FiniteVolume}) where {D,S,L}
    k=S*exp(-t)
    ndimcount(static(D),L*k/(2*π))/(L)^D 
end


function alphaflow(propagator, dpropagator,dpropagator2, dalpha, ddalpha,t,::FlowEquation{D,S,O,L}) where {D,S,O,L}
    k=S*exp(-t)
    k^2*0.5*(3*dpropagator^2 + 4*dpropagator2*dalpha + 2*(ddalpha+ dalpha^2)*propagator^2)
end 


struct SemiDiscrete{T,X}
    spatial_discretization::X
    equation_of_motion::T
end

#function SemiDiscrete(x,F)
#    SemiDiscrete{typeof(F),typeof(x)}(x,F)
#end


function (ff::SemiDiscrete)(du,u,p,t)

    x,flowequation=ff.spatial_discretization,ff.equation_of_motion
    dx=Δx(x)
    invdx=1/dx
    invdx2=1/(dx^2)
    prefactor=Prefactor(t, flowequation)
       
    i=firstindex(x)
            #derivataprima= (u[i+1]- u[i-1])/(2*dx)
            #derivataseconda= (u[i+1]+u[i-1] -2* u[i])/(dx^2)
            propagatore_i = - Propagator(  (u[1,i]+ u[1,i])*invdx,t, flowequation ) 
           # propagatore_i = - flowequation(  (u[1,i]+ u[1,i])*invdx,t ) 
            propagatore_i_1 = - Propagator(  (u[1,i+1]- u[1,i])*invdx,t, flowequation ) #derivata sx
           # propagatore_i_1 = - flowequation(  (u[1,i+1]- u[1,i])*invdx,t ) #derivata sx
            dpropagatore_i = ( propagatore_i_1 -  propagatore_i )*invdx
            dpropagatore2_i = ( propagatore_i_1^2 -  propagatore_i^2 )*invdx
            du[1,i]=prefactor*(potential_flow((u[1,i+1]- u[1,i])*invdx,t,flowequation) -potential_flow((u[1,i]+u[1,i])*invdx,t,flowequation))*invdx

            tau=u[2,i]
            dtau=(u[2,i]-u[2,i])*invdx
            ddtau=(u[2,i+1]+u[2,i] - 2* u[2,i])*invdx2

           # print( tauflow_conservative(propagatore_i, dpropagatore_i,dpropagatore2_i,tau,dtau,ddtau,t,flowequation))
            # tauflow_conservative(propagatore_i, dpropagatore_i,dpropagatore2_i,tau,dtau,ddtau,t,flowequation)
            du[2,i]=prefactor*alphaflow(propagatore_i, dpropagatore_i,dpropagatore2_i,dtau,ddtau,t,flowequation)
            # alpha ##du[2,i]=.5*Prefactor(t, flowequation)*alphaflow(propagatore_i, dpropagatore_i,dpropagatore2_i,tau,dtau,ddtau,flowequation)
       for i in interior(x)
            #derivataprima= (u[i+1]- u[i-1])/(2*dx)
            #derivataseconda= (u[i+1]+u[i-1] -2* u[i])/(dx^2)

            propagatore_i = - Propagator(  (u[1,i]- u[1,i-1])*invdx, t, flowequation ) 
            propagatore_i_1 = - Propagator(  (u[1,i+1]- u[1,i])*invdx, t, flowequation ) 
            dpropagatore_i = ( propagatore_i_1 -  propagatore_i )*invdx
            dpropagatore2_i = ( propagatore_i_1^2 - propagatore_i^2 )*invdx
           
            du[1,i]=prefactor*(potential_flow((u[1,i+1]- u[1,i])*invdx,t,flowequation) -potential_flow((u[1,i]- u[1,i-1])*invdx,t,flowequation))*invdx
            #du[1,i]=(tauflow((u[1,i+1]- u[1,i])*invdx,t) -tauflow((u[1,i]- u[1,i-1])*invdx,t))*invdx
            
            dtau=(u[2,i]-u[2,i-1])*invdx
            ddtau=(u[2,i+1]+u[2,i-1] -2* u[2,i])*invdx2
            tau=u[2,i]

            du[2,i]=prefactor*alphaflow(propagatore_i, dpropagatore_i,dpropagatore2_i,dtau,ddtau,t,flowequation)
            ## alpha #du[2,i]=.5*Prefactor(t, flowequation)*alphaflow(propagatore_i, dpropagatore_i,dpropagatore2_i,tau,dtau,ddtau,flowequation)
       end
   
    

       i=lastindex(x)
        ## in the last index we use bg derivative 
        
        du[1,i]= prefactor*(potential_flow((u[1,i]- u[1,i-1])*invdx,t,flowequation) -potential_flow((u[1,i-1]- u[1,i-2])*invdx,t,flowequation))*invdx

        du[2,i]=du[2,i-1]


end



N=100
d=3
Λ=1
Rmax=2

prova=SemiDiscrete(OriginInterval(N,Rmax),FlowEquation(d,Λ,8))

U=zeros(2,length(prova.spatial_discretization.x))

for i in axes(U,2)
   # U[i] = .5*prova.spatial_discretization.x[i]-.1*prova.spatial_discretization.x[i]^3 +.05*prova.spatial_discretization.x[i]^5
   U[1,i] = -.074*prova.spatial_discretization.x[i]+0.5*prova.spatial_discretization.x[i]^3 
   U[2,i]= 1
end

tspan = (0,1000)

prob = ODEProblem(prova,U,tspan,0)



#  Rosenbrock23() Kvaerno5() Karp4() FBDF()

sol = solve(prob,FBDF();progress = true,
progress_steps = 100)



plot_time(sol,tspan,1;Rmax=2,Rmin=0,N=10)
plot_time(sol,tspan,2;Rmax=2,Rmin=0,N=10,funcy=x->exp(-x))

plot_time(sol,tspan,2;Rmax=2,Rmin=0,N=10,funcy=x->x)







function plot_time(sol,tspan,index;Rmax=1,N=10,Rmin=0,funcx=x->x,funcy=x->x)
    nex=filter(x->x<Rmax&&x>Rmin,prova.spatial_discretization.x)
    range=eachindex(nex)
    p= scatter(funcx.(nex),funcy.(sol(tspan[1])[index,range]))
    for i in tspan[1]+(tspan[2]- tspan[1])/N:(tspan[2]- tspan[1])/N:tspan[2]
        scatter!(p,funcx.(nex),funcy.( sol(i)[index,range] ),legend=false)
    end
p
end








function temperaturescan(;N=100,d=3,Λ=1,Rmax=1,mass=-0.2)
    prova=SemiDiscrete(OriginInterval(N,Rmax),FlowEquation(d,Λ))
    U=similar(prova.spatial_discretization.x)
    for i in eachindex(U)
        # U[i] = .5*prova.spatial_discretization.x[i]-.1*prova.spatial_discretization.x[i]^3 +.05*prova.spatial_discretization.x[i]^5
        U[i] = mass*prova.spatial_discretization.x[i]+0.5*prova.spatial_discretization.x[i]^3 
    end
    tspan = (0,25)

    prob = ODEProblem(prova,U,tspan,0)

    sol = solve(prob,Kvaerno5();progress = true,
    progress_steps = 100)
    

    return interpolatezero(sol,prova)
end

function masssacan(;min=0,mfin=-0.3,nmass=50)

    massrange=collect(range(min,mfin,nmass))
    phimin=zeros(nmass)
    for i in ProgressBar(1:nmass) #wrap any iterator
        phimin[i]=temperaturescan(;N=100,d=3,Λ=1,Rmax=1,mass=massrange[i])
    end
    return (phimin, massrange)
end

(phimin, massrange) = masssacan()

scatter(massrange,phimin)

result[1][30][end]
function find_zero_L(x)
    @inbounds for i in 1:length(x)-1
        if x[i]<0 &&x[i+1]>0
        return ((i,i+1),(x[i],x[i+1])  )
        end
        return ((1,1),(0,0)  )
    end
end





function find_zero_R(x)
    @inbounds for i in length(x):-1:2
        if x[i]>0 &&x[i-1]<0
        return ((i-1,i),(x[i-1],x[i])  )
        end
        return ((1,1),(0,0)  )
    end
end

function interpolatezero(sol,prova)
((i_minus,i_plus), (V_minus, V_plus ) )=find_zero_L(sol[end])
    if i_minus==i_plus
        return 0
    end
    phi_minus=prova.spatial_discretization.x[i_minus]
    phi_plus=prova.spatial_discretization.x[i_plus]

    phi  = V_minus /(V_plus-V_minus)*(phi_plus-phi_minus) +phi_minus

    return  phi 
end