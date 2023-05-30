#using Logging: global_logger
#using TerminalLoggers: TerminalLogger
#global_logger(TerminalLogger())
include("discretization.jl")
using SpecialFunctions
using ForwardDiff
using LinearAlgebra
using Interpolations;
using DifferentialEquations
using SparseArrays
using Plots
using StaticNumbers
using Measurements

using Dierckx
using Roots

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
            
       end
   
    

       i=lastindex(x)
        ## in the last index we use bg derivative 
        
        du[1,i]= prefactor*(potential_flow((u[1,i]- u[1,i-1])*invdx,t,flowequation) -potential_flow((u[1,i-1]- u[1,i-2])*invdx,t,flowequation))*invdx
        du[2,i]=du[2,i-1]


end


# condizioni iniziali 
function initial_condition(mass,lambda,spatial_discretization )
    U=zeros(2,length(spatial_discretization.x))
     
    for i in axes(U,2)
   
        U[1,i] =  -mass*spatial_discretization.x[i]+lambda*spatial_discretization.x[i]^3 
        U[2,i]= 0
    end
    return U
end



function plot_time(sol,tspan,index;Rmax=1,N=10,Rmin=0,funcx=x->x,funcy=x->x)
    nex=filter(x->x<Rmax&&x>Rmin,prova.spatial_discretization.x)
    range=eachindex(nex)
    p= scatter(funcx.(nex),funcy.(sol(tspan[1])[index,range]))
    for i in tspan[1]+(tspan[2]- tspan[1])/N:(tspan[2]- tspan[1])/N:tspan[2]
        scatter!(p,funcx.(nex),funcy.( sol(i)[index,range] ),legend=false)
    end
p
end




function temperaturescan(;N,d=3,Λ=10,Rmax,mass, tfin, L, tre)
        prova=SemiDiscrete(OriginInterval(N,Rmax),FlowEquation(d,Λ,L)) # finite V
       # prova=SemiDiscrete(OriginInterval(N,Rmax),FlowEquation(d,Λ))
       # U = zeros(length(prova.spatial_discretization.x))
        U=zeros(2,length(prova.spatial_discretization.x))
    
        for i in axes(U,2)
       # U[i] = .5*prova.spatial_discretization.x[i]-.1*prova.spatial_discretization.x[i]^3 +.05*prova.spatial_discretization.x[i]^5
         U[1,i] = -mass*prova.spatial_discretization.x[i]+prova.spatial_discretization.x[i]^3 
         U[2,i]= 1
        end
        tspan = (0,tfin)
        prob = ODEProblem(prova,U,tspan,0)
        sol = solve(prob,Rosenbrock23();progress = true,progress_steps = 100)
       
        return interpolatezero(sol,prova, tfin, tre)
end



function temperaturescan_infV(;N=200,d=3,Λ=10,Rmax=2.0,mass=-0.9772, tfin=2,abstol=1e-9,reltol=1e-10,method=FBDF())
    
    line=OriginInterval(N,Rmax)
    prova=SemiDiscrete(line,FlowEquation(d,Λ)) 
 
    U=initial_condition(mass,1.,line)
    tspan = (0,tfin)
    prob = ODEProblem(prova,U,tspan,0)
    sol = solve(prob,method,
    save_everystep=false,save_start=false,save_end=true,abstol = abstol, reltol = reltol)
    result=sol[1]
   # print(length(result[1,:]))
    f_root =find_zero_R(result[1,:], Rmax, N)
   # @show f_root
    return f_root#(line[ileft]+line[iright])/2-Δx(line)/2

end


function X_scan_infV(;N=200,d=3,Λ=10,Rmax=2.0,mass=-0.9772, tfin=2,abstol=1e-9,reltol=1e-10,method=FBDF())
    
    line=OriginInterval(N,Rmax)
    prova=SemiDiscrete(line,FlowEquation(d,Λ)) 
 
    U=initial_condition(mass,1.,line)
    tspan = (0,tfin)
    prob = ODEProblem(prova,U,tspan,0)
    sol = solve(prob,method,save_everystep=false,save_start=false,save_end=true,abstol = abstol, reltol = reltol)
    result=sol[1]
   # print(length(result[1,:]))
    #f_root =find_zero_X(result[1,:],result[2,:] , Rmax, N)
    f_root = zero_sign(result[1,:], Rmax, N)
    plot( result[2,:]  )
   # @show f_root
    return result[2,f_root]

end



function X_plot(;N=200,d=3,Λ=10,Rmax=2.0,mass=-0.9772, tfin=2,abstol=1e-9,reltol=1e-10,method=FBDF())
    
    line=OriginInterval(N,Rmax)
    prova=SemiDiscrete(line,FlowEquation(d,Λ)) 
 
    U=initial_condition(mass,1.,line)
    tspan = (0,tfin)
    prob = ODEProblem(prova,U,tspan,0)
    sol = solve(prob,method,save_everystep=false,save_start=false,save_end=true,abstol = abstol, reltol = reltol)
    result=sol[1]
    #f_root =find_zero_X(result[1,:],result[2,:] , Rmax, N)
    f_root = zero_sign(result[1,:], Rmax, N)
    p1 =plot( 1 ./result[2,:]  )
   # @show f_root
    return p1

end


function Uprime_scan_infV(;N=200,d=3,Λ=10,Rmax=2.0,mass=-0.9772, tfin=2,abstol=1e-9,reltol=1e-10,method=FBDF())
    
    line=OriginInterval(N,Rmax)
    prova=SemiDiscrete(line,FlowEquation(d,Λ)) 
 
    U=initial_condition(mass,1.0,line)
    tspan = (0,tfin)
    prob = ODEProblem(prova,U,tspan,0)
    sol = solve(prob,method,save_everystep=false,save_start=false,save_end=true,abstol = abstol, reltol = reltol)
    result=sol[1]
  #  @show (result[1,:])
    f_root =find_uprime_0(result[1,:] , Rmax, N)
    
    #f_root = find_uprime_mass(result[1,:], Rmax, N)
    #@show f_root
    return f_root

end

# this function interpolates with Spline 
function find_zero_R(x, Rmax, N)
    line=OriginInterval(N,Rmax).x
    spl = Spline1D(line, x)
    a=2.
    try
        a = fzero(spl, 0,Rmax)
    catch error
        a = 0 
   #     @show a 
    end
   # @show a
    return a
end


# this function checks when the argument x changes sign and returns the index corresponding to the right point 
function zero_sign(x,Rmax, N)
    line=OriginInterval(N,Rmax).x
    for i in 1:length(x)-1
        if x[i]*x[i+1] < 0 
            return i+1
        end 
    end  
    return 1
end 

##############################################################################
function find_uprime_0(x, Rmax, N)
    line=OriginInterval(N,Rmax).x
    dx= line[2]-line[1]
    spl = Spline1D(line, x; k=3, bc="nearest")
    a=2.
    try
        a = fzero(spl, 0,Rmax)
    catch error
        a = 0 
   #     @show a 
    end
    @show a
    return derivative(spl, a; nu=1) # nu is the order of the derivative 
end
##################################################################
function find_uprime_mass(x, Rmax, N)
    
    line=OriginInterval(N,Rmax).x
    dx= line[2]-line[1]
     for i in 2:length(x)
        if x[i]>0 && x[i-1]<=0
          #  return (-x[i+6]+x[i+7])/dx
          @show( (-x[i+6]+x[i+7])/dx)
            return (-x[i+10]+x[i+11])/dx
        end  
    end
    @show( (-x[1]+x[2])/dx)
    return (x[2]-x[1])/dx
 end

# This function is used to find the value of X in the zero of Uprime
function find_zero_X(x,y, Rmax, N)
    line=OriginInterval(N,Rmax).x
    spl = Spline1D(line, x; k=1, bc="nearest")
    sply = Spline1D(line, exp.(y); k=1, bc="nearest")
    
    a=2.
    try
        a = fzero(spl, 0,Rmax)
    catch error
        a = 0 
   #     @show a 
    end
   # @show a
    return sply(a)
end


function mass_to_final_solution(mass,std,traj;N=800,d=3,Λ=10,Rmax=2.0, tfin=50,abstol=1e-9,reltol=1e-9,method=FBDF(),L=16)
    line=OriginInterval(N,Rmax)
    prova=SemiDiscrete(OriginInterval(N,Rmax),FlowEquation(d,Λ,L))
    U=initial_condition(mass,1.0,line)
    tspan = (0,tfin)
    prob = ODEProblem(prova,U,tspan,0)
    Random.seed!(123)
    function prob_func(prob,i,repeat)
        newu0=initial_condition(rand(Normal(mass,std )),1.0,line)
        @show i
        remake(prob,u0=newu0)   
    end

    ensemble_prob = EnsembleProblem(prob,prob_func=prob_func
    ,output_func = (sol,i) -> (sol[:,:,end],false),)
    #sim = solve(ensemble_prob,FBDF(),EnsembleThreads(),trajectories=2) # , abstol = 1e-5, reltol = 1e-5,
    sim = solve(ensemble_prob,FBDF(),trajectories=traj,abstol = abstol, reltol = reltol)
end


function to_final_soolution(;N=800,d=3,Λ=10,Rmax=2.0,mass=0.9278124, tfin=50,abstol=1e-9,reltol=1e-9,method=FBDF(),L=16)
    @show L
    line=OriginInterval(N,Rmax)
    prova=SemiDiscrete(OriginInterval(N,Rmax),FlowEquation(d,Λ,L))
    U=initial_condition(mass,1.0,line)
    tspan = (0,tfin)
    prob = ODEProblem(prova,U,tspan,0)
    sol = solve(prob,method,
    save_everystep=false,save_start=false,save_end=true,abstol = abstol, reltol = reltol)
    return sol[1] 
end


function tominimum(;N=200,d=3,Λ=10,Rmax=2.0,mass=-0.9772, tfin=2,abstol=1e-3,reltol=1e-3,method=FBDF())
    
    line=OriginInterval(N,Rmax)
    prova=SemiDiscrete(line,FlowEquation(d,Λ)) 
 
    U=initial_condition(mass,1.0,line)
    tspan = (0,tfin)
    prob = ODEProblem(prova,U,tspan,0)
    sol = solve(prob,method,
    save_everystep=false,save_start=false,save_end=true,abstol = abstol, reltol = reltol)
    result=sol[1]
    ((ileft,iright), (Vf,Vi))=find_zero_R(result[1,:],0)
    if  ileft == iright 
        return measurement(0,Δx(line))
    else 
        return measurement((line[ileft]+line[iright])/2,Δx(line))
    end 
end
function inerrorbar(x,y)
      z = abs(x - y)
      if ( Measurements.value(z)<  Measurements.uncertainty(z)   )
          return true 
      end
      return false
  end 

# a : one boundary
# b : the other boundary
# assumes that one side is identical to zero (thats why special bisection)
function special_bisection(f::Function, a::T, b::T; tol=1.0e-5, maxIter=100) where {T<:Number}
    fa = f(a)
    fb = f(b)
    iter = zero(maxIter)
    c = (a+b)/2

    @show (a,b)

    # while abs(b-a) > tol && iter < maxIter
    while max(abs(fa), abs(fb)) > tol && iter < maxIter
        iter += 1
        c = (a+b)/2
        fc = f(c)

        if (iszero(fa) && iszero(fc)) || (iszero(fb) && !iszero(fc))
            a = c
            fa = fc
        else
            b = c
            fb = fc
        end

        @show (a,b)
    end

    if iter == maxIter
        println("Maximum iterations reached in `special_bisection`, be careful!")
    end

    return c
end

# a : one boundary
# b : the other boundary
# assumes different signs of f(a) and f(b)
function bisection(f::Function, a::T, b::T; tol=1.0e-5, maxIter=100) where {T<:Number}
    fa = f(a)
    fb = f(b)
    iter = zero(maxIter)
    c = (a+b)/2

    @show (a,b)

    # while abs(b-a) > tol && iter < maxIter
    while max(abs(fa), abs(fb)) > tol && iter < maxIter
        iter += 1
        c = (a+b)/2
        fc = f(c)

        if iszero(fc)
            break
        end

        if fa * fc > 0
            a = c
            fa = fc
        else
            b = c
            fb = fc
        end

        @show (a,b)
    end

    if iter == maxIter
        println("Maximum iterations reached in `special_bisection`, be careful!")
    end

    return c
end


# mass_to_minimum(0.91, 0.94; N=150,d=3,Λ=10,Rmax=2.0, tfin=10.0,abstol=1e-9,reltol=1e-9,method=FBDF())
# mass_to_minimum2(0.91, 0.94; N=150,d=3,Λ=10,Rmax=2.0, tfin=10.0,abstol=1e-9,reltol=1e-9,method=FBDF())

function mass_to_minimum2(mass_a, mass_b; N=200,d=3,Λ=10,Rmax=2.0, tfin=2,abstol=1e-3,reltol=1e-3,method=FBDF())
    Δ = Δx(OriginInterval(N,Rmax)) # copy pasted from above, look here for bug if `tominimum` was changed
    f = x->tominimum(;N=N,d=d,Λ=Λ,Rmax=Rmax, mass = x, tfin=tfin,abstol=abstol,reltol=reltol,method=method).val
    tol = 10e-1 * Δ

    minMass = special_bisection(f, mass_a, mass_b, tol=tol)
    println("Estimated critical bare mass, now error estimation.")
    upperMass = bisection(x->f(x)-Δ, minMass, max(mass_a, mass_b), tol=tol)

    measurement(minMass, upperMass - minMass)
end



function mass_to_minimum(massin=0.9772, massfin=0.9775;N=200,d=3,Λ=10,Rmax=2.0, tfin=2,abstol=1e-3,reltol=1e-3,method=FBDF())
    if massin < massfin 
        mass_L = massin
        mass_R = massfin
    elseif massin > massfin 
        mass_L = massfin
        mass_R = massin
    end 
    aver = (mass_L + mass_R)/2
    phi_L = tominimum(;N=N,d=d,Λ=Λ,Rmax=Rmax, mass = mass_L, tfin=tfin,abstol=abstol,reltol=reltol,method=method)


    phi_R = tominimum(;N=N,d=d,Λ=Λ,Rmax=Rmax, mass = mass_R, tfin=tfin,abstol=abstol,reltol=reltol,method=method)
    phi_av = tominimum(;N=N,d=d,Λ=Λ,Rmax=Rmax, mass = aver, tfin=tfin,abstol=abstol,reltol=reltol,method=method)
    println("..............")
    @show phi_L
    @show phi_R
    @show mass_R
    @show mass_L
  #  @show mass_R
  #  @show tfin
    L = Measurements.value(phi_L)
    R = Measurements.value(phi_R)
    A = Measurements.value(phi_av)
    
    if inerrorbar(phi_L, phi_R) 
       # print( phi_L, phi_R, phi_av)
        return (mass_L, mass_R)
    elseif (L == 0. && R >0.)

        if A> 0.
            mass_to_minimum(mass_L, aver;N=N,d=d,Λ=Λ,Rmax=Rmax, tfin=tfin,abstol=abstol,reltol=reltol,method=method)
        else  
            mass_to_minimum(aver, mass_R;N=N,d=d,Λ=Λ,Rmax=Rmax, tfin=tfin,abstol=abstol,reltol=reltol,method=method)
        end 
    elseif  (L >0.)
        throw(DomainError("cambia m_L") ) 
    elseif  (R ==0.)
        throw(DomainError("cambia m_R") )   
    end 
end 

function find_zero_phi(x)
    for i in 2:length(x)
        if x[i]>0 && x[i-1]<=0
            return ((i-1,i),(x[i-1],x[i])  )
        end  
    end
    return ((1,1),(0,0)  )
end

function find_zero_phi(x,discretization)
    for i in 2:length(x)
        if x[i]>0 && x[i-1]<=0
            return ((discretization[i-1],discretization[i]),(x[i-1],x[i])  )
        end  
    end
    return ((1,1),(0,0)  )
end




##############
function interpolatezero(sol,prova,t, tre)
((i_minus,i_plus), (V_minus, V_plus ) )=find_zero_R(sol(t)[1,:], tre)
    if i_minus==i_plus
        return 0
    end
    phi_minus=prova.spatial_discretization.x[i_minus]
    phi_plus=prova.spatial_discretization.x[i_plus]

    phi  = V_minus /(V_plus-V_minus)*(phi_plus-phi_minus) +phi_minus

    return  phi
end

#####################################################

function mass_scan_infV(;min,mfin,nmass, tfin)

    massrange=collect(range(min,mfin,nmass))
   # print(massrange)
    phimin=zeros(nmass)
    for i in 1:nmass 
        phimin[i]=temperaturescan_infV(N=400, Rmax = 2,  mass=massrange[i], tfin=tfin, tre=0.0)
    end
    return (phimin, massrange)
end

function mass_scan(;min,mfin,nmass, tfin, l)

    massrange=collect(range(min,mfin,nmass))
   # print(massrange)
    phimin=zeros(nmass)
    for i in 1:nmass 
        phimin[i]=temperaturescan(N=400, Rmax = 2,  mass=massrange[i], tfin=tfin, L=l, tre=0.0)
    end
    return (phimin, massrange)
end

function volumescan(N=400, Rmax=2, mass=.07, tfin=5) 
    ls = [16,20,24,28,32,36, 48, 60]
    
    a = zeros(length(ls))
    for i in 1:(length(ls))
     
        a[i] = temperaturescan(;N,d=3,Λ=1,Rmax,mass, tfin, L=ls[i], tre=0.001)
    end
    return ls, a
end

function G1(sol, t_fin, a)
    x,flowequation=prova.spatial_discretization,prova.equation_of_motion
    dx=Δx(x)
    x_c = find_zero_R(sol(t_fin)[1,:] , 0.0)
    print(x_c)
    U2_c = (sol(t_fin)[1,x_c[1][1] + a+1] - sol(t_fin)[1,x_c[1][1]+a])/dx

    X_crit = sol(t_fin)[2,x_c[1][1]+a]
    return X_crit, U2_c
end


function G_1cent(sol, t_fin)
    x,flowequation=prova.spatial_discretization,prova.equation_of_motion
    dx=Δx(x)
    U2_c = (1/12*sol(t_fin)[1,3] +2/3*sol(t_fin)[1,2] + 2/3*sol(t_fin)[1,1] +1/12*sol(t_fin)[1,2] )/dx
    X_crit = sol(t_fin)[2,1]

    return X_crit, U2_c
end

function G_1cent(sol,N,Rmax)

    x=OriginInterval(N,Rmax)
    dx=Δx(x)
    U2_c = (1/12*sol[1,3] +2/3*sol[1,2] + 2/3*sol[1,1] +1/12*sol[1,2] )/dx
    X_crit = sol[2,1]

    return X_crit, U2_c
end


function L_G1(N, t_fin, Λ, mass_crit, ls)
    
    a = zeros(length(ls))
    b = zeros(length(ls))
    for i in 1:(length(ls))
     
        
        prova=SemiDiscrete(OriginInterval(N,Rmax),FlowEquation(d,Λ,ls[i]))
        U=zeros(2,length(prova.spatial_discretization.x))
        for i in axes(U,2)
        U[1,i] = -mass_crit*prova.spatial_discretization.x[i]+prova.spatial_discretization.x[i]^3 
        U[2,i]= 0
        end
        tspan = (0,t_fin)
        prob = ODEProblem(prova,U,tspan,0)
        print(prova.spatial_discretization.x)
        sol = solve(prob,FBDF(), abstol = 1e-5, reltol = 1e-5)
        println(ls[i] )
        a[i] = G_1cent(sol, t_fin)[1]
        b[i] = G_1cent(sol, t_fin)[2]
    end
    return ls, a, b
end
function timescan(N, L, time_array, mass)
   
    X_t = zeros(length(time_array))
    U_t = zeros(length(time_array))
    for i in 1:length(time_array)
         
            
            prova=SemiDiscrete(OriginInterval(N,Rmax),FlowEquation(d,Λ))
            U=zeros(2,length(prova.spatial_discretization.x))
            for i in axes(U,2)
            U[1,i] = -mass*prova.spatial_discretization.x[i]+prova.spatial_discretization.x[i]^3 
            U[2,i]= 1
            end
            tfin = time_array[i]
            tspan = (0,tfin)
            prob = ODEProblem(prova,U,tspan,0)
            print(prova.spatial_discretization.x)
            sol = solve(prob,FBDF();progress = true,
            progress_steps = 100)
            X_t[i] = G1(sol, t_fin, 1)[1]
            U_t[i] = G1(sol, t_fin, 1)[2]
    end
    return X_t, U_t
    
end

function G1_t(t_G, a)
    y = zeros(length(t_G))
    for i in 1:(length(y))
        y[i]=exp(-t_G[i]*a)
    end
    return y
end

function plot_G1(t, a, L)
    p =scatter(t,G1_t(t, a[1]), label=L[1], xlabel="t", ylabel="G_1(t)",title = "Correlation" )
    for i in 2:length(a)
        scatter!(t,G1_t(t, a[i]) , label=L[i])
    end
p
end

# Plot 2PF rescaled by the critical exponent z
function plot_G1_rescaled(t, a, L, z)
    p =scatter(t/L[1]^z,(G1_t(t, a[1])), label=L[1], xlabel="t/L^z", ylabel="G_1(t)",title =z )
    for i in 2:length(a)
        scatter!(t/L[i]^z, (G1_t(t, a[i]))  , label=L[i])
    end
p
end



function TPFplot(mass, exponent)
    ls = [ 16,24, 32, 48,56,64]
    ls, xx, uu= L_G1(800, 50, 10, mass, ls)
    tau = uu./xx
   # scatter(log.(ls), -log.(tau))

    t_G =range(0,stop=1000,step=10)
    plot_G1_rescaled(t_G, abs.(tau), ls,exponent)
   # scatter!(t_G, g1inf)
   return 0
end


function Binder_Cumulant(sol, t_fin )
    #find the zero of U_prime 
    x_c = find_zero_R(sol(t_fin)[1,:] , 0.0)
    print(x_c)#\
    x = x_c[1][1]+2
    t_fin=10
    U_2prime_i = (sol(t_fin)[1, x + 1] - sol(t_fin)[1, x])/dx
    U_2prime_i1 = (sol(t_fin)[1, x + 2] - sol(t_fin)[1, x+1])/dx
    U_2prime_i2 = (sol(t_fin)[1, x + 3] - sol(t_fin)[1, x+2])/dx
    U_2prime_i3 = (sol(t_fin)[1, x + 4] - sol(t_fin)[1, x+3])/dx
    U_4prime = (-U_2prime_i3 + 4*U_2prime_i2 - 5*U_2prime_i1 + 2*U_2prime_i)/dx^2
    return U_2prime_i,U_4prime,  1 - U_4prime*U_2prime_i^2/(3)
end
