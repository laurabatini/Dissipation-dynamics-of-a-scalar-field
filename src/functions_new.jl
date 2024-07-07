# old file 

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

@inline SolidAngle(::FlowEquation{D,S,O,InfiniteVolume}) where {D,S,O}= O


@inline volume(::FlowEquation{D,S,O,FiniteVolume}) where {D,S,O}= O^D

@inline length_(::FlowEquation{D,S,O,FiniteVolume}) where {D,S,O}= O

@inline volume(::FlowEquation{D,S,O,InfiniteVolume} where {D,S,O} )  = InfiniteVolume()


# Create a FlowEquation object with a specific dimension and scale
function FlowEquation(d, Lambda)
    FlowEquation{d, Lambda, Ω(d), InfiniteVolume}()
end

# Create a FlowEquation object with a specific dimension, scale, and length
function FlowEquation(d, Lambda, L)
    if isa(L, Real)
        FlowEquation{d, Lambda, L, FiniteVolume}()
    elseif isa(L, InfiniteVolume)
        FlowEquation{d, Lambda, Ω(d), InfiniteVolume}()
    else
        throw(ArgumentError("The length of the box is infinity or a specific length"))
    end
end

# Compute the potential flow for a given derivative and time
function potential_flow(du, t, ::FlowEquation{D, S, O, L}) where {D, S, O, L}
    k = S * exp(-t)
    -k^2 / (du + k^2)
end

# Compute the propagator for a given derivative and time
function Propagator(du, t, flow::FlowEquation{D, S, O, T}) where {D, S, O, T}
    k = S * exp(-t)
    1 / (du + k^2)
end

# Compute the prefactor for an infinite volume flow equation

function Prefactor(t, flow::FlowEquation{D, S, O, InfiniteVolume}) where {D, S, O}
    k = S * exp(-t)
    O / (4 * pi * k)
end

# Compute the prefactor for a finite volume flow equation
function Prefactor(t, flow::FlowEquation{D, S, L, FiniteVolume}) where {D, S, L}
    k = S * exp(-t)
    ndimcount(static(D), L * k / (2 * π)) / L^D
end

# Compute the alphaflow for a given derivative and time
function alphaflow(propagator, dpropagator, dpropagator2, dalpha, ddalpha, t, ::FlowEquation{D, S, O, L}) where {D, S, O, L}
    k = S * exp(-t)
    k^2 * 0.5 * (3 * dpropagator^2 + 4 * dpropagator2 * dalpha + 2 * (ddalpha + dalpha^2) * propagator^2)
end

function alphaflow(propagator, dpropagator, dpropagator2, dalpha,dalpha_c, ddalpha, t, ::FlowEquation{D, S, O, L}) where {D, S, O, L}
    k = S * exp(-t)
    k^2 * 0.5 * (3 * dpropagator^2 + 4 * dpropagator2 * dalpha + 2 * (ddalpha + dalpha_c^2) * propagator^2)
end

# Define the SemiDiscrete type
struct SemiDiscrete{T, X}
    spatial_discretization::X
    equation_of_motion::T
end


# Define the call operator for the SemiDiscrete type
function (ff::SemiDiscrete)(du, u, p, t)
    x, flowequation = ff.spatial_discretization, ff.equation_of_motion
    dx = Δx(x)
    invdx = 1 / dx
    invdx2 = 1 / (dx^2)
    prefactor = Prefactor(t, flowequation)
    
    # Compute du for the first grid point
    du[1, 1] = prefactor * (potential_flow((u[1, 2] - u[1, 1]) * invdx, t, flowequation) - potential_flow((u[1, 1] + u[1, 1]) * invdx, t, flowequation)) * invdx
    
    # Compute du for the second grid point
    tau = u[2, 1]
    dtau = (u[2, 1] - u[2, 1]) * invdx
    ddtau = (u[2, 2] + u[2, 1] - 2 * u[2, 1]) * invdx2
    
    du[2, 1] = prefactor * alphaflow(-Propagator((u[1, 1] + u[1, 1]) * invdx, t, flowequation), 
                                      (Propagator((u[1, 2] - u[1, 1]) * invdx, t, flowequation) - 
                                       Propagator((u[1, 1] + u[1, 1]) * invdx, t, flowequation)) * invdx, 
                                      (Propagator((u[1, 2] - u[1, 1]) * invdx, t, flowequation)^2 - 
                                       Propagator((u[1, 1] + u[1, 1]) * invdx, t, flowequation)^2) * invdx, 
                                      dtau, ddtau, t, flowequation)
    
    # Compute du for the interior grid points
    for i in 2:length(x)-1
        propagator_i = -Propagator((u[1, i] - u[1, i-1]) * invdx, t, flowequation)
        propagator_i_1 = -Propagator((u[1, i+1] - u[1, i]) * invdx, t, flowequation)
        dpropagator_i = (propagator_i_1 - propagator_i) * invdx
        dpropagator2_i = (propagator_i_1^2 - propagator_i^2) * invdx
        
        du[1, i] = prefactor * (potential_flow((u[1, i+1] - u[1, i]) * invdx, t, flowequation) - 
                                potential_flow((u[1, i] - u[1, i-1]) * invdx, t, flowequation)) * invdx
        
        dtau = (u[2, i] - u[2, i-1]) * invdx
        ddtau = (u[2, i+1] + u[2, i-1] - 2 * u[2, i]) * invdx2
        tau = u[2, i]
        
        du[2, i] = prefactor * alphaflow(propagator_i, dpropagator_i, dpropagator2_i, dtau, ddtau, t, flowequation)

       
    end
    
    # Compute du for the last grid point
    i = length(x)
    
    du[1, i] = prefactor * (potential_flow((u[1, i] - u[1, i-1]) * invdx, t, flowequation) - 
                            potential_flow((u[1, i-1] - u[1, i-2]) * invdx, t, flowequation)) * invdx
    du[2, i] = du[2, i-1]
end


# Define the call operator for the SemiDiscrete type
function not_symmetric(du, u, p, t)
    x, flowequation = p #ff.spatial_discretization, ff.equation_of_motion
    dx = Δx(x)
    invdx = 1 / dx
    invdx2 = 1 / (dx^2)
    prefactor = Prefactor(t, flowequation)
    
    
    # Compute du for the interior grid points
    for i in 2:length(x)-1
        propagator_i = -Propagator((u[1, i] - u[1, i-1]) * invdx, t, flowequation)
        propagator_i_1 = -Propagator((u[1, i+1] - u[1, i]) * invdx, t, flowequation)
        dpropagator_i = (propagator_i_1 - propagator_i) * invdx
        dpropagator2_i = (propagator_i_1^2 - propagator_i^2) * invdx
        
        du[1, i] = prefactor * (potential_flow((u[1, i+1] - u[1, i]) * invdx, t, flowequation) - 
                                potential_flow((u[1, i] - u[1, i-1]) * invdx, t, flowequation)) * invdx
    
    
        if x[i]>0 
            dtau = (u[2, i] - u[2, i-1]) * invdx
        else 
            dtau = (u[2, i+1] - u[2, i]) * invdx
        end
        dtaucent=( u[2, i+1] - u[2, i-1]) * invdx
        ddtau = (u[2, i+1] + u[2, i-1] - 2 * u[2, i]) * invdx2
        tau = u[2, i]
        
        du[2, i] = prefactor * alphaflow(propagator_i, dpropagator_i, dpropagator2_i, dtaucent,dtaucent, ddtau, t, flowequation)

        #@show x[i] dpropagator2_i, dpropagator_i

    end
    
    # Compute du for the last grid point
    i = length(x)
    
    du[1, i] = du[1, i-1]
    du[2, i] = du[2, i-1]

    i=1 
    du[1, i] = du[1, i+1]
    du[2, i] = du[2, i+1]
end


# Define initial conditions for the system
function initial_condition(mass, lambda, spatial_discretization)
    U = zeros(2, length(spatial_discretization.x))
     
    for i in axes(U, 2)
        U[1, i] = -mass * spatial_discretization.x[i] + lambda * spatial_discretization.x[i]^3
        U[2, i] = 0
    end
    
    return U
end


# Find the value of Uprime at x=0
function find_uprime_0(x, Rmax, N)
    line = OriginInterval(N, Rmax).x
    dx = line[2] - line[1]
    spl = Spline1D(line, x; k=3, bc="nearest")
    a = 0.0
    try
        a = fzero(spl, 0, Rmax)
    catch error
        a = 0.0
        # Handle the error or show a warning message if needed
        println("Error: $error")
    end
    return derivative(spl, a; nu=1) # nu is the order of the derivative
end

# Find the value of X at the zero of Uprime
function find_zero_X(x, y, Rmax, N)
    line = OriginInterval(N, Rmax).x
    spl = Spline1D(line, x; k=1, bc="nearest")
    sply = Spline1D(line, exp.(y); k=1, bc="nearest")
    
    a = 2.0
    try
        a = fzero(spl, 0, Rmax)
    catch error
        a = 0 
   #     @show a 
    end
   # @show a
    return sply(a)
end

# Perform a bisection search on a function with special conditions
# when one side is identical to zero
function special_bisection(f::Function, a::T, b::T; tol=1.0e-5, maxIter=100) where {T<:Number}
    fa = f(a)
    fb = f(b)
    iter = zero(maxIter)
    c = (a+b)/2

    @show (a,b)

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

# Perform a bisection search on a function with different signs at the boundaries
function bisection(f::Function, a::T, b::T; tol=1.0e-5, maxIter=100) where {T<:Number}
    fa = f(a)
    fb = f(b)
    iter = zero(maxIter)
    c = (a+b)/2

    @show (a,b)

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
