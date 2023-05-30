
abstract type SpatialDiscretization end



struct Interval{N,R1,R2,DX,T} <: SpatialDiscretization
    x::Vector{T}
end

function Interval(n,r1,r2) 
    data=[r1+(r2-r1)/(2*n)+(i-1)*(r2-r1)/(n) for i in 1:n]
    T=eltype(data)
    Interval{n,r1,r2,(r2-r1)/n,T}(data)
end

function SymmetricInterval(n,r2) 
    r1=-r2
    data=[r1+(r2-r1)/(2*n)+(i-1)*(r2-r1)/(n) for i in 1:n]
    T=eltype(data)
    Interval{n,r1,r2,(r2-r1)/n,T}(data)
end

function OriginInterval(n,r2) 
    r1=0
    data=[r1+(r2-r1)/(2*n)+(i-1)*(r2-r1)/(n) for i in 1:n]
    T=eltype(data)
    Interval{n,r1,r2,(r2-r1)/n,T}(data)
end

function PeriodicInterval(n) 
    r1=0
    r2=2*pi
    data=[r1+(i-1)*(r2-r1)/(n) for i in 1:n]
    T=eltype(data)
    Interval{n,r1,r2,(r2-r1)/n,T}(data)
end

function InTheOriginInterval(n,r2)
    r1=0
    data=[r1+(i-1)*(r2-r1)/(n) for i in 1:n]
    T=eltype(data)
    Interval{n,r1,r2,(r2-r1)/n,T}(data)
end


@inline Base.eltype(::Type{Interval{N,R1,R2,DX,T}}) where {N,R1,R2,DX,T} =T
@inline Base.length(::Interval{N,R1,R2,DX,T}) where {N,R1,R2,DX,T} =N
@inline Base.size(::Interval{N,R1,R2,DX,T})  where {N,R1,R2,DX,T} =(N,)

@inline Base.IndexStyle(::Type{<:Interval}) = IndexLinear()
@inline Base.getindex(S::Interval, i::Int) = S.x[i]
@inline Base.firstindex(S::Interval) = 1
@inline secondindex(S::Interval) = 2
@inline Base.lastindex(S::Interval{N,R1,R2,DX,T}) where {N,R1,R2,DX,T}=N 
@inline secondtolastindex(S::Interval{N,R1,R2,DX,T}) where {N,R1,R2,DX,T}=N-1


@inline interior(S::Interval{N,R1,R2,DX,T}) where {N,R1,R2,DX,T} =  2:N-1

@inline interior2(S::Interval{N,R1,R2,DX,T}) where {N,R1,R2,DX,T} =  3:N-2

@inline interior2_left(S::Interval{N,R1,R2,DX,T}) where {N,R1,R2,DX,T} =  2:N-2

@inline Δx(S::Interval{N,R1,R2,DX,T}) where {N,R1,R2,DX,T} = DX







