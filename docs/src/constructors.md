# Constructors

Constructing `MwlsObject` directly is not recommended and helper functions described below should be used instead.

```@docs
mwlsCll(::Array{T, 2}, ::Real, ::Function) where {T <: Real}
mwlsCll(::Array{T, 2}, ::Real, ::Function; ::Int = 1, ::Int = 2) where {T <: Real}
mwlsKd(inputs::Array{T, N}, outputs::Array{U}, EPS::Real, weightFunc::Function) where {T <: Real, U <: Real, N}
mwlsKd(inputs::Array{T, N}, outputs::Array{U}, EPS::Real, weightFunc::Function; leafsize::Int = 10, maxDegree::Int = 2) where {T <: Real, U <: Real, N}
```
