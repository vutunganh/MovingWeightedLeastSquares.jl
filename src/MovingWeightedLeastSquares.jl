__precompile__()
module MovingWeightedLeastSquares

const Point = Vector{Float64}

using DataFrames
using CSV
using DynamicPolynomials
using MultivariatePolynomials

export parseInputPoints
include("io.jl")
include("polynomial-generator.jl")
export generateRandomData
include("sample-generator.jl")
include("3d-version.jl")

export 
end
