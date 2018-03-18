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
export WlsObject, wls
include("wls-object.jl")
include("wls-ops.jl")

end

