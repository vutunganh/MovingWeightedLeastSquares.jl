__precompile__()
module MovingWeightedLeastSquares

const Point = Vector{Float64}

using DynamicPolynomials
using MultivariatePolynomials
using NearestNeighbors
using Plots

export parseInputPoints
include("io.jl")
include("polynomial-generator.jl")
export generateRandomData
include("sample-generator.jl")
export WlsObject, wls, plotWls
include("wls-object.jl")
include("wls-ops.jl")
export MwlsObject, mwls, calcMwlsCoefficients
include("mwls-object.jl")
include("mwls-ops.jl")

end

