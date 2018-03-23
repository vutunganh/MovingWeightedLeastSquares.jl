__precompile__()
module MovingWeightedLeastSquares

const Point = Vector{Float64}

using DynamicPolynomials
using MultivariatePolynomials
using Plots

export readData
include("io.jl")
include("polynomial-generator.jl")
export generateRandomData
include("sample-generator.jl")
export WlsObject, wls, plotWls
include("wls-object.jl")
include("wls-ops.jl")

end

