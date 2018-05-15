__precompile__()
module MovingWeightedLeastSquares

const Point = Vector{T} where {T <: Real}

using DynamicPolynomials
using MultivariatePolynomials
using NearestNeighbors
using DataStructures
using Plots

include("polynomial-generator.jl")
export WlsObject, wls, calcWlsCoefficients, plotWls
include("wls-object.jl")
include("wls-ops.jl")
export CellLinkedList, cllAdd!, cllRemove!, cllModify!, cllIteratedCells, cllInrange
include("cell-linked-list.jl")
export MwlsObject, mwlsNaive, mwlsKd, mwlsCll
include("mwls-object.jl")
export calcMwlsCoefficients, calcDiffMwlsPolys, mwlsDiff
include("mwls-ops.jl")

end

