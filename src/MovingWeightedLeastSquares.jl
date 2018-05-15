__precompile__()
module MovingWeightedLeastSquares

const Point = Vector{T} where {T <: Real}

using DynamicPolynomials
using MultivariatePolynomials
using NearestNeighbors
using DataStructures

include("polynomial-generator.jl")
export WlsObject, wls, calcWlsCoefficients, plotWls
include("wls-object.jl")
include("wls-ops.jl")
export CellLinkedList, cllAdd!, cllRemove!, cllModify!, cllIteratedCells, cllInrange
include("cell-linked-list.jl")
export MwlsObject,
       mwlsNaive, mwlsKd, mwlsCll,
       calcMwlsCoefficients, calcDiffMwlsPolys, mwlsDiff
include("mwls-object.jl")
include("mwls-ops.jl")

end

