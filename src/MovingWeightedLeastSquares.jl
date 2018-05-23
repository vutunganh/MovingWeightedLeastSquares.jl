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
export CellLinkedList, cllAdd!, cll_add!, cllRemove!, cll_remove!, cllModify!, cll_modify!, cllInrange, cll_inrange
include("cell-linked-list.jl")
export MwlsObject, mwlsNaive, mwls_naive, mwlsKd, mwls_kd, mwlsCll, mwls_cll
include("mwls-object.jl")
export calcMwlsCoefficients, calc_mwls_coefficients,
       approximate, calcDiffMwlsPolys,
       mwls_diff_polys, mwlsDiff, mwls_diff
include("mwls-ops.jl")

end
