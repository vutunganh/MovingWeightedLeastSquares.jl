# Weighted least squares in 3D space.
# Approximation will be done with polynomials.
# Assume that approximated function will be in form f: R^2 -> R
# The primary purpose of this program is to find all mistakes, that will
# happen in the real program and solve them early.

include("point.jl")
include("polynomial-generator.jl")

Base.__precompile__()
module MWLS3D
using Point3D
using PolynomialGenerator
using ArgParse
using DataFrames
using CSV
using DynamicPolynomials
using MultivariatePolynomials

export wls, parseInputPoints, wlsObject

"Reads the input and returns an array of points"
function parseInputPoints(inputFilename::String)
  parsedInput = CSV.read(inputFilename)
  pointCount::Int64 = size(parsedInput, 1)
  dimensionCount::Int64 = size(parsedInput, 2)
  toReturn::Vector{Point} = []
  for i in 1:pointCount
    p = parsedInput[i, :]
    apt = [p[1, j] for j in 1:dimensionCount]
    push!(toReturn, Point(apt))
  end
  return toReturn
end

function dd(p::Point, q::Point)
  sqrt((p.x - q.x)^2 + (p.y - q.y)^2)
end

struct WlsObject
  vars::Vector{PolyVar{true}}
  b::Array{Monomial{true}}
  data::Vector{Point}
  EPS::Float64
  wfun
  basis
end

function WlsObject(vars, b, data, EPS, wfun)
  WlsObject(vars, b, data, EPS, wfun, b * transpose(b))
end

"""
this function approximates the WlsObject at inputPoint
"""
function (w::WlsObject)(inputPoint::Point)
  m = length(w.b)
  firstTerm = zeros(m, m)
  secondTerm = zeros(m)

  for p in w.data
    dist = w.wfun(dd(inputPoint, p), w.EPS)
    for i in 1:m
      for j in 1:m
        firstTerm[i, j] += dist * subs(w.basis[i, j], w.vars => (p.x, p.y))
      end
      secondTerm[i] += dist * subs(w.b[i], w.vars => (p.x, p.y)) * p.z
    end
  end

  result = firstTerm \ secondTerm
end

"""
this function creates a wls object

wfun is a weighting function, that should be used
it should take one parameter
"""
function wls(samplePoints::Vector{Point}, variableCount::Int64, dimensionCount::Int64, EPS::Float64, wfun)
  vars::Vector{PolyVar{true}}, b::Array{Monomial{true}} = gen(variableCount, dimensionCount)
  return WlsObject(vars, b, samplePoints, EPS, wfun)
end
end

