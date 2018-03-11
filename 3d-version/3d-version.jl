# Weighted least squares in 3D space.
# Approximation will be done with polynomials.
# Assume that approximated function will be in form f: R^2 -> R
# The primary purpose of this program is to find all mistakes, that will
# happen in the real program and solve them early.

include("point.jl")
include("polynomial-generator.jl")

__precompile__()
module MWLS3D
using Point3D
using PolynomialGenerator
using ArgParse
using DataFrames
using CSV
using DynamicPolynomials
using MultivariatePolynomials

export wls, parseInputPoints

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

function weightingFunction(d, EPS)
  return 1 / (d^2 + EPS^2)
end

function dd(p::Point, q::Point)
  sqrt((p.x - q.x)^2 + (p.y - q.y)^2)
end

function wls(samplePoints::Vector{Point}, EPS::Float64, dimensionCount::Int64, variableCount::Int64, inputPoint::Point)
  println(inputPoint)
  vars::Vector{PolyVar{true}}, b::Array{Monomial{true}} = gen(variableCount, dimensionCount)
  basis = b * transpose(b)
  println(basis)
  firstTerm = zeros(length(b), length(b))
  secondTerm = zeros(length(b))

  for s in samplePoints
    distance = weightingFunction(dd(inputPoint, s), EPS)
    for i in 1:length(b)
      for j in 1:length(b)
        firstTerm[i, j] += (distance * subs(basis[i, j], vars => (s.x, s.y)))
      end
      secondTerm[i] += (distance * subs(b[i], vars => (s.x, s.y)) * s.z)
    end
  end
  result = firstTerm \ secondTerm
  println(inv(firstTerm))
  println(secondTerm)
  return result
end
end

