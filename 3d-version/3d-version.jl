# Weighted least squares in 3D space.
# Approximation will be done with polynomials.
# Assume that approximated function will be in form f: R^2 -> R
# The primary purpose of this program is to find all mistakes, that will
# happen in the real program and solve them early.

include("point.jl")
include("polynomial-generator.jl")
using Point3D
using PolynomialGenerator
using ArgParse
using DataFrames
using CSV
using DynamicPolynomials
using MultivariatePolynomials

function parseCommandline()
  s = ArgParseSettings()
  @add_arg_table s begin
    "input"
      help = "a csv file with sample points"
      required = true
      arg_type = String
  end

  return parse_args(s)
end

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

function parseQuery(query::String)
  return map(x -> parse(Float64, x), split(chomp(query)))
end

function weightingFunction(d, EPS)
  return 1 / (d^2 + EPS^2)
end

function dd(p::Point, q::Point)
  return sqrt((p.x - q.x)^2 + (p.y - q.y)^2)
end

function wls(samplePoints::Vector{Point}, inputPoint::Point)
  vars::Vector{PolyVar{true}}, b::Vector{Monomial{true}} = gen(2, 2) # TODO: parametrize this
  m = length(b)
  n = length(samplePoints)
  A = [weightingFunction(dd(inputPoint, samplePoints[j]), 0.0001) * Float64(subs(b[i], vars => (samplePoints[j].x, samplePoints[j].y))) for j in 1:n, i in 1:m] # just beautiful
  firstTerm = inv(transpose(A) * A)
  println(firstTerm)
end

# TODO: dysfunctional cli prompt
function main()
  options = parseCommandline()
  samplePoints = parseInputPoints(options["input"])
  while true
    print("Input x and y > ")
    input = readline()
    if input == ""
      break
    end
    parsedInput = parseQuery(input)
    push!(parsedInput, 0)
    inputPoint = Point(parsedInput[1], parsedInput[2], parsedInput[3])
    wls(samplePoints, inputPoint)
  end
  println("")
  println("Bye.")
end

main()

