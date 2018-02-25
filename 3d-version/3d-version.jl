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
using JSON
using DynamicPolynomials
using MultivariatePolynomials

function parseCommandline()
  s = ArgParseSettings()
  @add_arg_table s begin
    "input"
      help = "input json"
      required = true
      arg_type = String
  end

  return parse_args(s)
end

"Reads the input and returns an array of points"
function parseInputPoints(inputFilename::String)
  inputStream = open(inputFilename, "r")
  rawInput = readstring(inputStream)
  close(inputStream)
  parsedInput = JSON.parse(rawInput)
  pointCount::Int64 = length(parsedInput["points"])
  parsedPoints = parsedInput["points"]
  toReturn::Vector{Point} = []
  for p in parsedPoints
    push!(toReturn, Point(p[1], p[2], p[3]))
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
  for p in samplePoints
    theta::Real = dd(inputPoint, p)
    bx = subs(sum(b), vars[1] => p.x, vars[2] => p.y)
    weightingFunction(theta, 0.0001) * bx * transpose(bx)
  end
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
  println("Bye.")
end

main()

