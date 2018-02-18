# Weighted least squares in 3D space.
# Approximation will be done with polynomials.
# Assume that approximated function will be in form f: R^2 -> R
# The primary purpose of this program is to find all mistakes, that will
# happen in the real program and solve them early.

include("point.jl")
using Point3D
using ArgParse
using JSON
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
  for i in 1:pointCount
    push!(toReturn, Point(parsedPoints[i][1], parsedPoints[i][2], parsedPoints[i][3]))
  end
  return toReturn
end

function parseQuery(query::String)
  return map(x -> parse(Float64, x), split(chomp(query)))
end

function wls(samplePoints::Vector{Point}, inputPoints)

end

#TODO: dysfunctional cli prompt
function main()
  options = parseCommandline()
  samplePoints = parseInputPoints(options["input"])
  while true
    print("Input x and y > ")
    input = readline(Base.STDIN)
    if eof(Base.STDIN)
      break
    end
    println(input)
    parsedInput = parseQuery(input)
    wls(samplePoints, parsedInput)
  end
  println("Bye.")
end

main()

