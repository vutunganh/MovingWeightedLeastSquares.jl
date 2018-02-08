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
function parseInput(inputFilename::String)
  inputStream = open(inputFilename, "r")
  rawInput = readstring(inputStream)
  close(inputStream)
  parsedInput = JSON.parse(rawInput)
  pointCount::Int64 = length(parsedInput["points"])
  parsedPoints = parsedInput["points"]
  toReturn::Vector{Point} = []
  for i in 1:1
    push!(toReturn, Point(parsedPoints[i]["x"], parsedPoints[i]["y"], parsedPoints[i]["z"]))
  end
  return toReturn
end

#TODO: dysfunctional cli prompt
function main()
  options = parseCommandline()
  parseInput(options["input"])
  while true
    print("Input x and y > ")
    input = readline(STDIN)
    println(chomp(input))
    println("")
    flush(STDIN)
    if eof(STDIN)
      println("")
      break
    end
  end
  println("Bye.")
end

main()

