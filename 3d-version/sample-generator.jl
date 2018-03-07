include("point.jl")

using DynamicPolynomials
using MultivariatePolynomials
using Base.Random
using ArgParse
using DataFrames
using CSV
using Point3D

@polyvar x y

"Uniformly generates a random float from a to b (b noninclusive)"
function randUniform(a, b)
  return a + rand() * (b - a)
end

"Generates a random polynomial"
function generatePolynomial()
  coefficients = [randUniform(-100.0, 100.0) for i in 1:6]
  return MultivariatePolynomials.polynomial([1, x, y, x^2, x*y, y^2], coefficients)
end

"Generates x and y coordinates randomly"
function generateSinglePoint()
  return Point(randUniform(-100.0, 100.0), randUniform(-100.0, 100.0), 0)
end

"Generates `n` random points without outputs"
function generatePoints(n::Int64)
  return [generateSinglePoint() for i in 1:n]
end

"Evaluates the polynomial `pol` at `p` "
function calcOutput(p, pol)
  return pol(x => p.x, y => p.y)
end

function parseCommandline()
  s = ArgParseSettings()
  @add_arg_table s begin
    "-n", "--amount"
      help = "amount of generated points"
      arg_type = Int
    "-o", "--output"
      help = "output file for generated points"
      arg_type = String
  end
  toReturn = parse_args(s)
  if toReturn["amount"] == nothing
    println("Specify the amount of generated points")
    exit(-1)
  end
  if toReturn["output"] == nothing
    println("Specify the output file")
    exit(-1)
  end

  return parse_args(s)
end

function main()  
  options = parseCommandline()

  randomPolynomial = generatePolynomial()
  randomPoints = generatePoints(options["amount"])
  for p in randomPoints
    p.z = calcOutput(p, randomPolynomial)
  end

  df = DataFrame()
  for f in 1:length(fieldnames(randomPoints[1]))
    col = [x[f] for x in randomPoints]
    df[Symbol("x", f)] = col
  end
  CSV.write(options["output"], df)
  println("Done")
end

main()
