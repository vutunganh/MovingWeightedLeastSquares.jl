include("point.jl")

using DynamicPolynomials
using MultivariatePolynomials
using Base.Random
using ArgParse
using JSON
using Point3D

@polyvar x y

"Uniformly generates a random float from a to b (b noninclusive)"
function randUniform(a, b)
  return a + rand() * (b - a)
end

"""
Generates a random polynomial
TODO: how to make this parametrizable by max degree, 
e.g. generatePolynomial(3) -> x^1, x^2, x^3, x^2*y, x*y^2,...?
btw `@polyvar a[1:n]` works
"""
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

  return parse_args(s)
end

function main()  
  options = parseCommandline()

  randomPolynomial = generatePolynomial()
  randomPoints = generatePoints(options["amount"])
  for p in randomPoints
    p.z = calcOutput(p, randomPolynomial)
  end

  toOutput = Dict("polynomial" => randomPolynomial,
                  "points" => randomPoints)
  jsonOutput = json(toOutput)
  if options["output"] == nothing
    outputStream = STDOUT
  else
    outputStream = open(options["output"], "w")
  end
  write(outputStream, jsonOutput)
  close(outputStream)
end

main()
