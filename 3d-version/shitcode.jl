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
  coefficients = [randUniform(-100.0, 100.0) for i in 1:2]
  @polyvar x y
  return MultivariatePolynomials.polynomial([x, y], coefficients)
end

"Generates x and y coordinates randomly"
function generateSinglePoint()
  return Point(randUniform(-100.0, 100.0), randUniform(-100.0, 100.0), 0)
end

"Generates `n` random points without outputs"
function generatePoints(n::Int64)
  return [generateSinglePoint() for i in 1:n]
end

function calcOutput(p, pol)
  tmp = pol(x => p.x, y => p.y)
  println(tmp)
  return tmp
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
      default = "STDOUT"
  end

  return parse_args(s)
end

function main()  
  options = parseCommandline()
  coefficients = [randUniform(-10, 10) for i in 1:2]
  randomPolynomial = polynomial([x, y], coefficients)
  println(randomPolynomial)
  randomPoints = [Point(randUniform(-10, 10), randUniform(-10, 10), 0) for i in 1:options["amount"]]
  calcOutput.(randomPoints, randomPolynomial)
  for i in 1:options["amount"]
    randomPoints[i].z = calcOutput(randomPoints[i], randomPolynomial)
  end
  for i in 1:options["amount"]
    println(randomPoints[i].x, " ", randomPoints[i].y, " ", randomPoints[i].z)
  end
  #= toOutput = Dict("polynomial" => randomPolynomial, =#
  #=                 "points" => randomPoints) =#
  #= jsonOutput = json(toOutput) =#
  #= outputStream = open(options["output"], "w") =#
  #= write(outputStream, jsonOutput) =#
end

main()
