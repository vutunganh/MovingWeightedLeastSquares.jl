# Weighted least squares in 3D space.
# Approximation will be done with polynomials.
# Assume that approximated function will be in form f: R^2 -> R
# The primary purpose of this program is to find all mistakes, that will
# happen in the real program and solve them early.

pointsStr = "points"

import Base: +, -, *, /
import JSON
using DynamicPolynomials
using MultivariatePolynomials
using Base.Random

"A geometrical point structure."
struct Point
  x::Float64
  y::Float64
  z::Float64
end

"Performs an element-wise addition"
function +(p::Point, q::Point)
  return Point(p.x + q.x, p.y + q.y, p.z + q.z)
end

"Performs an element-wise subtraction"
function -(p::Point, q::Point)
  return Point(p.x - q.x, p.y - q.y, p.z - q.z)
end

"Multiplies each element by a constant"
function *(p::Point, a::Float64)
  return Point(a * p.x, a * p.y, a * p.z)
end

"Divides each element by a constant"
function /(p::Point, a::Float64)
  return Point(a / p.x, a / p.y, a / p.z)
end

"Calculates the Euclidean distance of two 3D points"
function dist(p::Point, q::Point)
  return sqrt((p.x - q.x) ^ 2 + (p.y - q.y) ^ 2 + (p.z - q.z) ^ 2)
end

"Calculates the Euclidean distance of a point from (0, 0, 0)"
function radial_distance(p::Point)
  return dist(p, Point(0,0,0))
end

"Uniformly generates a random number from a to b (b noninclusive)"
function randUniform(a = -100, b = 100)
  a + rand() * (b - a)
end

"""
Generates a random polynomial
TODO: how to make this parametrizable by max degree, 
e.g. generatePolynomial(3) -> x^1, x^2, x^3?
btw `@ polyvar a[1:n]` works
"""
function generatePolynomial()
  coefficients = [randUniform() for i in 1:6]
  @polyvar x y
  polynomial([1, x, y, x^2, x*y, y^2], coefficients)
end

"Generates x and y coordinates randomly"
function generateInput()
  Point(randUniform(), randUniform(), 0)
end

"Generates a `n` random points without outputs"
function generatePoints(n::Int64)
  [generateInput() for i in 1:n]
end

"Calculates outputs for generated points"
function calculateOutputs(points, polynomial)
  for p in points
    p.z = polynomial(x => p.x, y => p.y)
end

"An example weighting function"
function w1(d::Float64, EPS::Float64)
  1 / (d ^ 2 + EPS ^ 2)
end

function wls(approx_point::Point, points, weighting_function)
end

function program()
  rawInput = JSON.parsefile(ARGS[1])
  points = rawInput[pointsStr]
  println(generatePoint(5))
end

if length(ARGS) == 1
  program()
else
  println("Usage: 3d-version.jl <file>")
end

