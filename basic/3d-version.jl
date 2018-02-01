# Weighted least squares in 3D space.
# Approximation will be done with polynomials.
# Assume that approximated function will be in form f: R^2 -> R
# The primary purpose of this program is to find all mistakes, that will
# happen in the real program.

pointsStr = "points"

import Base: +, -, *, /
import JSON

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

function radial_distance(p::Point)
  return dist(p, Point(0,0,0))
end

"""
Generates the polynomial basis
This should be much simpler in the N-dimensional case
dimensions will be indexed by a number
such as x -> 0, y -> 1, ...
"""
function gen_base()
  return (
          1,
          (p::Point -> p.x),
          (p::Point -> p.y),
          (p::Point -> p.x ^ 2),
          (p::Point -> p.x * p.y),
          (p::Point -> p.y ^ 2)
 )
end

"An example weighting function"
function w1(d::Float64, EPS::Float64)
  1 / (d ^ 2 + EPS ^ 2)
end

function wls(approx_point::Point, points, weighting_function)
end

function program()
  print(ARGS)
  rawInput = JSON.parsefile(ARGS[1])
  points = rawInput[pointsStr]

end

if length(ARGS) == 1
  program()
else
  println("Usage: 3d-version.jl <file>")
end

