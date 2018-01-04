# Weighted least squares in 3D space.
# Approximation will be done with polynomials.
# Assume that approximated function will be in form f: R^2 -> R

import Base.+
import Base.-
import Base.*
import Base./

# A geometrical point structure.
struct Point
  x::Float64
  y::Float64
  z::Float64
end

# Performs an element-wise addition
function +(p::Point, q::Point)
  return Point(p.x + q.x, p.y + q.y, p.z + q.z)
end

# Performs an element-wise subtraction
function -(p::Point, q::Point)
  return Point(p.x - q.x, p.y - q.y, p.z - q.z)
end

# Multiplies each element by a constant
function *(p::Point, a::Float64)
  return Point(a * p.x, a * p.y, a * p.z)
end

# Divides each element by a constant
function /(p::Point, a::Float64)
  return Point(a / p.x, a / p.y, a / p.z)
end

# Calculates the Euclidean distance of two 3D points
function dist(p::Point, q::Point)
  return sqrt((p.x - q.x) ^ 2 + (p.y - q.y) ^ 2 + (p.z - q.z) ^ 2)
end

function radial_distance(p::Point)
  return dist(p, Point(0,0,0))
end