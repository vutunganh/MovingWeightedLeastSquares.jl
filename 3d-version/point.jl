module Point3D

import Base: +, -, *, /, getindex

export Point, +, -, *, /, dist, radial_distance, getindex

"A geometrical point structure."
mutable struct Point
  x::Float64
  y::Float64
  z::Float64

  Point(a::Real, b::Real, c::Real) = new(a, b, c)
  Point(a::Array{<:Real, 3}) = new(a[1], a[2], a[3])
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

"Getindex for points, 1 -> x, 2 -> y, 3 -> z"
function getindex(p::Point, i::Int)
  @boundscheck if i < 1 || i > 3
    throw(BoundsError())
  end

  if i == 1
    return p.x
  elseif i == 2
    return p.y
  else
    return p.z
  end
end

end
