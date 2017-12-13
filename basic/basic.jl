#!/usr/bin/julia

struct Point
  dimensions
  point::Vector{Float64}
  Point(n) = new(n,zeros(n))
  Point(arr::Vector{Float64}) = point=arr
end



a = Point([1,2,3,4,5.2])
print(a)
