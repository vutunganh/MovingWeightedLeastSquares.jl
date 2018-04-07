struct CellLinkedList{D}
  grid::Array{LinkedList{Point}, D}
  EPS::Real
  offsets::Vector{Real}
  index::Function
end

"""
Each point is in a single column.

`EPS` is the edge length of each hypercube of the cell linked list.
Also `EPS` will be used by default, when searching for neighbors.
"""
function CellLinkedList(data::Array{T, 2}, EPS::Real) where {T <: Real}
  dim = size(data, 1)
  cnt = size(data, 2)
  maxs = [maximum(data[i, :]) for i in 1:dim]
  mins = [minimum(data[i, :]) for i in 1:dim]
  cubes::Vector{Int} = [Int(ceil(maxs[i] / EPS) - floor(mins[i] / EPS)) for i in 1:dim]
  println(cubes)
  grid::Array{Nil{Point}, dim} = fill(nil(Point), Tuple(cubes))

  index = (p::Point) -> begin
    offsets = -mins
    Int.(ceil.((p + offsets) / EPS)) + 1
  end

  for i in 1:cnt
    println("===")
    cur = data[:, i]
    pos = index(cur)
    println(cur)
    println(pos)
    grid[pos...] = cons(cur, grid[pos...])
  end

  return CellLinkedList(grid, EPS, offsets, index)
end

