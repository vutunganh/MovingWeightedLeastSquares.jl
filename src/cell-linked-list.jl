const GridType = Array{LinkedList{T}} where {T <: Int}

struct CellLinkedList{D}
  data::Array{T, 2} where {T <: Real}
  grid::Array{LinkedList{T}, D} where {T <: Int}
  EPS::Real
  maxs::Vector{T} where {T <: Real}
  mins::Vector{T} where {T <: Real}
end

function cllCubeCount(maxs::Vector{Real}, mins::Vector{Real}, EPS::Real)
  return [Int(ceil(maxs[i] / EPS) - floor(mins[i] / EPS)) for i in 1:dim]
end

function cllCreateEmptyGrid(cubes::Vector{Int})
  return fill(nil(Int), Tuple(cubes))
end

function cllIndex(cll::CellLinkedList, pt::Point)
  return Int.(floor.((pt - cll.mins) / cll.EPS)) + 1
end

function cllIndex(mins::Vector{T}, EPS::Real, pt::Point) where {T <: Real}
  return Int.(floor.((pt - mins) / EPS)) + 1
end

function cllRecreateGrid!(oldGrid::GridType, newGrid::GridType, mins::Vector{T}, EPS::Real) where {T <: Real}
  for l in oldGrid
    l == nil() && continue

    pos = cllIndex(mins, EPS, head(l))
    newGrid[pos...] = l
  end
end

function add!(cll::CellLinkedList, pt::Point)
  size(pt, 1) != size(cll.data, 1) && error("Dimension mismatch")
  nMax = max.(pt, cll.maxs)
  nMin = min.(pt, cll.mins)
  if nMax != cll.maxs || nMin != cll.mins
    nCubes::Vector{Int} = [Int]
  end
  pos = cll.index(pt)
  hcat(cll.data, pt)
  cll.grid[pos...] = cons(size(cll.data, 2), cll.grid[pos...])
end

"""
Each point is in a single column.

`EPS` is the edge length of each hypercube of the cell linked list.
Also `EPS` will be used by default, when searching for neighbors.
"""
function CellLinkedList(data::Array{T, 2}, EPS::T) where {T <: Real}
  dim = size(data, 1)
  cnt = size(data, 2)
  maxs = [maximum(data[i, :]) for i in 1:dim]
  mins = [minimum(data[i, :]) for i in 1:dim]
  grid::Array{LinkedList{Int}, dim} = cllCreateEmptyGrid(cllCubeCount(maxs, mins, EPS))

  for i in 1:cnt
    cur = data[:, i]
    pos = index(cur)
    grid[pos...] = cons(i, grid[pos...])
  end

  return CellLinkedList(data, grid, EPS, index, maxs, mins)
end

function remove(cll::CellLinkedList, idx::Int)
  if idx < 1 || idx > size(data, 2)
    throw(BoundsError(cll.data, idx))
  end

  pos = cll.index(data[:, idx])
  cur = cll.grid[pos...]
  while cur != nil()
    if cur == idx

  end
end

