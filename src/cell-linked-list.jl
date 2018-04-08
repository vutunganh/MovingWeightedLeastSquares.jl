const GridType = Array{LinkedList{T}} where {T <: Int}

struct CellLinkedList{D}
  data::Array{T, 2} where {T <: Real}
  grid::Array{LinkedList{T}, D} where {T <: Int}
  EPS::Real
  maxs::Vector{T} where {T <: Real}
  mins::Vector{T} where {T <: Real}
end

function cllCubeCount(maxs::Vector{T}, mins::Vector{T}, EPS::Real) where {T <: Real}
  size(maxs, 1) != size(mins, 1) && throw(DimensionMismatch())
  return [Int(ceil(maxs[i] / EPS) - floor(mins[i] / EPS)) for i in 1:size(maxs, 1)]
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

function cllRecreateGrid!(oldGrid::GridType,
                          newGrid::GridType,
                          mins::Vector{T},
                          EPS::Real) where {T <: Real}
  for l in oldGrid
    l == nil() && continue

    pos = cllIndex(mins, EPS, head(l))
    newGrid[pos...] = l
  end
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
    pos = cllIndex(mins, EPS, cur)
    grid[pos...] = cons(i, grid[pos...])
  end

  return CellLinkedList(data, grid, EPS, maxs, mins)
end

function cllAdd!(cll::CellLinkedList, pt::Point)
  size(pt, 1) != size(cll.data, 1) && throw(DimensionMismatch())
  nMax = max.(pt, cll.maxs)
  nMin = min.(pt, cll.mins)
  if nMax != cll.maxs || nMin != cll.mins
    newGrid::GridType = cllCreateEmptyGrid(cllCubeCount(nMax, nMin, cll.EPS))
    cllRecreateGrid!(cll.grid, newGrid, nMin, cll.EPS)
    cll.grid = newGrid
    cll.maxs = nMax
    cll.mins = nMin
  end
  pos = cll.index(pt)
  hcat(cll.data, pt)
  cll.grid[pos...] = cons(size(cll.data, 2), cll.grid[pos...])
end

function cllRemoveIdx(cll::CellLinkedList, idx::Int, pos)
  h = cll.grid[pos...]
  toDel = tail(h)
  if head(h) == idx
    cll.grid[pos...] = toDel
  else 
    while toDel != nil()
      if head(toDel) == idx
        h.tail = tail(toDel)
        break
      else
        h = tail(h)
        toDel = tail(toDel)
      end
    end
  end
end

"""
Only removes it from the linked list, because array deletions are expensive.
Maybe add a deletion counter and perform a deletion when `counter > sqrt(length(data))`
"""
function cllRemove!(cll::CellLinkedList, idx::Int)
  if idx < 1 || idx > size(data, 2)
    throw(BoundsError(cll.data, idx))
  end

  pos = cll.index(data[:, idx])
  cllRemoveIdx(cll, idx, pos)
end

function cllModify!(cll::CellLinkedList, idx::Int, newPt::Point)
  oldPos = cllIndex(cll, data[:, idx])
  if oldPos == newPos
    return nothing
  end
  cllRemoveIdx(cll, idx, oldPos)

  nMax = max.(newPt, cll.maxs)
  nMin = min.(newPt, cll.mins)
  if nMax != cll.maxs || nMin != cll.mins
    newGrid::GridType = cllCreateEmptyGrid(cllCubeCount(nMax, nMin, cll.EPS))
    cllRecreateGrid!(cll.grid, newGrid, nMin, cll.EPS)
    cll.grid = newGrid
    cll.maxs = nMax
    cll.mins = nMin
  end
  newPos = cllIndex(cll, newPt)
  
  cll.grid[newPos...] = cons(idx, cll.grid[newPos...])
end

function cllInrange(cll::CellLinkedList, pt::Point, d::Real = cll.EPS)
  # TODO: pt is outside grid
  pos = Tuple(cllIndex(cll, pt))
  cnt = d / cll.EPS
  from = CartesianIndex(pos) - cnt
  to = CartesianIndex(pos) + cnt
  res::Vector{Int} = []
  for c in CartesianRange(from, to)
    for p in c
      if norm(data[p] - pt) < d + 1e-9
        push!(res, p)
      end
    end
  end
  
  return res
end

