"""
Cell linked list splits the vector space a regular grid with edge length `EPS`.


# Attributes
- `data::Array{T, 2} where {T <: Real}`: an array holding the original data
- `grid::Array{LinkedList{T}} where {T <: Int}`: the gridded space
- `EPS::Real`: the edge length of each grid cube
- `maxs::Vector{T} where {T <: Real}`: stores the maximal coordinate over all data
- `mins::Vector{T} where {T <: Real}`: stores the minimal coordinate over all data
"""
mutable struct CellLinkedList
  data::Array{T, 2} where {T <: Real}
  grid::Array{LinkedList{T}} where {T <: Int}
  EPS::Real
  maxs::Vector{T} where {T <: Real}
  mins::Vector{T} where {T <: Real}
end

"""
Calculates the amount of cubes in the grid.
A cell linked list with the resulting amount of cubes will be able to hold any vector `v` where `v .<= maxs` and `v .>= mins`.
"""
function cllCubeCount(maxs::Vector{T}, mins::Vector{T}, EPS::Real) where {T <: Real}
  size(maxs, 1) != size(mins, 1) && throw(DimensionMismatch())
  return [Int(ceil(maxs[i] / EPS) - floor(mins[i] / EPS)) for i in 1:size(maxs, 1)] + 1
end

"""
Creates a grid of empty linked lists.
"""
function cllCreateEmptyGrid(cubes::Vector{Int})
  return fill(nil(Int), Tuple(cubes))
end

"""
Finds point's index in the grid of `cll`.
"""
function cllIndex(cll::CellLinkedList, pt::Point)
  return cllIndex(cll.mins, cll.EPS, pt)
end

"""
Finds point's index in the grid given by `mins` and `EPS`.
"""
function cllIndex(mins::Vector{T}, EPS::Real, pt::Point) where {T <: Real}
  return Int.(floor.((pt - mins) / EPS)) + 1
end

"""
Some operations on the CLL might increase the size of the gridded space.
Instead of recreating a new CLL, we reuse the already created lists assuming that `EPS` hasn't changed.
"""
function cllRecreateGrid!(cll::CellLinkedList,
                          newGrid::Array{LinkedList{T}}) where {T <: Int}
  for l in cll.grid
    l == nil() && continue

    pos = cllIndex(cll.mins, cll.EPS, cll.data[:, head(l)])
    newGrid[pos...] = l
  end
  return nothing
end

"""
Each point is in a single column.

# Arguments
- `data::Array{T} where {T <: Real}`: a 2d array with data, where each vector is in a single column,
- `EPS::Real`: the edge length of the regular grid.
"""
function CellLinkedList(data::Array{T, 2}, EPS::Real) where {T <: Real}
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

"""
Adds a point to a CLL.
"""
function cllAdd!(cll::CellLinkedList, pt::Point)
  size(pt, 1) != size(cll.data, 1) && throw(DimensionMismatch())
  nMax = max.(pt, cll.maxs)
  nMin = min.(pt, cll.mins)
  if nMax != cll.maxs || nMin != cll.mins
    newGrid::Array{LinkedList{Int}, size(pt, 1)} = cllCreateEmptyGrid(cllCubeCount(nMax, nMin, cll.EPS))
    cllRecreateGrid!(cll, newGrid)
    cll.grid = newGrid
    cll.maxs = nMax
    cll.mins = nMin
  end
  pos = cllIndex(cll, pt)
  cll.data = hcat(cll.data, pt)
  cll.grid[pos...] = cons(size(cll.data, 2), cll.grid[pos...])
  return nothing
end

"""
Removes a point from the its cell.  It remains in `data`.
"""
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
  return nothing
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
  return nothing
end

"""
Changes `idx`th data to new coordinates.
"""
function cllModify!(cll::CellLinkedList, idx::Int, newPt::Point)
  oldPos = cllIndex(cll, cll.data[:, idx])
  newPos = cllIndex(cll, newPt)
  if oldPos == newPos
    return nothing
  end
  cllRemoveIdx(cll, idx, oldPos)

  nMax = max.(newPt, cll.maxs)
  nMin = min.(newPt, cll.mins)
  if nMax != cll.maxs || nMin != cll.mins
    newGrid::Array{LinkedList{Int}, size(newPt, 1)} =
      cllCreateEmptyGrid(cllCubeCount(nMax, nMin, cll.EPS))
    cllRecreateGrid!(cll, newGrid)
    cll.grid = newGrid
    cll.maxs = nMax
    cll.mins = nMin
    newPos = cllIndex(cll, newPt)
  end
  
  cll.grid[newPos...] = cons(idx, cll.grid[newPos...])
  return nothing
end

"""
Obtains the indices of `cll.data` of points, that are within `d` from `pt`.
"""
function cllInrange(cll::CellLinkedList, pt::Point, d::Real = cll.EPS)
  pos = Tuple(cllIndex(cll, pt))
  cnt::Int = ceil(d / cll.EPS)
  from = CartesianIndex(pos .- cnt)
  to = CartesianIndex(pos .+ cnt)
  res::Vector{Int} = []
  oneT = pos .- pos .+ 1

  for c in CartesianRange(from, to)
    if true in (size(cll.grid) .< c.I)
      continue
    end
    if true in (c.I .< oneT)
      continue
    end
    for p in cll.grid[c]
      if norm(cll.data[:, p] - pt) < d + d / 1e6 # <=
        push!(res, p)
      end
    end
  end
  
  return res
end

