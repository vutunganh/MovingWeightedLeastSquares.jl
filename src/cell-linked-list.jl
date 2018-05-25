"""
Cell linked list splits the vector space a regular grid with edge length `EPS`.

# Member variables
- `data::Array{T, 2} where {T <: Real}`: an array holding the original data
- `grid::Array{LinkedList{T}} where {T <: Int}`: the gridded space
- `EPS::Real`: the edge length of each grid cube
- `maxs::Vector{T} where {T <: Real}`: stores the maximal coordinate over all data
- `mins::Vector{T} where {T <: Real}`: stores the minimal coordinate over all data
- `prevquery::Real`: skip `cll_neighbor_cells` in `cll_inrange` if `Int(ceil(dist / edge)) == Int(ceil(prevquery / edge))`
- `dirs::Vector{T}`: cached vector of neighbor cells
"""
struct CellLinkedList
  data::Array{T, 2} where {T <: Real}
  grid::Array{LinkedList{T}} where {T <: Int}
  EPS::Real
  maxs::Vector{T} where {T <: Real}
  mins::Vector{T} where {T <: Real}
  prevquery::Real
  dirs::Vector{Int}
end

"""
Calculates the amount of cells in the grid.
A cell linked list with the resulting amount of cells will be able to hold any vector `v` where `v .<= maxs` and `v .>= mins`.
"""
function cll_cubecount(maxs::Vector{T}, mins::Vector{T}, EPS::Real) where {T <: Real}
  size(maxs, 1) != size(mins, 1) && throw(DimensionMismatch())
  return [Int(floor(maxs[i] / EPS) - floor(mins[i] / EPS)) for i in 1:size(maxs, 1)] + 1
end

cllCubeCount(a...) = warn("`cllCubeCount` is deprecated, use `cll_cubecount` instead.")

"""
Creates a grid of empty linked lists.
"""
function cll_initgrid(cells::Vector{T}) where {T <: Integer}
  return fill(nil(T), Tuple(cells))
end

cllCreateEmptyGrid(a...) = warn("`cllCreateEmptyGrid` is deprecated, use `cll_initgrid` instead.")

"""
Finds point's index in the grid of `cll`.
"""
function cll_index(cll::CellLinkedList, pt::Point)
  return cll_index(cll.mins, cll.EPS, pt)
end

cllIndex(a...) = warn("`cllIndex` is deprecated, use `cll_index` instead.")

"""
Finds point's index in the grid given by `mins` and `EPS`.
"""
function cll_index(mins::Vector{T}, EPS::Real, pt::Point) where {T <: Real}
  return Int.(floor.((pt - mins) / EPS)) + 1
end

"""
Some operations on the CLL might increase the size of the gridded space.
Instead of recreating a new CLL, we reuse the already created lists assuming that `EPS` hasn't changed.
"""
function cll_reusegrid!(cll::CellLinkedList, newgrid::Array{LinkedList{T}}) where {T <: Integer} # FIXME: Integer?
  for l in cll.grid
    l == nil() && continue
    pos = cll_index(cll.mins, cll.EPS, cll.data[:, head(l)])
    newgrid[pos...] = l
  end
  return nothing
end

cllCreateEmptyGrid!(a...) = warn("`cllCreateEmptyGrid` is deprecated, use `cll_reusegrid` instead.")

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
  grid::Array{LinkedList{Int}, dim} = cll_initgrid(cll_cubecount(maxs, mins, EPS))

  @inbounds for i in 1:cnt
    cur = data[:, i]
    pos = cll_index(mins, EPS, cur)
    grid[pos...] = cons(i, grid[pos...])
  end

  return CellLinkedList(data, grid, EPS, maxs, mins, 0., Vector{Int}())
end

"""
Adds a point to a CLL.
"""
function cll_add!(cll::CellLinkedList, pt::Point)
  size(pt, 1) != size(cll.data, 1) && throw(DimensionMismatch())
  newmax = max.(pt, cll.maxs)
  newmin = min.(pt, cll.mins)
  if newmax != cll.maxs || newmin != cll.mins
    newgrid::Array{LinkedList{Int}, size(pt, 1)} = cll_initgrid(cll_cubecount(newmax, newmin, cll.EPS))
    cll_reusegrid!(cll, newgrid)
    cll.grid = newgrid
    cll.maxs = newmax
    cll.mins = newmin
  end
  pos = cll_index(cll, pt)
  cll.data = hcat(cll.data, pt)
  cll.grid[pos...] = cons(size(cll.data, 2), cll.grid[pos...])
  return nothing
end

cllAdd!(a...) = warn("`cllAdd!` is deprecated, use `cll_add` instead.")

"""
Removes a point from the its cell.  It remains in `data`.
"""
function cll_removeindex!(cll::CellLinkedList, idx::Integer, pos)
  h = cll.grid[pos...]
  todelete = tail(h)
  if head(h) == idx
    cll.grid[pos...] = todelete
  else 
    while todelete != nil()
      if head(todelete) == idx
        h.tail = tail(todelete)
        break
      else
        h = tail(h)
        todelete = tail(todelete)
      end
    end
  end
  return nothing
end

cllRemoveIdx!(a...) = warn("`cllRemoveIdx!` is deprecated, use `cll_removeindex!` instead.")

"""
Only removes it from the linked list, because array deletions are expensive.
Maybe add a deletion counter and perform a deletion when `counter > sqrt(length(data))`
"""
function cll_remove!(cll::CellLinkedList, idx::Integer)
  if idx < 1 || idx > size(data, 2)
    throw(BoundsError(cll.data, idx))
  end

  pos = cll.index(data[:, idx])
  cll_removeindex!(cll, idx, pos)
  return nothing
end

cllRemove!(a...) = warn("`cllRemove!` is deprecated, use `cll_remove!` instead")

"""
Changes `idx`th data to new coordinates.
"""
function cll_modify!(cll::CellLinkedList, idx::Integer, newpt::Point)
  oldPos = cll_index(cll, cll.data[:, idx])
  newPos = cll_index(cll, newpt)
  if oldPos == newPos
    return nothing
  end
  cll_removeindex!(cll, idx, oldPos)

  newmax = max.(newpt, cll.maxs)
  newmin = min.(newpt, cll.mins)
  if newmax != cll.maxs || newmin != cll.mins
    newgrid::Array{LinkedList{Int}, size(newpt, 1)} =
      cll_initgrid(cll_cubecount(newmax, newmin, cll.EPS))
    cll_reusegrid!(cll, newgrid)
    cll.grid = newgrid
    cll.maxs = newmax
    cll.mins = newmin
    newPos = cll_index(cll, newpt)
  end
  
  cll.grid[newPos...] = cons(idx, cll.grid[newPos...])
  return nothing
end

cllModify!(a...) = warn("`cllModify!` is deprecated, use `cll_modify!` instead.")

"""
Precalculates the cells, which need to be checked in cll_inrange.

### Example
`cll_neighbor_cells(e::Real, d::Real)`
"""
function cll_neighbor_cells(edge::Real, dist::Real, dim::Integer)
  dim < 1 && error("dimension has to be at least 1, but it was $dim")
  (dist < 0 || dist < 1e-9) && error("distance has to be positive, but it was $dist")
  (edge < 0 || edge < 1e-9) && error("edge length has to be positive, but it was $edge")
  r = Int(ceil(dist / edge)) + 1
  to = ntuple(i -> r, dim)
  from = .- to
  res::Vector{Vector{Int}} = []
  for c in CartesianRange(CartesianIndex(from), CartesianIndex(to))
    n = sqrt(sum(map(x -> x^2, c.I)))
    if r < n - n / 1e8
      continue
    end
    push!(res, collect(c.I))
  end
  return res
end

cllNeighborCells(a...) = warn("`cllNeighborCells` is deprecated, use `cll_neighbor_cells` instead.")

"""
Obtains the indices of `cll.data` of points, that are within `d` from `pt`.
"""
function cll_inrange(cll::CellLinkedList, pt::Point, d::Real = cll.EPS)
  if Int(ceil(d / cll.EPS)) == Int(ceil(cll.prevquery) / cll.EPS)
    return cll_inrange(cll, pt, cll.dirs, d)
  else
    neighbors = cll_neighbor_cells(cll.EPS, d, size(cll.mins, 1))
    #FIXME: cache prev query
    return cll_inrange(cll, pt, neighbors, d)
  end
end

function cll_inrange(cll::CellLinkedList, pt::Point, neighbor_dirs, d::Real = cll.EPS)
  ptCell = cll_index(cll, pt)
  res::Vector{Int} = []
  for c in neighbor_dirs
    cur = ptCell + c
    if any(cur .< 1)
      continue
    end
    if any(cur .> collect(size(cll.grid)))
      continue
    end
    for p in cll.grid[cur...]
      if norm(cll.data[:, p] - pt) < d + d / 1e8
        push!(res, p)
      end
    end
  end
  return res
end

cllInrange(a...) = warn("`cllInrange` is deprecated, use `cll_inrange` instead.")
