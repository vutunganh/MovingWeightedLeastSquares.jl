"""
Return the indices of sample inputs in `obj` within `dist` from `inPt`.
"""
function getInrangeData(obj::MwlsNaiveObject, inPt::Point, dist::Real = obj.EPS)
  res::Vector{Int} = []

  @inbounds for p in 1:size(obj.inputs, 2)
    if norm(obj.inputs[:, p] - inPt) < dist + dist / 1e8
      push!(res, p)
    end
  end

  return res
end

"""
Returns the indices of sample inputs in `obj::MwlsKdObject` within `dist` from `inPt`.
"""
function getInrangeData(obj::MwlsKdObject, inPt::Point, dist::Real = obj.EPS)
  return inrange(obj.tree, inPt, dist)
end

"""
Returns the indices of sample inputs in `obj::MwlsCllObject` within `dist` from `inPt`.
"""
function getInrangeData(obj::MwlsCllObject, inPt::Point, dist::Real = obj.EPS)
  return cllInrange(obj.cll, inPt, dist)
end

"""
    `calcMwlsCoefficients(obj::MwlsObject, inPt::Point, dist::Float64)`

Calculates the coefficients of the linear combination of basis functions for each dimension of output vectors.
"""
function calcMwlsCoefficients(obj::MwlsObject, inPt::Point, dist::Real = obj.EPS)
  m = size(obj.b, 1)
  outputDim = size(obj.outputs, 1)
  firstTerm = zeros(m, m)
  secondTerm = zeros(outputDim, m)

  data = getInrangeData(obj, inPt, dist)

  for p in data
    curPt = obj.inputs[:, p]
    w = obj.weightFunc(norm(inPt - curPt), obj.EPS)
    @inbounds for j in 1:m
      @inbounds for i in 1:m
        firstTerm[i, j] += w * obj.matrix[i, j](obj.vars => curPt)
      end
      secondTerm[:, j] += w * obj.b[j](obj.vars => curPt) * obj.outputs[:, p]
    end
  end

  # probably singular matrix
  # returns zeroes
  if cond(firstTerm, 1) >= 1e14 || cond(firstTerm, 2) >= 1e14
    return zeros(m, outputDim)
  end

  result = firstTerm \ transpose(secondTerm)
  return result
end

"""
    `approximate(obj::MwlsObject, inPt::Point, dist::Real = obj.EPS)`

Approximates the `MwlsObject` at `inputPoint`.
If for each sample input ``x_i``, ``\|\text{inPt} - x_i\|`` is greater than `dist`, then the result of the weight function is zero.
"""
function approximate(obj::MwlsObject, inPt::Point, dist::Real = obj.EPS)
  cs = calcMwlsCoefficients(obj, inPt, dist)
  poly = [polynomial(cs[:, i], obj.b) for i in 1:size(cs, 2)]
  res = [p(obj.vars => inPt) for p in poly]
  return length(res) == 1 ? res[1] : res
end

function (obj::MwlsNaiveObject)(inPt::Point, dist::Real = obj.EPS)
  approximate(obj, inPt, dist)
end

function (obj::MwlsNaiveObject)(inPt::Real, dist::Real = obj.EPS)
  return obj([inPt], dist)
end

function (obj::MwlsKdObject)(inPt::Point, dist::Real = obj.EPS)
  approximate(obj, inPt, dist)
end

function (obj::MwlsKdObject)(inPt::Real, dist::Real = obj.EPS)
  return obj([inPt, 0], dist)
end

function (obj::MwlsCllObject)(inPt::Point, dist::Real = obj.EPS)
  approximate(obj, inPt, dist)
end

function (obj::MwlsCllObject)(inPt::Real, dist::Real = obj.EPS)
  return obj([inPt, 0], dist)
end

"""
    `calcDiffMwlsPolys(obj::MwlsObject, inPt::Point, dirs::NTuple{N, Int64}; dist::Real = obj.EPS) where {N}`
"""
function calcDiffMwlsPolys(obj::MwlsObject, inPt::Point, dirs::NTuple{N, Int64}; dist::Real = obj.EPS) where {N}
  cs = calcMwlsCoefficients(obj, inPt, dist)
  poly = [polynomial(cs[:, i], obj.b) for i in 1:size(cs, 2)]
  @inbounds for j in 1:length(poly)
    @inbounds for i in 1:length(dirs)
      poly[j] = differentiate(poly[j], obj.vars[i], dirs[i])
    end
  end
  return poly
end

function mwlsDiff(obj::MwlsNaiveObject, inPt::Real, dirs::Int; dist::Real = obj.EPS)
  mwlsDiff(obj, [inPt], ntuple(x -> dirs, 1); dist = dist)
end

function mwlsDiff(obj::MwlsNaiveObject, inPt::Point, dirs::Int; dist::Real = obj.EPS)
  mwlsDiff(obj, inPt, ntuple(x -> dirs, 1); dist = dist)
end

function mwlsDiff(obj::MwlsNaiveObject, inPt::Real, dirs::NTuple{N, Int}; dist::Real = obj.EPS) where {N}
  mwlsDiff(obj, [inPt], dirs, dist)
end

function mwlsDiff(obj::MwlsObject, inPt::Real, dirs::Int; dist::Real = obj.EPS)
  mwlsDiff(obj, [inPt, 0], ntuple(x -> dirs, 1); dist = dist)
end

"this function exists, because writing a 1d point is cumbersome"
function mwlsDiff(obj::MwlsObject, inPt::Real, dirs::NTuple{N, Int}; dist::Real = obj.EPS) where {N}
  mwlsDiff(obj, [inPt, 0], ntuple(x -> dirs, 1); dist = dist)
end

"this function exists, because writing a tuple literal with a single element is difficult"
function mwlsDiff(obj::MwlsObject, inPt::Point, dirs::Int; dist::Real = obj.EPS)
  mwlsDiff(obj, inPt, ntuple(x -> dirs, 1); dist = dist)
end

"""
    `mwlsDiff(obj::MwlsObject, inPt::Point, dirs::NTuple{N, Int64}; dist = obj.EPS)`

Calculates the approximated derivative at `inPt`, where `x[i]` is differentiated `dirs[i]` times.
"""
function mwlsDiff(obj::MwlsObject, inPt::Point, dirs::NTuple{N, Int}; dist::Real = obj.EPS) where {N}
  N != length(obj.vars) && error("Mismatch between tuple size and amount of variables")
  pl = calcDiffMwlsPolys(obj, inPt, dirs; dist = dist)
  res = [p(obj.vars => inPt) for p in pl]
  return length(res) == 1 ? res[1] : res
end
