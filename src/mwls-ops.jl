"""
return the indices of points in range for a MwlsNaiveObject
"""
function getInrangeData(obj::MwlsNaiveObject, inPt::Point, dist::Real = obj.EPS)
  res::Vector{Int} = []

  for p in 1:size(obj.inputs, 2)
    if norm(obj.inputs[:, p] - inPt) < dist + dist * 1/e-6
      push!(res, p)
    end
  end

  return res
end

"""
returns the indices of points in range for a MwlsKdObject
"""
function getInrangeData(obj::MwlsKdObject, inPt::Point, dist::Real = obj.EPS)
  return inrange(obj.tree, inPt, dist)
end

"""
returns the indices of points in range for a MwlsCllObject
"""
function getInrangeData(obj::MwlsCllObject, inPt::Point, dist::Real = obj.EPS)
  return cllInrange(obj.cll, inPt, dist)
end

"""
    `calcMwlsCoefficients(obj::MwlsObject, inPt::Point, dist::Float64)`

Calculates the coefficients of the linear combination of the basis for each vector of outputs.
"""
function calcMwlsCoefficients(obj::MwlsObject, inPt::Point, dist::Real = obj.EPS)
  m = size(obj.b, 1)
  outputDim = size(obj.outputs, 1)
  firstTerm = zeros(m, m)
  secondTerm = zeros(outputDim, m)

  data = getInrangeData(obj, inPt, dist)
  if length(data) < length(obj.b)
    return transpose(secondTerm)
  end

  for p in data
    curPt = obj.inputs[:, p]
    w = obj.weightFunc(norm(inPt - curPt), obj.EPS)
    for j in 1:m
      for i in 1:m
        firstTerm[i, j] += w * obj.matrix[i, j](obj.vars => curPt)
      end
      tmp = w * obj.b[j](obj.vars => curPt) * obj.outputs[:, p]
      secondTerm[:, j] += tmp
    end
  end

  result = firstTerm \ transpose(secondTerm)
  return result
end

function approximate(obj::MwlsObject, inPt::Point, dist = obj.EPS)
  cs = calcMwlsCoefficients(obj, inPt, dist)
  poly = [polynomial(cs[:, i], obj.b) for i in 1:size(cs, 2)]
  res = [p(obj.vars => inPt) for p in poly]
  return length(res) == 1 ? res[1] : res
end

function (obj::MwlsNaiveObject)(inPt::Point, dist = obj.EPS)
  approximate(obj, inPt, dist)
end

"""
    `(obj::MwlsObject)(inputPoint::Point, dist::Float64)`

Approximates the `MwlsObject` at `inputPoint`.
`dist` determines the method's distance threshold around `inPt`.
"""
function (obj::MwlsKdObject)(inPt::Point, dist = obj.EPS)
  approximate(obj, inPt, dist)
end

function (obj::MwlsKdObject)(inPt::Real, dist = obj.EPS)
  return obj([inPt, 0], dist)
end

"""
    `(obj::MwlsObject)(inputPoint::Point, dist::Float64)`

Approximates the `MwlsObject` at `inputPoint`.
`dist` determines the method's distance threshold around `inPt`.
"""
function (obj::MwlsCllObject)(inPt::Point, dist = obj.EPS)
  approximate(obj, inPt, dist)
end

function (obj::MwlsCllObject)(inPt::Real, dist = obj.EPS)
  return obj([inPt, 0], dist)
end

"""
    `calcDiffMwlsPolys(obj::MwlsObject, inPt::Point, dirs::NTuple{N, Int64}; dist = obj.EPS) where {N}`

Calculates the differentiated polynomials in each direction.
"""
function calcDiffMwlsPolys(obj::MwlsObject, inPt::Point, dirs::NTuple{N, Int64}; dist = obj.EPS) where {N}
  cs = calcMwlsCoefficients(obj, inPt, dist)
  poly = [polynomial(cs[:, i], obj.b) for i in 1:size(cs, 2)]
  for j in 1:length(poly)
    for i in 1:length(dirs)
      poly[j] = differentiate(poly[j], obj.vars[i], dirs[i])
    end
  end
  return poly
end

"this function exists, because writing a 1d point is cumbersome"
function diff(obj::MwlsObject, inPt::Real, dirs::Int; dist = obj.EPS)
  diff(obj, [inPt, 0], Tuple(dirs); dist = dist)
end

"this function exists, because writing a tuple literal with a single element is difficult"
function diff(obj::MwlsObject, inPt::Point, dirs::Int; dist = obj.EPS)
  diff(obj, inPt, Tuple(dirs); dist = dist)
end

"""
    `diff(obj::MwlsObject, inPt::Point, dirs::NTuple{N, Int64}; dist = obj.EPS)`

Calculates the approximated derivative at `inPt`, where `x[i]` is differentiated `dirs[i]` times.
"""
function diff(obj::MwlsObject, inPt::Point, dirs::NTuple{N, Int}; dist = obj.EPS) where {N}
  N != length(obj.vars) && error("Mismatch between tuple size and amount of variables")
  pl = calcDiffMwlsPolys(obj, inPt, dirs; dist = dist)
  res = [p(obj.vars => inPt) for p in pl]
  return length(res) == 1 ? res[1] : res
end
