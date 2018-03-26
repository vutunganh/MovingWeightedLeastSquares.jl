function getInrangeData(obj::MwlsObject, inPt::Point, dist::Float64 = obj.EPS)
  return inrange(obj.tree, inPt, dist)
end

"""
    `calcMwlsCoefficients(obj::MwlsObject, inPt::Point, dist::Float64)`

Calculates the coefficients of the linear combination of the basis.
"""
function calcMwlsCoefficients(obj::MwlsObject, inPt::Point, dist::Float64 = obj.EPS)
  m = length(obj.b)
  firstTerm = zeros(m, m)
  secondTerm = zeros(m)

  data = getInrangeData(obj, inPt, dist)
  if length(data) == 0
    return secondTerm
  end

  for p in data
    curPt = obj.inputs[p, :]
    w = obj.weightFunc(norm(inPt - curPt), obj.EPS)
    for i in 1:m
      for j in 1:m
        firstTerm[i, j] += w * obj.matrix[i, j](obj.vars => curPt)
      end
      secondTerm[i] += w * obj.b[i](obj.vars => curPt) * obj.outputs[p]
    end
  end

  result = firstTerm \ secondTerm
  return result
end

"""
    `(obj::MwlsObject)(inputPoint::Point, dist::Float64)`

Approximates the `MwlsObject` at `inputPoint`.
`dist` determines the method's distance threshold around `inPt`.
"""
function (obj::MwlsObject)(inPt::Point, dist = obj.EPS)
  if size(inPt, 1) == 1
    inPt = [inPt; 0]
  end
  c = calcMwlsCoefficients(obj, inPt, dist)
  poly = polynomial(c, obj.b)
  return poly(obj.vars => inPt)
end

"""
    `calcDiffMwlsPolys(obj::MwlsObject, inPt::Point, dirs::NTuple{N, Int64}; dist = obj.EPS) where {N}`

Calculates the differentiated polynomials in each direction.
"""
function calcDiffMwlsPolys(obj::MwlsObject, inPt::Point, dirs::NTuple{N, Int64}; dist = obj.EPS) where {N}
  c = calcMwlsCoefficients(obj, inPt, dist)
  poly = polynomial(c, obj.b)
  return [differentiate(poly, obj.vars[i], dirs[i]) for i in 1:length(obj.vars)]
end

"this function exists, because writing a tuple literal with a single element is difficult"
function diff(obj::MwlsObject, inPt::Point, dirs::Int64; dist = obj.EPS)
  diff(obj, inPt, Tuple(dirs); dist = dist)
end

"""
    `diff(obj::MwlsObject, inPt::Point, dirs::NTuple{N, Int64}; dist = obj.EPS)`

Calculates the approximated derivative at `inPt`, where `x[i]` is differentiated `dirs[i]` times.
"""
function diff(obj::MwlsObject, inPt::Point, dirs::NTuple{N, Int64}; dist = obj.EPS) where {N}
  if length(inPt) == 1
    inPt = [inPt; 0]
  end
  N != length(obj.vars) && error("Mismatch between tuple size and amount of variables")
  pl = calcDiffMwlsPolys(obj, inPt, dirs; dist = dist)
  return [p(obj.vars => inPt) for p in pl]
end
