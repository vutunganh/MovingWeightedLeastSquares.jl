function getInrangeData(obj::MwlsObject, inPt::Point, dist::Float64)
  return inrange(obj.tree, inPt, dist)
end

"""
    `calcMwlsCoefficients(obj::MwlsObject, inPt::Point, dist::Float64)`

Calculates the coefficients of the linear combination of the basis.
"""
function calcMwlsCoefficients(obj::MwlsObject, inPt::Point, dist::Float64)
  m = length(obj.b)
  firstTerm = zeros(m, m)
  secondTerm = zeros(m)

  data = getInrangeData(obj, inPt, dist)

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
`Dist` determines the method's distance threshold around `inPt`.
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
    `calcDiffMwlsSingle(obj::MwlsObject, pt::Point, dir::Int, times::Int; dist = obj.EPS)`

Calculates `obj`'s approximated polynomial's `times`-th derivative by `var`-th variable.
"""
function calcDiffMwlsSingle(obj::MwlsObject, pt::Point, var::Int, times::Int; dist = obj.EPS)
  c = calcMwlsCoefficients(obj, pt, dist)
  poly = polynomial(c, obj.b)
  return differentiate(poly, obj.vars[var], times)
end

"""
    `calcDiffMwlsPolys(obj::MwlsObject, pt::Point, dirs::Tuple; dist = obj.EPS)`

Calculates 
"""
function calcDiffMwlsPolys(obj::MwlsObject, pt::Point, dirs; dist = obj.EPS)
  if length(dirs) != size(obj.inputs, 2)
    error("Dimension mismatch")
  end
  c = calcMwlsCoefficients(obj, pt, dist)
  poly = polynomial(c, obj.b)
  return [differentiate(poly, obj.vars[i], dirs[i]) for i in 1:length(obj.vars)]
end

function diff(obj::MwlsObject, pt::Point, dirs::Int; dist = obj.EPS)
  diff(obj, pt, (dirs, 0); dist = dist)
end

function diff(obj::MwlsObject, pt::Point, dirs::Tuple; dist = obj.EPS)
  if size(pt, 1) == 1
    pt = [pt; 0]
  end
  pl = calcDiffMwlsPolys(obj, pt, dirs; dist = dist)
  println(pl)
  return [p(obj.vars => pt) for p in pl]
end
