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
    `(mo::MwlsObject)(inputPoint::Point, dist::Float64)`

Approximates the `MwlsObject` at `inputPoint`.
`Dist` determines the method's distance threshold around `inPt`.
"""
function (obj::MwlsObject)(inPt::Point, dist::Float64)
  if size(inPt, 1) == 1
    inPt = [inPt; 0]
  end
  c = calcMwlsCoefficients(obj, inPt, dist)
  poly = polynomial(c, obj.b)
  return poly(obj.vars => inPt)
end
