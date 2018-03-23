function getInrangeData(obj::MwlsObject, inPt::Point, dist::Float64)
  return inrange(obj.tree, inPt, dist)
end


function getMwlsCoefficients(obj::MwlsObject, inPt::Point, dist::Float64)
  m = length(obj.b)
  firstTerm = zeros(m, m)
  secondTerm = zeros(m)

  data = getInrangeData(obj, inPt, dist)

  for p in data
    println(obj.inputs[p, :])
    curPt = obj.inputs[p, :]
    w = obj.weightFunction(norm(inPt - curPt), obj.EPS)
    for i in 1:m
      for j in 1:m
        firstTerm[i, j] += w * obj.matrix[i, j](obj.vars => curPt)
      end
      secondTerm[i] += w * obj.b[i](obj.vars => curPt) * obj.outputs[p]
    end
  end

  println(firstTerm)
  println(secondTerm)
  result = firstTerm \ secondTerm
  return result
end

"""
    `(mo::MwlsObject)(inputPoint::Point, dist::Float64)`

Approximates the `MwlsObject` at `inputPoint`.
`Dist` determines the method's distance threshold around `inPt`.
"""
function (obj::MwlsObject)(inPt::Point, dist::Float64)
  data = inrange(obj.tree, inPt, dist)

  c = getMwlsCoefficients(obj, inPt, dist)
  poly = polynomial(c, obj.b)
  return poly(obj.vars => inPt)
end
