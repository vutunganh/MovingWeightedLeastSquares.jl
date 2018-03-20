"""
    `(wo::WlsObject)(inputPoint::Point, value::Bool)`

Approximates the `WlsObject` at `inputPoint`.
If `value` is true, the approximated value will be outputted instead of the coefficients.
"""
function (obj::WlsObject)(inPt::Point, value::Bool)
  m = length(obj.b)
  datalen = size(obj.inputs, 1)
  firstTerm = zeros(m, m)
  secondTerm = zeros(m)

  for p in 1:datalen
    curPt = view(obj.inputs, p, :)
    w = obj.weightFunction(norm(inPt - curPt), obj.EPS)
    for i in 1:m
      for j in 1:m
        firstTerm[i, j] += w * obj.matrix[i, j](obj.vars => curPt)
      end
      secondTerm[i] += w * obj.b[i](obj.vars => curPt) * obj.outputs[p]
    end
  end

  coefficients = firstTerm \ secondTerm
  if value
    poly = polynomial(coefficients, obj.b)
    return poly(obj.vars => inPt)
  else
    return coefficients
  end
end

