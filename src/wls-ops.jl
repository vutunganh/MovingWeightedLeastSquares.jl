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


function plotWls(obj::WlsObject)
  inputDim = size(obj.inputs, 2)
  println(inputDim)
  if inputDim < 1 || inputDim > 2
    error("Cannot plot $inputDim dimensional graph")
  end
  x = [linspace(-2, 2, 21) for i in 1:inputDim]
  println(x)
  println(length(x))
  if inputDim == 2
    return plot(x[1], x[2], [obj([a, b], true) for b in x[1], a in x[2]], linetype=:surface)
  else
    return plot(x[1], [obj([a], true) for a in x[1]])
  end
end