"""
# Attributes
- `vars::Vector{PolyVar{true}}`: variables of the polynomial, these are constructed automatically if constructed by the `wls` function
- `b::Vector{Monomial{true}}`: basis of the polynomial, these are constructed automatically if constructed by the `wls` function
- `inputs`: the vector or array view of input points
- `outputs`: the vector or array view of output scalars
- `EPS::Float64`: ε of the method
- `weightFunction::Function`: weighting function of the method
- `matrix`: the result of `b * transpose(b)`, this is constructed automatically if constructed by the `wls` function
"""
struct WlsObject
  vars::Vector{PolyVar{true}}
  b::Vector{Monomial{true}}
  inputs
  outputs
  EPS::Float64
  weightFunction
  matrix
end

function WlsObject(vars, b, inputs, outputs, EPS, wfun)
  WlsObject(vars, b, inputs, outputs, EPS, wfun, b * transpose(b))
end

"""
    wls(inputs::Vector{Point}, outputs::Vector{Float64}, maxDegree::Int64, EPS::Float64, wfun::Function)

# Arguments
- `inputs`: a 2d array where each point is on a single row,
- `outputs`: an array of output scalars,
- `maxDegree::Int64`: the maximal degree of each polynomial term in the method,
- `EPS::Float64`: ε of the method,
- `wfun::Function`: wfun should be in form `wfun(distance::Float64, EPS::Float64)`.
"""
function wls(inputs, outputs, maxDegree::Int64, EPS::Float64, wfun::Function)
  inputDim = size(inputs, 2)
  vars::Vector{PolyVar{true}}, b::Vector{Monomial{true}} = generateMonomials(inputDim, maxDegree)
  return WlsObject(vars, b, inputs, outputs, EPS, wfun)
end

