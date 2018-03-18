"""
# Attributes
- `vars::Vector{PolyVar{true}}`: variables of the polynomial, these are constructed automatically if constructed by the `wls` function
- `b::Vector{Monomial{true}}`: basis of the polynomial, these are constructed automatically if constructed by the `wls` function
- `inputs::Vector{Point}`: the vector of input points
- `outputs::Vector{Float64}`: the vector of output scalars
- `inputDimension::Int`: dimension of each input in `inputs`, this is constructed automatically if constructed by the `wls` function
- `EPS::Float64`: ε of the method
- `weightFunction::Function`: weighting function of the method
- `matrix`: the result of `b * transpose(b)`, this is constructed automatically if constructed by the `wls` function
"""
struct WlsObject
  vars::Vector{PolyVar{true}}
  b::Vector{Monomial{true}}
  inputs::Vector{Point}
  outputs::Vector{Float64}
  inputDimension::Int
  EPS::Float64
  weightFunction
  matrix
end

function WlsObject(vars, b, inputs, outputs, EPS, wfun)
  WlsObject(vars, b, inputs, outputs, length(inputs[1]) - 1, EPS, wfun, b * transpose(b))
end

"""
    wls(inputs::Vector{Point}, outputs::Vector{Float64}, maxDegree::Int64, EPS::Float64, wfun::Function)

# Arguments
- `inputs::Vector{Point}`: a vector of input points,
- `outputs::Vector{Point}`: a vector of output points,
- `maxDegree::Int64`: the maximal degree of each polynomial term in the method,
- `EPS::Float64`: ε of the method,
- `wfun::Function`: wfun should be in form `wfun(distance::Float64, EPS::Float64)`.
"""
function wls(inputs::Vector{Point}, outputs::Vector{Float64}, maxDegree::Int64, EPS::Float64, wfun::Function)
  inputDimension = length(inputs[1])
  vars::Vector{PolyVar{true}}, b::Vector{Monomial{true}} = generateMonomials(inputDimension, maxDegree)
  return WlsObject(vars, b, inputs, outputs, EPS, wfun)
end

