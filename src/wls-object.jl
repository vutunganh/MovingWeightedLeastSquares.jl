"""
# User provided attributes
- `inputs`: the vector or array view of input points
- `outputs`: the vector or array view of output scalars
- `EPS::Float64`: ε of the method
- `weightFunc::Function`: weighting function of the method

# Automatically created attributes
If created by the `wls` function, the attributes in this section are created automatically.

- `vars::Vector{PolyVar{true}}`: variables of the polynomial,
- `b::Vector{Monomial{true}}`: basis of the polynomial,
- `matrix`: the result of `b * transpose(b)`.
"""
struct WlsObject
  inputs
  outputs
  EPS::Float64
  weightFunc::Function
  vars::Vector{PolyVar{true}}
  b::Vector{Monomial{true}}
  matrix
end

function WlsObject(inputs, outputs, EPS, weightFunc, vars, b)
  WlsObject(inputs, outputs, EPS, weightFunc, vars, b, b * transpose(b))
end

"""
    wls(inputs::Array{Float64, 2}, outputs::Vector{Float64}, maxDegree::Int, EPS::Float64, weightFunc::Function)

# Arguments
- `inputs`: a 2d array where each point is on a single row,
- `outputs`: an array of output scalars,
- `maxDegree::Int64`: the maximal degree of each polynomial term in the method,
- `EPS::Float64`: ε of the method,
- `weightFunc::Function`: weighting function of the method. It should be in form of `(distance, EPS) -> Float64`.
"""
function wls(inputs, outputs, maxDegree::Int, EPS::Float64, weightFunc::Function)
  inputDim = size(inputs, 2)
  vars::Vector{PolyVar{true}}, b::Vector{Monomial{true}} = generateMonomials(inputDim, maxDegree)
  return WlsObject(inputs, outputs, EPS, weightFunc, vars, b)
end

