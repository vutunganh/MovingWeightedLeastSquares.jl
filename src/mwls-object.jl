"""
# User provided attributes
- `inputs`: the vector or array view of input points,
- `outputs`: the vector or array view of output scalars,
- `EPS::Float64`: ε of the method,
- `weightFunc::Function`: weighting function of the method.

# Automatically created attributes
If created by the `mwls` function, the attributes in this section are created automatically.

- `vars::Vector{PolyVar{true}}`: variables of the polynomial,
- `b::Vector{Monomial{true}}`: basis of the polynomial,
- `matrix`: the result of `b * transpose(b)`,
- `tree`: a k-d tree for the neareast neighbor search.
"""
struct MwlsObject
  inputs
  outputs
  EPS::Float64
  weightFunc::Function
  vars::Vector{PolyVar{true}}
  b::Vector{Monomial{true}}
  matrix
  tree
end

function MwlsObject(inputs, outputs, EPS, leafSize, weightFunc, vars, b)
  return MwlsObject(inputs, outputs, EPS, weightFunc, vars, b, b * transpose(b), KDTree(transpose(inputs); leafsize = leafSize))
end

"""
mwls(inputs::Array{Float64, 2}, outputs::Vector{Float64}, maxDegree::Int, leafSize::Int, EPS::Float64, weightFunc::Function)

Euclidean metric is used by default.

# Arguments
- `inputs`: a 2d array where each point is on a single row,
- `outputs`: an array of output scalars,
- `maxDegree::Int64`: the maximal degree of each polynomial term in the method,
- `leafSize::Int64`: the size of the leaves in the kd-tree,
- `EPS::Float64`: ε of the method,
- `weightFunc::Function`: weighting function of the method. It should be in form `(distance, EPS) -> Float64`.
"""
function mwls(inputs, outputs, maxDegree::Int, EPS::Float64, weightFunc::Function; leafSize=10)
  inputDim = size(inputs, 2)
  if inputDim == 1
    inputs = hcat(inputs, zeros(size(inputs, 1)))
  end
  vars::Vector{PolyVar{true}}, b::Vector{Monomial{true}} = generateMonomials(inputDim, maxDegree)
  return MwlsObject(inputs, outputs, EPS, leafSize, weightFunc, vars, b)
end

