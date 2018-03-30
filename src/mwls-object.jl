"""
# User provided attributes
- `inputs::Array{Float64, 2}`: a 2d array of input points where each point is on a single row,
- `outputs::Vector{Float64}`: a vector of output scalars,
- `EPS::Float64`: ε of the method (default distance threshold for neighbor search),
- `weightFunc::Function`: weighting function of the method.

# Automatically created attributes
If created by the `mwls` function, the attributes in this section are created automatically.

- `vars::Vector{PolyVar{true}}`: variables of the polynomial,
- `b::Vector{Monomial{true}}`: the basis of the polynomial,
- `matrix`: the result of `b * transpose(b)`,
- `tree`: a k-d tree for neareast neighbor search.
"""
struct MwlsObject
  inputs::Array{Float64, 2}
  outputs::Vector{Float64}
  EPS::Float64
  weightFunc::Function
  vars::Vector{PolyVar{true}}
  b::Vector{Monomial{true}}
  matrix::Array{Term{true, Int64}, 2}
  tree::KDTree
end

function MwlsObject(inputs, outputs, EPS, weightFunc, vars, b; leafSize = 10)
  return MwlsObject(inputs, outputs, EPS, weightFunc, vars, b, b * transpose(b), KDTree(transpose(inputs); leafsize = leafSize))
end

"""
mwls(inputs::Array{Float64, 2}, outputs::Vector{Float64}, maxDegree::Int, leafSize::Int, EPS::Float64, weightFunc::Function)

Euclidean metric is used by default.

# Arguments
- `inputs`: a 2d array of input points where each point is on a single row,
- `outputs`: a vector of output scalars,
- `EPS::Float64`: ε of the method (default distance threshold for neighbor search),
- `weightFunc::Function`: weighting function of the method. It should be in form `(distance, EPS) -> Float64`.

# Keyword arguments
- `leafSize::Int64`: the size of the leaves in the kd-tree, 10 by default.
- `maxDegree::Int64`: the maximal degree of each polynomial term in the method, 2 by default.
"""
function mwls(inputs, outputs, EPS::Float64, weightFunc::Function; leafSize::Int = 10, maxDegree::Int = 10)
  length(outputs) != size(inputs, 1) && error("The amount of inputs and outputs differs.")
  inputDim = size(inputs, 2)
  if inputDim == 1
    inputs = hcat(inputs, zeros(size(inputs, 1)))
  end
  vars::Vector{PolyVar{true}}, b::Vector{Monomial{true}} = generateMonomials(inputDim, maxDegree)
  return MwlsObject(inputs, outputs, EPS, weightFunc, vars, b; leafSize = leafSize)
end

