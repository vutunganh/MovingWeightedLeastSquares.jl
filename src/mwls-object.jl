"""
# Attributes
- `vars::Vector{PolyVar{true}}`: variables of the polynomial, created automatically if constructed by the `mwls` function,
- `b::Vector{Monomial{true}}`: basis of the polynomial, created automatically if constructed by the `mwls` function,
- `inputs`: the vector or array view of input points,
- `outputs`: the vector or array view of output scalars,
- `EPS::Float64`: ε of the method,
- `weightFunction::Function`: weighting function of the method,
- `matrix`: the result of `b * transpose(b)`, created automatically if constructed by the `mwls` function,
- `tree`: a k-d tree for the neareast neighbor search, created automatically if constructed by `mwls` function.
"""
struct MwlsObject
  vars::Vector{PolyVar{true}}
  b::Vector{Monomial{true}}
  inputs
  outputs
  EPS::Float64
  weightFunction
  matrix
  tree
end

function MwlsObject(vars, b, inputs, outputs, leafSize, EPS, wfun)
  MwlsObject(vars, b, inputs, outputs, EPS, wfun, b * transpose(b), KDTree(transpose(inputs)))
end

"""
    mwls(inputs::Vector{Point}, outputs::Vector{Float64}, maxDegree::Int64, leafSize::Int64, EPS::Float64, wfun::Function)

Euclidean metric is used by default.

# Arguments
- `inputs`: a 2d array where each point is on a single row,
- `outputs`: an array of output scalars,
- `maxDegree::Int64`: the maximal degree of each polynomial term in the method,
- `leafSize::Int64`: the size of the leaves in the kd-tree
- `EPS::Float64`: ε of the method,
- `wfun::Function`: wfun should be in form `wfun(distance::Float64, EPS::Float64)`.
"""
function mwls(inputs, outputs, maxDegree::Int64, leafSize::Int64, EPS::Float64, wfun::Function)
  inputDim = size(inputs, 2)
  vars::Vector{PolyVar{true}}, b::Vector{Monomial{true}} = generateMonomials(inputDim, maxDegree)
  return MwlsObject(vars, b, inputs, outputs, leafSize, EPS, wfun)
end

