"""
# User provided attributes
- `inputs::Array{Float64, 2}`: a 2d array of input points where each point is in a single column,
- `outputs::Vector{Float64}`: a 2d array of of outputs where each output is in a single column,
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
  inputs::Array{Real, 2}
  outputs::Array{Real, N} where {N} # accepts both scalar and vector outputs
  EPS::Real
  weightFunc::Function
  vars::Vector{PolyVar{true}}
  b::Array{Monomial{true}}
  matrix::Array{Term{true, Int}, 2}
  tree::KDTree
end

function MwlsObject(inputs, outputs, EPS, weightFunc, vars, b; leafSize = 10)
  t = transpose(inputs)
  o = transpose(outputs)
  return MwlsObject(t, o, EPS, weightFunc,
                    vars, b,
                    b * transpose(b), KDTree(t; leafsize = leafSize))
end

"""
mwls(inputs::Array{Float64, 2}, outputs::Vector{Float64}, maxDegree::Int, leafSize::Int, EPS::Float64, weightFunc::Function)

Euclidean metric is used.

# Arguments
- `inputs`: a 2d array of input points where each point is on a single row,
- `outputs`: a vector of output scalars,
- `EPS::Float64`: ε of the method (default distance threshold for neighbor search),
- `weightFunc::Function`: weighting function of the method. It should be in form `(distance, EPS) -> Float64`.

# Keyword arguments
- `leafSize::Int64`: the size of the leaves in the kd-tree, 10 by default.
- `maxDegree::Int64`: the maximal degree of each polynomial term in the method, 2 by default.
"""
function mwls(inputs::Array{T, N},
              outputs::Array{T, M},
              EPS::Real, weightFunc::Function;
              leafSize::Int = 10, maxDegree::Int = 2) where {T <: Real, N, M}
  size(outputs, 1) != size(inputs, 1) &&
    error("The amount of inputs and outputs differs.")

  vars::Vector{PolyVar{true}}, b::Vector{Monomial{true}} =
    generateMonomials(size(inputs, 2), maxDegree)
  return MwlsObject(N == 1 ? hcat(inputs, zeros(size(inputs, 1))) : inputs,
                    #= M == 1 ? hcat(outputs, zeros(size(outputs, 1))) : outputs, =#
                    outputs,
                    EPS, weightFunc,
                    vars, b,
                    leafSize = leafSize)
end

function mwls(input::Array{T, 2},
              EPS::Real, weightFunc::Function;
              outputDim::Int = 1, leafSize::Int = 10,
              maxDegree::Int = 2) where {T <: Real}
  width = size(input, 2)
  inputEnd = width - outputDim
  outputStart = inputEnd + 1
  mwls(input[:, 1:inputEnd],
       input[:, outputStart:width],
       EPS, weightFunc,
       leafSize = leafSize, maxDegree = maxDegree)
end

