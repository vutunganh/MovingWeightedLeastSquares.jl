abstract type MwlsObject end

"""
# User provided member variables
- `inputs::Array{Float64, 2}`: a 2d array of input points where each point is in a single column,
- `outputs::Array{Float64, N}`: a 2d array or a vector of output points where each output is in a single column,
- `EPS::Float64`: ε of the method (cell edge length and the default distance for range search),
- `weightFunc::Function`: weighting function θ of the method.

# Automatically created member variables
If created by the `mwlsNaive` function, the member variables in this section are created automatically.

- `vars::Vector{PolyVar{true}}`: variables of the polynomial,
- `b::Vector{Monomial{true}}`: the basis of the polynomial,
- `matrix`: the result of `b * transpose(b)`.
"""
struct MwlsNaiveObject <: MwlsObject
  inputs::Array{Real, 2}
  outputs::Array{Real, N} where {N} # accepts both scalar and vector outputs
  EPS::Real
  weightFunc::Function
  vars::Vector{PolyVar{true}}
  b::Vector{Monomial{true}}
  matrix::Array{Term{true, Int}, 2}
end

function MwlsNaiveObject(inputs, outputs, EPS, weightFunc, vars, b)
  t = transpose(inputs)
  o = transpose(outputs)
  return MwlsNaiveObject(t, o, EPS, weightFunc, vars, b, b * transpose(b))
end

"""
Documentation is left out on purpose to discourage the use of MwlsNaiveObject.
K-d trees are always faster and cell linked lists are faster on datasets larger than 20.
`mwlsNaive` has the same signature as `mwlsCll`.
"""
function mwlsNaive(inputs::Array{T, N}, outputs::Array{U, M},
                   EPS::Real, weightFunc::Function;
                   maxDegree::Int = 2) where {T <: Real, U <: Real, N, M}
  size(outputs, 1) != size(inputs, 1) &&
    error("The amount of inputs and outputs differs")

  vars::Vector{PolyVar{true}}, b :: Vector{Monomial{true}} =
    generateMonomials(size(inputs, 2), maxDegree)

  return MwlsNaiveObject(inputs, outputs, EPS, weightFunc, vars, b)
end

"""
Documentation is left out on purpose to discourage the use of MwlsNaiveObject.
K-d trees are always faster and cell linked lists are faster on datasets larger than 20.
`mwlsNaive` has the same signature as `mwlsCll`.
"""
function mwlsNaive(input::Array{T, 2}, EPS::Real, weightFunc::Function;
                   outputDim::Int = 1, maxDegree::Int = 2) where {T <: Real}
  width = size(input, 2)
  inputEnd = width - outputDim
  outputStart = inputEnd + 1
  return mwlsNaive(input[:, inputEnd == 1 ? 1 : 1:inputEnd],
                   input[:, outputStart == width ? width : outputStart:width],
                   EPS, weightFunc,
                   maxDegree = maxDegree)
end

"""
# User provided member variables
- `inputs::Array{Float64, 2}`: a 2d array of input points where each point is in a single column,
- `outputs::Array{Float64, N}`: a 2d array or a vector of output points where each output is in a single column,
- `EPS::Float64`: ε of the method (cell edge length and the default distance for range search),
- `weightFunc::Function`: weighting function θ of the method.

# Automatically created member variables
If created by the `mwlsCll` function, the member variables in this section are created automatically.

- `vars::Vector{PolyVar{true}}`: variables of the polynomial,
- `b::Vector{Monomial{true}}`: the basis of the polynomial,
- `matrix`: the result of `b * transpose(b)`,
- `tree`: a k-d tree for range search.
"""
struct MwlsKdObject <: MwlsObject
  inputs::Array{Real, 2}
  outputs::Array{Real, N} where {N} # accepts both scalar and vector outputs
  EPS::Real
  weightFunc::Function
  vars::Vector{PolyVar{true}}
  b::Vector{Monomial{true}}
  matrix::Array{Term{true, Int}, 2}
  tree::KDTree
end

function MwlsKdObject(inputs, outputs, EPS, weightFunc, vars, b;
                      leafSize::Int = 10)
  t = transpose(inputs)
  o = transpose(outputs)
  return MwlsKdObject(t, o, EPS, weightFunc,
                      vars, b, b * transpose(b),
                      KDTree(t, leafsize = leafSize))
end

"""
`mwlsKd(inputs::Array{T, N}, outputs::Array{U}, EPS::Real, weightFunc::Function) where {T <: Real, U <: Real, N}`
`mwlsKd(inputs::Array{T, N}, outputs::Array{U}, EPS::Real, weightFunc::Function; leafsize::Int = 10, maxDegree::Int = 2) where {T <: Real, U <: Real, N}`

Creates `MwlsKdObject` from sample input and sample output data, the cutoff distance ε and a weighting function θ.

# Arguments
- `inputs`: a 2d array of input points where each point is on a single row,
- `outputs`: a 2d array or a vector of output scalars where each output is on a single row,
- `EPS::Real`: ε of the method (default distance threshold for neighbor search),
- `weightFunc::Function`: weighting function of the method. It should be in form `(distance, EPS) -> Float64`.

# Keyword arguments
- `leafSize::Int`: the size of the leaves in the k-d-tree, 10 by default.
- `maxDegree::Int`: the maximal degree of polynomials used for approximation, 2 by default.
"""
function mwlsKd(inputs::Array{T, N}, outputs::Array{U},
                EPS::Real, weightFunc::Function;
                leafSize::Int = 10, maxDegree::Int = 2) where {T <: Real, U <: Real, N}
  size(outputs, 1) != size(inputs, 1) &&
    error("The amount of inputs and outputs differs.")

  vars::Vector{PolyVar{true}}, b::Vector{Monomial{true}} =
    generateMonomials(size(inputs, 2), maxDegree)
  return MwlsKdObject(N == 1 ? hcat(inputs, zeros(size(inputs, 1))) : inputs,
                      outputs, EPS, weightFunc,
                      vars, b, leafSize = leafSize)
end

"""
  `mwlsKd(input::Array{T, 2}, EPS::Real, weightFunc::Function) where {T <: Real}`
  `mwlsKd(input::Array{T, 2}, EPS::Real, weightFunc::Function; outputDim::Int = 1, leafSize::Int = 10, maxDegree::Int = 2) where {T <: Real}`

In this `mwlsKd` function, the sample input and sample output data are passed in a single array.
It is assumed that each pair of input and output is on a single row.
Dimension of the output is specified with kwarg `outputDim`.
"""
function mwlsKd(input::Array{T, 2}, EPS::Real, weightFunc::Function;
                outputDim::Int = 1, leafSize::Int = 10, maxDegree::Int = 2) where {T <: Real}
  width = size(input, 2)
  inputEnd = width - outputDim
  outputStart = inputEnd + 1
  return mwlsKd(input[:, inputEnd == 1 ? 1 : 1:inputEnd],
                input[:, outputStart == width ? width : outputStart:width],
                EPS, weightFunc,
                leafSize = leafSize, maxDegree = maxDegree)
end

"""
# User provided member variables
- `inputs::Array{Float64, 2}`: a 2d array of input points where each point is in a single column,
- `outputs::Array{Float64, N}`: a 2d array or a vector of output points where each output is in a single column,
- `EPS::Float64`: ε of the method (cell edge length and the default distance for range search),
- `weightFunc::Function`: weighting function θ of the method.

# Automatically created member variables
If created by the `mwlsCll` function, the member variables in this section are created automatically.

- `vars::Vector{PolyVar{true}}`: variables of the polynomial,
- `b::Vector{Monomial{true}}`: the basis of the polynomial,
- `matrix`: the result of `b * transpose(b)`,
- `cll`: a cell linked list for range search.
"""
struct MwlsCllObject <: MwlsObject
  inputs::Array{Real, 2}
  outputs::Array{Real, N} where {N} # accepts both scalar and vector outputs
  EPS::Real
  weightFunc::Function
  vars::Vector{PolyVar{true}}
  b::Vector{Monomial{true}}
  matrix::Array{Term{true, Int}, 2}
  cll::CellLinkedList
end

function MwlsCllObject(inputs, outputs, EPS, weightFunc, vars, b)
  t = transpose(inputs)
  o = transpose(outputs)
  return MwlsCllObject(t, o, EPS, weightFunc,
                       vars, b, b * transpose(b),
                       CellLinkedList(t, EPS))
end

"""
    `mwlsCll(inputs::Array{T, N}, outputs::Array{U}, EPS::Real, weightFunc::Function) where {T <: Real, U <: Real, N}`
    `mwlsCll(inputs::Array{T, N}, outputs::Array{U}, EPS::Real, weightFunc::Function; maxDegree::Int = 2) where {T <: Real, U <: Real, N}`

Creates `MwlsCllObject` from sample input and sample output data, the cutoff distance ε and a weighting function θ.

# Arguments
- `inputs`: a 2d array of input points where each point is on a single row,
- `outputs`: a 2d array or a vector of output scalars where each output is on a single row,
- `EPS::Real`: ε of the method (cell edge length and the default distance for range search),
- `weightFunc::Function`: weighting function θ of the method. It should be in form `(distance between two vectors, EPS) -> Float64`.

# Keyword arguments
- `maxDegree::Int`: the maximal degree of polynomials used for approximation, 2 by default.
"""
function mwlsCll(inputs::Array{T, N}, outputs::Array{U},
                 EPS::Real, weightFunc::Function;
                 maxDegree::Int = 2) where {T <: Real, U <: Real, N}
  size(outputs, 1) != size(inputs, 1) &&
    error("The amount of inputs and outputs differs.")

  vars::Vector{PolyVar{true}}, b::Vector{Monomial{true}} =
    generateMonomials(size(inputs, 2), maxDegree)

  return MwlsCllObject(N == 1 ? hcat(inputs, zeros(size(inputs,1))) : inputs,
                       outputs, EPS, weightFunc, vars, b)
end

"""
    `mwlsCll(input::Array{T, 2}, EPS::Real, weightFunc::Function) where {T <: Real}`
    `mwlsCll(input::Array{T, 2}, EPS::Real, weightFunc::Function; outputDim::Int = 1, maxDegree::Int = 2) where {T <: Real}`

In this `mwlsCll` function, the sample input and sample output data are passed in a single array.
It is assumed that each pair of input and output is on a single row.
Dimension of the output is specified with kwarg `outputDim`.
"""
function mwlsCll(input::Array{T, 2}, EPS::Real, weightFunc::Function;
                 outputDim::Int = 1, maxDegree::Int = 2) where {T <: Real}
  width = size(input, 2)
  inputEnd = width - outputDim
  outputStart = inputEnd + 1
  return mwlsCll(input[:, inputEnd == 1 ? 1 : 1:inputEnd],
                 input[:, outputStart == width ? width : outputStart:width],
                 EPS, weightFunc, maxDegree = maxDegree)
end
