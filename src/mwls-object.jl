abstract type MwlsObject end

"""
# User provided attributes
- `inputs::Array{Float64, 2}`: a 2d array of input points where each point is in a single column,
- `outputs::Array{Float64, N}`: a 2d array or a vector of of outputs where each output is in a single column,
- `EPS::Float64`: ε of the method (default distance threshold for neighbor search),
- `weightFunc::Function`: weighting function of the method.

# Automatically created attributes
If created by the `mwls` function, the attributes in this section are created automatically.

- `vars::Vector{PolyVar{true}}`: variables of the polynomial,
- `b::Vector{Monomial{true}}`: the basis of the polynomial,
- `matrix`: the result of `b * transpose(b)`,
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

function mwlsNaive(inputs::Array{T, N}, outputs::Array{U, M},
                   EPS::Real, weightFunc::Function;
                   maxDegree::Int = 2) where {T <: Real, U <: Real, N, M}
  size(outputs, 1) != size(inputs, 1) &&
    error("The amount of inputs and outputs differs")

  vars::Vector{PolyVar{true}}, b :: Vector{Monomial{true}} =
    generateMonomials(size(inputs, 2), maxDegree)

  return MwlsNaiveObject(inputs, outputs, EPS, weightFunc, vars, b)
end

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
# User provided attributes
- `inputs::Array{Float64, 2}`: a 2d array of input points where each point is in a single column,
- `outputs::Array{Float64, N}`: a 2d array or a vector of of outputs where each output is in a single column,
- `EPS::Float64`: ε of the method (default distance threshold for neighbor search),
- `weightFunc::Function`: weighting function of the method.

# Automatically created attributes
If created by the `mwls` function, the attributes in this section are created automatically.

- `vars::Vector{PolyVar{true}}`: variables of the polynomial,
- `b::Vector{Monomial{true}}`: the basis of the polynomial,
- `matrix`: the result of `b * transpose(b)`,
- `tree`: a k-d tree for neareast neighbor search.
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
Creates MwlsKdObject from inputs, outputs, ϵ and a weighting function. Euclidean metric is used.

# Arguments
- `inputs`: a 2d array of input points where each point is on a single row,
- `outputs`: a 2d array or a vector of output scalars where each output is on a single row,
- `EPS::Float64`: ε of the method (default distance threshold for neighbor search),
- `weightFunc::Function`: weighting function of the method. It should be in form `(distance, EPS) -> Float64`.

# Keyword arguments
- `leafSize::Int64`: the size of the leaves in the kd-tree, 10 by default.
- `maxDegree::Int64`: the maximal degree of each polynomial term in the method, 2 by default.
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
In this mwlsKd function, the inputs and outupts are passed in a single array.
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
# User provided attributes
- `inputs::Array{Float64, 2}`: a 2d array of input points where each point is in a single column,
- `outputs::Array{Float64, N}`: a 2d array or a vector of of outputs where each output is in a single column,
- `EPS::Float64`: ε of the method (default distance threshold for neighbor search),
- `weightFunc::Function`: weighting function of the method.

# Automatically created attributes
If created by the `mwls` function, the attributes in this section are created automatically.

- `vars::Vector{PolyVar{true}}`: variables of the polynomial,
- `b::Vector{Monomial{true}}`: the basis of the polynomial,
- `matrix`: the result of `b * transpose(b)`,
- `cll`: a cell linked list for nearest neighbor search.
"""
mutable struct MwlsCllObject <: MwlsObject
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
Creates MwlsCllObject from inputs, outputs, ϵ and a weighting function. Euclidean metric is used.

# Arguments
- `inputs`: a 2d array of input points where each point is on a single row,
- `outputs`: a 2d array or a vector of output scalars where each output is on a single row,
- `EPS::Float64`: ε of the method (grid edge length and default distance threshold for neighbor search),
- `weightFunc::Function`: weighting function of the method. It should be in form `(distance, EPS) -> Float64`.

# Keyword arguments
- `maxDegree::Int64`: the maximal degree of each polynomial term in the method, 2 by default.
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
In this mwlsCll function, the inputs and outupts are passed in a single array.
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

