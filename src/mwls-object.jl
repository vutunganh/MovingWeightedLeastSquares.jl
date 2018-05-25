abstract type MwlsObject end

"""
# User provided member variables
- `inputs::Array{Float64, 2}`: a 2d array of input points where each point is in a single column,
- `outputs::Array{Float64, N}`: a 2d array or a vector of output points where each output is in a single column,
- `EPS::Float64`: ε of the method (cell edge length and the default distance for range search),
- `weightfunc::Function`: weighting function θ of the method.

# Automatically created member variables
If created by the `mwls_naive` function, the member variables in this section are created automatically.

- `vars::Vector{PolyVar{true}}`: variables of the polynomial,
- `b::Vector{Monomial{true}}`: the basis of the polynomial,
- `matrix`: the result of `b * transpose(b)`.
"""
struct MwlsNaiveObject <: MwlsObject
  inputs::Array{Real, 2}
  outputs::Array{Real, N} where {N} # accepts both scalar and vector outputs
  EPS::Real
  weightfunc::Function
  vars::Vector{PolyVar{true}}
  b::Vector{Monomial{true}}
  matrix::Array{Term{true, Int}, 2}
end

function MwlsNaiveObject(inputs, outputs, EPS, weightfunc, vars, b)
  t = transpose(inputs)
  o = transpose(outputs)
  return MwlsNaiveObject(t, o, EPS, weightfunc, vars, b, b * transpose(b))
end

"""
Documentation is left out on purpose to discourage the use of MwlsNaiveObject.
K-d trees are always faster and cell linked lists are faster on datasets larger than 20.
`mwls_naive` has the same signature as `mwls_cll`.
"""
function mwls_naive(inputs::Array{T, N}, outputs::Array{U, M},
                    EPS::Real, weightfunc::Function;
                    maxDegree::Integer = 2) where {T <: Real, U <: Real, N, M}
  size(outputs, 1) != size(inputs, 1) &&
    error("The amount of inputs and outputs differs")

  vars::Vector{PolyVar{true}}, b :: Vector{Monomial{true}} =
    genmonomials(size(inputs, 2), maxDegree)

  return MwlsNaiveObject(inputs, outputs, EPS, weightfunc, vars, b)
end

"""
Documentation is left out on purpose to discourage the use of MwlsNaiveObject.
K-d trees are always faster and cell linked lists are faster on datasets larger than 20.
`mwls_naive` has the same signature as `mwls_cll`.
"""
function mwls_naive(input::Array{T, 2}, EPS::Real, weightfunc::Function;
                    outputDim::Integer = 1, maxDegree::Integer = 2) where {T <: Real}
  width = size(input, 2)
  inputEnd = width - outputDim
  outputStart = inputEnd + 1
  return mwls_naive(input[:, inputEnd == 1 ? 1 : 1:inputEnd],
                    input[:, outputStart == width ? width : outputStart:width],
                    EPS, weightfunc,
                    maxDegree = maxDegree)
end

mwlsNaive(a...;b...) = warn("`mwlsNaive` is deprecated, use `mwls_naive` instead.")

"""
# User provided member variables
- `inputs::Array{Float64, 2}`: a 2d array of input points where each point is in a single column,
- `outputs::Array{Float64, N}`: a 2d array or a vector of output points where each output is in a single column,
- `EPS::Float64`: ε of the method (cell edge length and the default distance for range search),
- `weightfunc::Function`: weighting function θ of the method.

# Automatically created member variables
If created by the `mwls_cll` function, the member variables in this section are created automatically.

- `vars::Vector{PolyVar{true}}`: variables of the polynomial,
- `b::Vector{Monomial{true}}`: the basis of the polynomial,
- `matrix`: the result of `b * transpose(b)`,
- `tree`: a k-d tree for range search.
"""
struct MwlsKdObject <: MwlsObject
  inputs::Array{Real, 2}
  outputs::Array{Real, N} where {N} # accepts both scalar and vector outputs
  EPS::Real
  weightfunc::Function
  vars::Vector{PolyVar{true}}
  b::Vector{Monomial{true}}
  matrix::Array{Term{true, Int}, 2}
  tree::KDTree
end

function MwlsKdObject(inputs, outputs, EPS, weightfunc, vars, b;
                      leafSize::Integer = 10)
  t = transpose(inputs)
  o = transpose(outputs)
  return MwlsKdObject(t, o, EPS, weightfunc,
                      vars, b, b * transpose(b),
                      KDTree(t, leafsize = leafSize))
end

"""
    mwls_kd(inputs::Array{T, N}, outputs::Array{U}, EPS::Real, weightfunc::Function) where {T <: Real, U <: Real, N}
    mwls_kd(inputs::Array{T, N}, outputs::Array{U}, EPS::Real, weightfunc::Function; leafsize::Int = 10, maxDegree::Int = 2) where {T <: Real, U <: Real, N}

Creates `MwlsKdObject` from sample input and sample output data, the cutoff distance ε and a weighting function θ.

# Arguments
- `inputs`: a 2d array or a vector of input points where each point is on a single row,
- `outputs`: a 2d array or a vector of output scalars where each output is on a single row,
- `EPS::Real`: ε of the method (default distance threshold for neighbor search),
- `weightfunc::Function`: weighting function of the method. It should be in form `(distance, EPS) -> Float64`.

# Keyword arguments
- `leafSize::Int`: the size of the leaves in the k-d-tree, 10 by default.
- `maxDegree::Int`: the maximal degree of polynomials used for approximation, 2 by default.
"""
function mwls_kd(inputs::Array{T, N}, outputs::Array{U},
                 EPS::Real, weightfunc::Function;
                 leafSize::Integer = 10, maxDegree::Integer = 2) where {T <: Real, U <: Real, N}
  size(outputs, 1) != size(inputs, 1) &&
    error("The amount of inputs and outputs differs.")

  vars::Vector{PolyVar{true}}, b::Vector{Monomial{true}} =
    genmonomials(size(inputs, 2), maxDegree)
  return MwlsKdObject(N == 1 ? hcat(inputs, zeros(size(inputs, 1))) : inputs,
                      outputs, EPS, weightfunc,
                      vars, b, leafSize = leafSize)
end

"""
    mwls_kd(input::Array{T, 2}, EPS::Real, weightfunc::Function) where {T <: Real}
    mwls_kd(input::Array{T, 2}, EPS::Real, weightfunc::Function; outputDim::Int = 1, leafSize::Int = 10, maxDegree::Int = 2) where {T <: Real}

In this `mwls_kd` function, the sample input and sample output data are passed in a single array.
It is assumed that each pair of input and output is on a single row.
Dimension of the output is specified with kwarg `outputDim`.
"""
function mwls_kd(input::Array{T, 2}, EPS::Real, weightfunc::Function;
                outputDim::Integer = 1, leafSize::Integer = 10, maxDegree::Int = 2) where {T <: Real}
  width = size(input, 2)
  inputEnd = width - outputDim
  outputStart = inputEnd + 1
  return mwls_kd(input[:, inputEnd == 1 ? 1 : 1:inputEnd],
                 input[:, outputStart == width ? width : outputStart:width],
                 EPS, weightfunc,
                 leafSize = leafSize, maxDegree = maxDegree)
end

mwlsKd(a...;b...) = warn("`mwlsKd` is deprecated, use `mwls_kd` instead.")

"""
# User provided member variables
- `inputs::Array{Float64, 2}`: a 2d array of input points where each point is in a single column,
- `outputs::Array{Float64, N}`: a 2d array or a vector of output points where each output is in a single column,
- `EPS::Float64`: ε of the method (cell edge length and the default distance for range search),
- `weightfunc::Function`: weighting function θ of the method.

# Automatically created member variables
If created by the `mwls_cll` function, the member variables in this section are created automatically.

- `vars::Vector{PolyVar{true}}`: variables of the polynomial,
- `b::Vector{Monomial{true}}`: the basis of the polynomial,
- `matrix`: the result of `b * transpose(b)`,
- `cll`: a cell linked list for range search.
"""
struct MwlsCllObject <: MwlsObject
  inputs::Array{Real, 2}
  outputs::Array{Real, N} where {N} # accepts both scalar and vector outputs
  EPS::Real
  weightfunc::Function
  vars::Vector{PolyVar{true}}
  b::Vector{Monomial{true}}
  matrix::Array{Term{true, Int}, 2}
  cll::CellLinkedList
end

function MwlsCllObject(inputs, outputs, EPS, weightfunc, vars, b)
  t = transpose(inputs)
  o = transpose(outputs)
  return MwlsCllObject(t, o, EPS, weightfunc,
                       vars, b, b * transpose(b),
                       CellLinkedList(t, EPS))
end

"""
    mwls_cll(inputs::Array{T, N}, outputs::Array{U}, EPS::Real, weightfunc::Function) where {T <: Real, U <: Real, N}
    mwls_cll(inputs::Array{T, N}, outputs::Array{U}, EPS::Real, weightfunc::Function; maxDegree::Int = 2) where {T <: Real, U <: Real, N}

Creates `MwlsCllObject` from sample input and sample output data, the cutoff distance ε and a weighting function θ.

# Arguments
- `inputs`: a 2d array of input points where each point is on a single row,
- `outputs`: a 2d array or a vector of output scalars where each output is on a single row,
- `EPS::Real`: ε of the method (cell edge length and the default distance for range search),
- `weightfunc::Function`: weighting function θ of the method. It should be in form `(distance between two vectors, EPS) -> Float64`.

# Keyword arguments
- `maxDegree::Int`: the maximal degree of polynomials used for approximation, 2 by default.
"""
function mwls_cll(inputs::Array{T, N}, outputs::Array{U},
                  EPS::Real, weightfunc::Function;
                  maxDegree::Integer = 2) where {T <: Real, U <: Real, N}
  size(outputs, 1) != size(inputs, 1) &&
    error("The amount of inputs and outputs differs.")

  vars::Vector{PolyVar{true}}, b::Vector{Monomial{true}} =
    genmonomials(size(inputs, 2), maxDegree)

  return MwlsCllObject(N == 1 ? hcat(inputs, zeros(size(inputs,1))) : inputs,
                       outputs, EPS, weightfunc, vars, b)
end

"""
    mwls_cll(input::Array{T, 2}, EPS::Real, weightfunc::Function) where {T <: Real}
    mwls_cll(input::Array{T, 2}, EPS::Real, weightfunc::Function; outputDim::Int = 1, maxDegree::Int = 2) where {T <: Real}

In this `mwls_cll` function, the sample input and sample output data are passed in a single array.
It is assumed that each pair of input and output is on a single row.
Dimension of the output is specified with kwarg `outputDim`.
"""
function mwls_cll(input::Array{T, 2}, EPS::Real, weightfunc::Function;
                  outputDim::Integer = 1, maxDegree::Integer = 2) where {T <: Real}
  width = size(input, 2)
  inputEnd = width - outputDim
  outputStart = inputEnd + 1
  return mwls_cll(input[:, inputEnd == 1 ? 1 : 1:inputEnd],
                  input[:, outputStart == width ? width : outputStart:width],
                  EPS, weightfunc, maxDegree = maxDegree)
end

mwlsCll(a...;b...) = warn("`mwlsCll` is deprecated, use `mwls_cll` instead.")
