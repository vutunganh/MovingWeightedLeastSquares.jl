# Weighted least squares in 3D space.
# Approximation will be done with polynomials.
# Assume that approximated function will be in form f: R^2 -> R
# The primary purpose of this program is to find all mistakes, that will
# happen in the real program and solve them early.

include("polynomial-generator.jl")

Base.__precompile__()
module MWLS3D
using PolynomialGenerator
using DataFrames
using CSV
using DynamicPolynomials
using MultivariatePolynomials

export wls, parseInputPoints, WlsObject, Point

const Point = Vector{Float64}

"""
reads an input csv file

returns input "points" and their respective outputs
"""
function parseInputPoints(inputFilename::String)
  parsedInput = CSV.read(inputFilename)
  pointCount::Int = size(parsedInput, 1)
  dimensions::Int = size(parsedInput, 2)
  inputs::Vector{Point} = []
  inputDimension::Int = dimensions - 1
  outputs::Vector{Float64} = []
  for i in 1:pointCount
    p = parsedInput[i, :]
    ipt = [p[1, j] for j in 1:inputDimension]
    push!(inputs, ipt)
    push!(outputs, p[1, dimensions])
  end
  return inputs, outputs
end

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
    `(wo::WlsObject)(inputPoint::Point, value::Bool)`

Approximates the `WlsObject` at `inputPoint`.
If `value` is true, the approximated value will be outputted instead of the coefficients.
"""
function (obj::WlsObject)(inPt::Point, value::Bool)
  m = length(obj.b)
  datalen = length(obj.inputs)
  firstTerm = zeros(m, m)
  secondTerm = zeros(m)

  for p in 1:datalen
    curPt = obj.inputs[p]
    w = obj.weightFunction(norm(inPt - curPt), obj.EPS)
    for i in 1:m
      for j in 1:m
        firstTerm[i, j] += w * obj.matrix[i, j](obj.vars => curPt)
      end
      secondTerm[i] += w * obj.b[i](obj.vars => curPt) * obj.outputs[p]
    end
  end

  coefficients = firstTerm \ secondTerm
  if value
    poly = sum(coefficients .* obj.b)
    return poly(obj.vars => inPt)
  else
    return coefficients
  end
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
  vars::Vector{PolyVar{true}}, b::Vector{Monomial{true}} = gen(inputDimension, maxDegree)
  return WlsObject(vars, b, inputs, outputs, EPS, wfun)
end
end

