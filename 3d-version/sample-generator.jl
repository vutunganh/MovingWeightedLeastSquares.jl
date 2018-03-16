include("polynomial-generator.jl")

using PolynomialGenerator
using DynamicPolynomials
using MultivariatePolynomials
using Base.Random
using ArgParse
using DataFrames
using CSV

export main

function randUniform(a, b)
  return a + rand() * (b - a)
end

"""
    `generatePoint(dims::Int, range::AbstractFloat, offset::AbstractFloat = 0)`

Generates a random `dims`-dimensional point represented as an array.

# Arguments
- `range`: specifies symmetric bounds for each dimension of the generated point,
- `offset`: offsets the center of `range` interval.
"""
function generatePoint(dims::Int, range::AbstractFloat, offset::AbstractFloat = 0.0)
  return range * rand(dims) - range / 2 + offset 
end

"Generates `n` random points without outputs. The rest of parameters are same as arguments for `generatePoint`."
function generatePoints(n::Int, dims::Int, range::AbstractFloat, offset::AbstractFloat = 0.0)
  return [generatePoint(dims, range, offset) for i in 1:n]
end

function writeToCsv(data::Vector{Vector{Float64}}, output::String)
  for i in 1:(length(data) - 1)
    if length(data[i]) != length(data[i + 1])
      error("Array dimension mismatch on row $i")
    end
  end

  df = DataFrame()
  for i in 1:(length(data[1]) - 1)
    df[Symbol("x", i)] = [d[i] for d in data]
  end
  df[Symbol("y")] = [d[length(data[1])] for d in data]
  CSV.write(output, df)
end

function parseCommandline()
  s = ArgParseSettings()
  @add_arg_table s begin
    "-n", "--amount"
      help = "amount of generated points"
      arg_type = Int
      required = true
    "-d", "--degree"
      help = "the degree (amount of variables) of the polynomial used for generation"
      arg_type = Int
      required = true
    "-s", "--spatial-dimensions"
      help = "the degree of each term in the polynomial used for generation"
      arg_type = Int
      required = true
    "-r", "--range"
      help = "a symmetric range of each dimension of inputs"
      arg_type = Float64
      required = true
    "-f", "--offset"
      help = "an offset of --range, 0 by default"
      arg_type = Float64
      default = 0.0
    "-p", "--polynomial-range"
      help = "a symmetric range of each coefficient of polynomial used for data generation"
      arg_type = Float64
      required = true
    "-t", "--polynomial-offset"
      help = "an offset of --polynomial-range, 0 by default"
      arg_type = Float64
      default = 0.0
    "-o", "--output"
      help = "output file for generated points"
      arg_type = String
      required = true
  end
  return parse_args(s)
end

function main()  
  options = parseCommandline()

  vars, mons = PolynomialGenerator.generate(options["degree"], options["spatial-dimensions"])
  coeffs =
  pal = polynomial(options["polynomial-range"] * rand(length(mons)) - options["polynomial-range"] / 2 + options["polynomial-offset"], mons)
  randomPoints = generatePoints(options["amount"], options["degree"], float(options["range"]), float(options["offset"]))
  for p in randomPoints
    push!(p, pal(vars => p))
  end

  writeToCsv(randomPoints, options["output"])
  println("Polynomial '$pal' was used for data generation")
  println("Done")
end

