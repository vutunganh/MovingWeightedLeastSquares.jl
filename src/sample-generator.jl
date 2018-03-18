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

function generateRandomData(amount::Int, polyDegree::Int, polyDim::Int, ptRange::Float64, ptOffset::Float64, polyRange::Float64, polyOffset::Float64, outputFile::String)
  vars, mons = generateMonomials(polyDegree, polyDim)
  coeffs = polyRange * rand(length(mons)) - (polyRange / 2 + polyOffset)
  poly = polynomial(coeffs, mons)
  res = generatePoints(amount, ptDegree, ptRange, ptOffset)
  for p in res
    push!(p, poly(vars => p))
  end
  writeToCsv(res, outputFile)
  println("Polynomial '$pal' was used for data generation")
  println("Done")
end

