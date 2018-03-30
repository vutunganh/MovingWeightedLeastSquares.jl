"""
    `generatePoint(dims::Int, range::AbstractFloat, offset::AbstractFloat = 0)`

Generates a random `dims`-dimensional point represented as an array.

d ∈ (offset - range / 2, offset + range / 2) for each dimension of the point.
"""
function generatePoint(dims::Int, range::AbstractFloat, offset::AbstractFloat = 0.0)
  return range * rand(dims) - range / 2 + offset 
end

"""
Generates `n` random points without outputs. The rest of parameters are same as parameters for `generatePoint`.
"""
function generatePoints(n::Int, dims::Int, range::AbstractFloat, offset::AbstractFloat = 0.0)
  return [generatePoint(dims, range, offset) for i in 1:n]
end

"""
# Arguments
- `amount`: the amount of points to be generated,
- `polyDegree`: the degree of polynomial used for data generation,
- `polyDim`: the spatial dimension of each term in the polynomial used for data generation,
- `ptRange` and `ptOffset`: the same as `generatePoint`,
- `polyRange` and `polyOffset`: the same as `generatePoint`, but for each term in the polynomial used for data generation,
- `outputFile`: name of the output file
"""
function generateRandomData(amount::Int, polyDegree::Int, polyDim::Int, ptRange::Float64, ptOffset::Float64, polyRange::Float64, polyOffset::Float64, outputFile::String)
  vars, mons = generateMonomials(polyDegree, polyDim)
  coeffs = polyRange * rand(length(mons)) - (polyRange / 2 + polyOffset)
  poly = polynomial(coeffs, mons)
  res = generatePoints(amount, polyDegree, ptRange, ptOffset)
  for p in res
    push!(p, poly(vars => p))
  end
  writecsv(outputFile, res)
  println("Polynomial '$poly' was used for data generation")
  println("Done")
end

