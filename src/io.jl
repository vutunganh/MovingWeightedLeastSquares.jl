"""
    `parseInputPoints(inputFilename::String)`

Parses a .csv file into a vector of inputs and outputs.
TODO: slice instead of copying
"""
function parseInputPoints(inputFilename::String)
  parsedInput = CSV.read(inputFilename)
  pointCount::Int = size(parsedInput, 1)
  dimensions::Int = size(parsedInput, 2)
  inputs = Vector{Point}(pointCount)
  inputDimension::Int = dimensions - 1
  outputs = Vector{Float64}(pointCount)
  for i in 1:pointCount
    p = parsedInput[i, :]
    ipt = [p[1, j] for j in 1:inputDimension]
    inputs[i] = ipt
    outputs[i] = p[1, dimensions]
  end
  return inputs, outputs
end

