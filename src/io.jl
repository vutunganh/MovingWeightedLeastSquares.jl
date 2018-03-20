"""
    `parseInputPoints(inputFilename::String)`

Parses a .csv file into a view of inputs and outputs.
"""
function parseInputPoints(inputFilename::String)
  parsedInput = readcsv(inputFilename; skipstart = 1)
  pointCount::Int = size(parsedInput, 1)
  dimensions::Int = size(parsedInput, 2)
  inputDim::Int = dimensions - 1
  inputs = view(parsedInput, :, 1:inputDim)
  outputs = view(parsedInput, :, dimensions)
  return inputs, outputs
end

