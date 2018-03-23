"""
    `readData(inputFilename::String)`

Parses a .csv file into an array of inputs and an array of outputs. Pass kwarg `skipStart = 0`, if skipping the header is not desired.
"""
function readData(inputFilename::String; skipStart = 1)
  parsedInput = readcsv(inputFilename; skipstart = skipStart)
  pointCount::Int = size(parsedInput, 1)
  dimensions::Int = size(parsedInput, 2)
  inputDim::Int = dimensions - 1
  inputs = parsedInput[:, 1:inputDim]
  outputs = parsedInput[:, dimensions]
  return inputs, outputs
end

