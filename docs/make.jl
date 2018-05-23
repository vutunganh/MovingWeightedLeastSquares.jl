using Documenter, MovingWeightedLeastSquares

makedocs(
  modules = [MovingWeightedLeastSquares],
  format = :html,
  sitename ="MovingWeightedLeastSquares.jl",
  authors = "Tung Anh Vu",
  pages = Any[
    "Home" => "index.md",
    "Constructors" => "constructors.md",
    "Approximation" => "approximation.md"
  ]
)

deploydocs(
  repo = "github.com/vutunganh/MovingWeightedLeastSquares.jl.git",
  target = "build",
  julia = "0.6",
  osname = "linux",
  deps = nothing,
  make = nothing
)
