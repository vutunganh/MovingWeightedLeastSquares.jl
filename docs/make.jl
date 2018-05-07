using Documenter, MovingWeightedLeastSquares

makedocs(format = :html,
         sitename ="MovingWeightedLeastSquares.jl",
         pages = [
                  "Introduction" => "index.md"]
        )

deploydocs(repo = "github.com/vutunganh/MovingWeightedLeastSquares.jl"
           julia = "0.6")
