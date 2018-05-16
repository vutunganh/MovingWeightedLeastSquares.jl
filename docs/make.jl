using Documenter, MovingWeightedLeastSquares

makedocs(modules = [MovingWeightedLeastSquares],
         format = :html,
         sitename ="MovingWeightedLeastSquares.jl",
         authors = "Tung Anh Vu",
         pages = Any[
           "Home" => "index.md",
           "Constructors" => "constructors.md",
           "Approximation" => "approximation.md"
         ]
        )

deploydocs(deps = Deps.pip("mkdocs", "python-markdown-math"),
           repo = "github.com/vutunganh/MovingWeightedLeastSquares.jl",
           julia = "0.6",
           osname = "linux")
