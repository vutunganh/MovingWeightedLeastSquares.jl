using Documenter, MovingWeightedLeastSquares

makedocs()
#= makedocs(format = :html, =#
#=          sitename ="MovingWeightedLeastSquares.jl", =#
#=          pages = [ =#
#=                   "Introduction" => "index.md"] =#
#=         ) =#

deploydocs(deps = Deps.pip("mkdocs", "python-markdown-math"),
           repo = "github.com/vutunganh/MovingWeightedLeastSquares.jl",
           julia = "0.6",
           osname = "linux")
