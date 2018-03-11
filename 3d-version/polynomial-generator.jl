__precompile__()
module PolynomialGenerator

using DynamicPolynomials
using MultivariatePolynomials

export gen

"""
Generates a random polynomial
@variableCount - x_1, x_2, ..., x_variableCount
@maxDegree - max degree of all terms
returns the polynomial variables and the polynomial itself

the amount of terms is given by (variableCount + maxDegree choose variableCount)
"""
function gen(variableCount::Int64, maxDegree::Int64)
  vars = @polyvar x[1:variableCount]
  if variableCount < 1
    throw("nice one")
  end
  toReturn::Vector{Monomial{true}} = [constantterm(1, sum(x))]

  for v in vars
    tmpMons::Vector{Monomial{true}} = []
    for done in toReturn
      for d in 1:maxDegree
        tmp::Monomial{true} = done * v^d
        if degree(tmp) > maxDegree
          break
        else
          push!(tmpMons, tmp)
        end
      end
    end
    for m in tmpMons
      push!(toReturn, m)
    end
  end

  return vars, toReturn
end

end

