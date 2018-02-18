module PolynomialGenerator

import DynamicPolynomials
import MultivariatePolynomials

"""
Generates a random polynomial
@variableCount - x_1, x_2, ..., x_variableCount
@maxDegree - max degree of all terms
returns the polynomial variables and the polynomial itself

note that the amount of variables is given by C(variableCount + maxDegree, variableCount), so don't overdo it
"""
function gen(variableCount, maxDegree::Int64)
  toReturn = []
  variables = DynamicPolynomials.@polyvar x[1:variableCount]
  for v in variables
    for d in 1:maxDegree
      push!(toReturn, v^d)
    end
  end
  for v in variables
    for done in toReturn
      if MultivariatePolynomials.nvariables(done) == 1 && MultivariatePolynomials.variables(done)[1] == v
        continue
      end
      for d in 1:maxDegree
        tmp = done * (v^d)
        if MultivariatePolynomials.degree(tmp) > maxDegree
          break
        else
          push!(toReturn, tmp)
        end
      end
    end
  end
  return variables, toReturn
end

end

