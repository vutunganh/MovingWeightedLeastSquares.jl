import DynamicPolynomials
import MultivariatePolynomials

function gen(variables, maxDegree::Int64)
  toReturn = []
  for v in variables
    for d in 1:maxDegree
      push!(toReturn, v^d)
    end
  end
  for v in variables
    println("done: $toReturn")
    for done in toReturn
      if MultivariatePolynomials.nvariables(done) == 1 && MultivariatePolynomials.variables(done)[1] == v
        continue
      end
      println("doing $done, var $v")
      for d in 1:maxDegree
        tmp = done * (v^d)
        if MultivariatePolynomials.degree(tmp) > maxDegree
          break
        else
          println("pushing $tmp")
          push!(toReturn, tmp)
        end
      end
    end
  end
  return toReturn
end

DynamicPolynomials.@polyvar x[1:3]
println(gen(x[1:3], 4))

