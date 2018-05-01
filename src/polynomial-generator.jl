"""
    `generate(variableCount::Int, maxDegree:Int)`

Generates an array of monomials.

# Arguments
- `variableCount::Int`: the amount of variables, aka dimension of variable,
- `maxDegree::Int`: maximal degree of each term in the polynomial, aka spatial dimension.

Returns the monomial variables and the monomials themselves.

In general, the amount of terms is given by `{variableCount + maxDegree \choose variableCount}`.
"""
function generateMonomials(variableCount::Int, maxDegree::Int)
  vars = @polyvar x[1:variableCount]
  if variableCount < 0 || maxDegree < 0
    error("Cannot generate a polynomial with negative dimensions")
  end
  monomials::Vector{Monomial{true}} = [constantterm(1, sum(vars))]

  for v in vars
    tmpMonomials::Vector{Monomial{true}} = []
    for done in monomials
      @inbounds for d in 1:maxDegree
        tmp::Monomial{true} = done * v^d
        if degree(tmp) > maxDegree
          break
        else
          push!(tmpMonomials, tmp)
        end
      end
    end
    for m in tmpMonomials
      push!(monomials, m)
    end
  end

  return vars, monomials
end

