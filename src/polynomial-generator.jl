"""
    generate(dim::Integer, maxdegree:Integer)

Generates an array of monomials.

# Arguments
- `dim::Int`: the amount of variables, aka dimension of the input data,
- `maxdegree::Int`: maximal degree of each term in the polynomial.

Returns the monomial variables and the monomials themselves.
"""
function generateMonomials(dim::Integer, maxdegree::Integer)
  vars = @polyvar x[1:dim]
  if dim < 0 || maxdegree < 0
    error("Cannot generate a polynomial with negative dimensions")
  end
  monomials::Vector{Monomial{true}} = [constantterm(1, sum(vars))]

  for v in vars
    tmpMonomials::Vector{Monomial{true}} = []
    for done in monomials
      @inbounds for d in 1:maxdegree
        tmp::Monomial{true} = done * v^d
        if degree(tmp) > maxdegree
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

