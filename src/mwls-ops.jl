# return the indices of sample inputs in `obj::MwlsNaiveObject` within `dist` from `pt`
function inrange_data(obj::MwlsNaiveObject, pt::Point, dist::Real = obj.EPS)
  res::Vector{Int} = []

  @inbounds for p in 1:size(obj.inputs, 2)
    if norm(obj.inputs[:, p] - pt) < dist + dist / 1e8
      push!(res, p)
    end
  end

  return res
end

# return the indices of sample inputs in `obj::MwlsNaiveObject` within `dist` from `pt`
function inrange_data(obj::MwlsKdObject, pt::Point, dist::Real = obj.EPS)
  return inrange(obj.tree, pt, dist)
end

# return the indices of sample inputs in `obj::MwlsNaiveObject` within `dist` from `pt`
function inrange_data(obj::MwlsCllObject, pt::Point, dist::Real = obj.EPS)
  return cll_inrange(obj.cll, pt, dist)
end

getInrangeData(a...;b...) = warn("`getInrangeData` is deprecated, use` inrange_data` instead.")

"""
    mwls_coefficients(obj::MwlsObject, pt::Point, dist::Real)

Calculates the coefficients of the linear combination of polynomials used for approximation.
This is done for each dimension of output data.

!!! note
    If the matrix in the system of linear equations used to find the coefficients is singular, then zero coefficients are returned!
"""
function mwls_coefficients(obj::MwlsObject, pt::Point, dist::Real = obj.EPS)
  m = size(obj.b, 1)
  outputDim = size(obj.outputs, 1)
  firstTerm = zeros(m, m)
  secondTerm = zeros(outputDim, m)

  data = inrange_data(obj, pt, dist)

  for p in data
    curPt = obj.inputs[:, p]
    w = obj.weightFunc(norm(pt - curPt), obj.EPS)
    @inbounds for j in 1:m
      @inbounds for i in 1:m
        firstTerm[i, j] += w * obj.matrix[i, j](obj.vars => curPt)
      end
      secondTerm[:, j] += w * obj.b[j](obj.vars => curPt) * obj.outputs[:, p]
    end
  end

  # probably singular matrix
  # returns zeroes
  if cond(firstTerm, 1) >= 1e14 || cond(firstTerm, 2) >= 1e14
    return zeros(m, outputDim)
  end

  result = firstTerm \ transpose(secondTerm)
  return result
end

calcMwlsCoefficients(a...;b...) = warn("`calcMwlsCoefficients` is deprecated, use `mwls_coefficients` instead.")

# does the actual approximation
"""
    approximate(obj::MwlsObject, pt::Point)
    approximate(obj::MwlsObject, pt::Point; dist::Real = obj.EPS)

This calculates the approximated value at `pt` for each dimension of output data.
The approximated value is returned.
"""
function approximate(obj::MwlsObject, pt::Point, dist::Real = obj.EPS)
  cs = mwls_coefficients(obj, pt, dist)
  poly = [polynomial(cs[:, i], obj.b) for i in 1:size(cs, 2)]
  res = [p(obj.vars => pt) for p in poly]
  return length(res) == 1 ? res[1] : res
end

function (obj::MwlsNaiveObject)(pt::Point, dist::Real = obj.EPS)
  approximate(obj, pt, dist)
end

function (obj::MwlsNaiveObject)(pt::Real, dist::Real = obj.EPS)
  return obj([pt], dist)
end

function (obj::MwlsKdObject)(pt::Point, dist::Real = obj.EPS)
  approximate(obj, pt, dist)
end

function (obj::MwlsKdObject)(pt::Real, dist::Real = obj.EPS)
  return obj([pt, 0], dist)
end

function (obj::MwlsCllObject)(pt::Point, dist::Real = obj.EPS)
  approximate(obj, pt, dist)
end

# TODO: huh
function (obj::MwlsCllObject)(pt::Real, dist::Real = obj.EPS)
  return obj([pt, 0], dist)
end

"""
    mwls_diff_polys(obj::MwlsObject, inPt::Point, dirs::NTuple{N, Integer}; dist::Real = obj.EPS) where {N}

Polynomials created by [mwls_coefficients](@ref) are differentiated according to dirs.
For an example if ``f`` is the polynomial used for approximation, then ``dirs = (1,2)`` returns ``\\frac{\\partial^3}{\\partial x_1 \\partial x_2^2}f``.
"""
function mwls_diff_polys(obj::MwlsObject, inPt::Point, dirs::NTuple{N, Integer}; dist::Real = obj.EPS) where {N}
  cs = mwls_coefficients(obj, inPt, dist)
  poly = [polynomial(cs[:, i], obj.b) for i in 1:size(cs, 2)]
  @inbounds for j in 1:length(poly)
    @inbounds for i in 1:length(dirs)
      poly[j] = differentiate(poly[j], obj.vars[i], dirs[i])
    end
  end
  return poly
end

calcDiffMwlsPolys(a...;b...) = warn("`calcDiffMwlsPolys` is deprecated, use `mwls_diff_polys` instead")

"""
    mwls_diff(obj::MwlsNaiveObject, inPt::Real, dirs::Integer; dist::Real = obj.EPS)
"""
function mwls_diff(obj::MwlsNaiveObject, inPt::Real, dirs::Integer; dist::Real = obj.EPS)
  mwls_diff(obj, [inPt], ntuple(x -> dirs, 1); dist = dist)
end

"""
    mwls_diff(obj::MwlsNaiveObject, inPt::Point, dirs::Integer; dist::Real = obj.EPS)
"""
function mwls_diff(obj::MwlsNaiveObject, inPt::Point, dirs::Integer; dist::Real = obj.EPS)
  mwls_diff(obj, inPt, ntuple(x -> dirs, 1); dist = dist)
end

"""
    mwls_diff(obj::MwlsNaiveObject, inPt::Real, dirs::NTuple{N, Integer}; dist::Real = obj.EPS) where {N}
"""
function mwls_diff(obj::MwlsNaiveObject, inPt::Real, dirs::NTuple{N, Integer}; dist::Real = obj.EPS) where {N}
  mwls_diff(obj, [inPt], dirs, dist)
end

"""
    mwls_diff(obj::MwlsObject, inPt::Real, dirs::Integer; dist::Real = obj.EPS)
"""
function mwls_diff(obj::MwlsObject, inPt::Real, dirs::Integer; dist::Real = obj.EPS)
  mwls_diff(obj, [inPt, 0], ntuple(x -> dirs, 1); dist = dist)
end

"""
    mwls_diff(obj::MwlsObject, inPt::Real, dirs::NTuple{N, Integer}; dist::Real = obj.EPS) where {N}
"""
function mwls_diff(obj::MwlsObject, inPt::Real, dirs::NTuple{N, Integer}; dist::Real = obj.EPS) where {N}
  mwls_diff(obj, [inPt, 0], ntuple(x -> dirs, 1); dist = dist)
end

"""
    mwls_diff(obj::MwlsObject, inPt::Point, dirs::Integer; dist::Real = obj.EPS)
"""
function mwls_diff(obj::MwlsObject, inPt::Point, dirs::Integer; dist::Real = obj.EPS)
  mwls_diff(obj, inPt, ntuple(x -> dirs, 1); dist = dist)
end

"""
    mwls_diff(obj::MwlsObject, inPt::Point, dirs::NTuple{N, Integer}; dist::Real = obj.EPS)

Calculates the approximated derivative at `inPt`, where `x[i]` is differentiated `dirs[i]` times.
"""
function mwls_diff(obj::MwlsObject, inPt::Point, dirs::NTuple{N, Integer}; dist::Real = obj.EPS) where {N}
  N != length(obj.vars) && error("Mismatch between tuple size and amount of variables")
  pl = mwls_diff_polys(obj, inPt, dirs; dist = dist)
  res = [p(obj.vars => inPt) for p in pl]
  return length(res) == 1 ? res[1] : res
end

mwlsDiff(a...;b...) = warn("`mwlsDiff` is deprecated, use `mwls_diff` instead")
