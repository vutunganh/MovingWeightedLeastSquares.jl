using MovingWeightedLeastSquares
using Base.Random
using Base.Test
using DynamicPolynomials, MultivariatePolynomials

include("sample-generator.jl")

function randRange(mean::Real, width::Real)
  return 
end

@testset "1d approximation" begin
  xs = [j * 0.1 for j in 0:100]
  fs = [sin(x) for x in xs]

  m = mwlsKd(xs, fs, 1, (d, e) -> exp(-d^2))
  cll = mwlsCll(xs, fs, 1, (d, e) -> exp(-d^2))

  testStatus = true
  for p in 0:0.1:10
    testStatus &= isapprox(m(p), sin(p), atol=0.1)
  end
  @test testStatus == true

  testStatus = true
  for p in 0:0.1:10
    testStatus &= isapprox(cll(p), sin(p), atol=0.1)
  end
  @test testStatus == true

  # not enough data for approximation, zero returned by default
  @test isapprox(m(-0.9), 0, atol=0.1)
  @test isapprox(cll(-0.9), 0, atol=0.1)
  @test isapprox(m(10.9), 0, atol=0.1)
  @test isapprox(cll(10.9), 0, atol=0.1)
end

@testset "single table data" begin
  xs = [j * 0.1 for j in 0:100]
  fs = [sin(x) for x in xs]
  data = hcat(xs, fs)

  m = mwlsKd(data, 1, (d, e) -> exp(-d^2))
  cll = mwlsCll(data, 1, (d, e) -> exp(-d^2))

  testStatus = true
  for p in 0:0.1:10
    testStatus &= isapprox(m(p), sin(p), atol=0.1)
  end
  @test testStatus == true

  testStatus = true
  for p in 0:0.1:10
    testStatus &= isapprox(cll(p), sin(p), atol=0.1)
  end
  @test testStatus == true

  # not enough data for approximation, zero returned by default
  @test isapprox(m(-0.9), 0, atol=0.1)
  @test isapprox(cll(-0.9), 0, atol=0.1)
  @test isapprox(m(10.9), 0, atol=0.1)
  @test isapprox(cll(10.9), 0, atol=0.1)
end

@testset "1d derivative" begin
  xs = [j * 0.1 for j in 0:100]
  fs = [sin(x) for x in xs]

  m = mwlsKd(xs, fs, 1, (d, e) -> exp(-d^2))
  cll = mwlsCll(xs, fs, 1, (d, e) -> exp(-d^2))

  testStatus = true
  for p in 0:0.1:10
    isapprox(mwlsDiff(m, p, 1), cos(p), atol=0.1)
  end
  @test testStatus == true

  testStatus = true
  for p in 0:0.1:10
    testStatus &= isapprox(mwlsDiff(cll, p, 1), cos(p), rtol=0.1)
  end
  @test testStatus == true

  # not enough data for approximation, zero returned by default
  @test isapprox(mwlsDiff(m, -0.9, 1), 0, atol=0.1)
  @test isapprox(mwlsDiff(cll, -0.9, 1), 0, atol=0.1)
  @test isapprox(mwlsDiff(m, 10.9, 1), 0, atol=0.1)
  @test isapprox(mwlsDiff(cll, 10.9, 1), 0, atol=0.1)
end
