using MovingWeightedLeastSquares
using Base.Random
using Base.Test
using Base.Iterators
using DynamicPolynomials, MultivariatePolynomials

include("sample-generator.jl")

@testset "1d approximation" begin
  xs = collect(0:0.1:10)
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
  xs = collect(0:0.1:10)
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
  xs = collect(0:0.1:10)
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

@testset "2d approximation" begin
  l = collect(-2:0.1:2);
  ins = transpose(hcat(collect.(collect(product(l, l))[:])...));
  func = (x, y) -> sin(x) + exp(-y^2)
  fs = [func(ins[i, 1], ins[i, 2]) for i in 1:size(ins, 1)]

  m = mwlsKd(ins, fs, 1, (d, eps) -> exp(-d^2))
  cll = mwlsCll(hcat(ins, fs), 1, (d, eps) -> exp(-d^2))
  testStatus = true
  for i in size(ins, 1)
    testStatus &= isapprox(m(ins[i, 1], ins[i, 2]),
                           func(ins[i, 1], ins[i, 2]),
                           atol = 0.2)
    if !testStatus
      @show ins[i, :]
      @show m(ins[i, 1], ins[i, 2])
      @show func(ins[i, 1], ins[i, 2])
      break
    end
  end
  @test testStatus == true

  for i in size(ins, 1)
    testStatus &= isapprox(cll(ins[i, 1], ins[i, 2]),
                           func(ins[i, 1], ins[i, 2]),
                           atol = 0.3)
    if !testStatus
      @show ins[i, :]
      @show m(ins[i, 1], ins[i, 2])
      @show func(ins[i, 1], ins[i, 2])
      break
    end
  end
  @test testStatus == true
end

