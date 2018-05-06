using MovingWeightedLeastSquares
using Base.Random
using Base.Test
using Base.Iterators
using DynamicPolynomials, MultivariatePolynomials

include("sample-generator.jl")

@testset "cell linked list range test" begin
  range = 1000
  dataset = range * rand(2, 100) - range / 2
  cll = CellLinkedList(dataset, 3.)

  testStatus = true
  for c in 1:100
    pt = 1.1 * range * rand(2) - range * 1.1 / 2
    r = range / 3
    cllRes = cllInrange(cll, pt, r)
    naiveRes = [p for p in 1:size(dataset, 2) if norm(pt - dataset[:, p]) < r + 1e-9]
    testStatus &= sort(cllRes) == sort(naiveRes)
    if !testStatus
      @show dataset
      @show pt
      @show r
      break
    end
  end
  @test testStatus == true
end

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
  l = collect(-2:0.2:2);
  ins = transpose(hcat(collect.(collect(product(l, l))[:])...));
  func = (x, y) -> sin(x) + exp(-y^2)
  fs = [func(ins[i, 1], ins[i, 2]) for i in 1:size(ins, 1)]

  m = mwlsKd(ins, fs, 2, (d, eps) -> exp(-d^2), maxDegree = 3)
  cll = mwlsCll(hcat(ins, fs), 2, (d, eps) -> exp(-d^2), maxDegree = 3)
  testStatus = true
  for i in size(ins, 1)
    testStatus &= isapprox(m([ins[i, 1], ins[i, 2]]),
                           func(ins[i, 1], ins[i, 2]),
                           atol = 0.2)
  end
  @test testStatus == true

  for i in size(ins, 1)
    testStatus &= isapprox(cll([ins[i, 1], ins[i, 2]]),
                           func(ins[i, 1], ins[i, 2]),
                           atol = 0.3)
  end
  @test testStatus == true
end

@testset "2d derivative" begin
  l = collect(-1.5:0.2:1.5);
  ins = transpose(hcat(collect.(collect(product(l, l))[:])...));
  func = (x, y) -> exp(-(x^2 + y^2));
  dfunc = (x, y) -> 4 * x * y * exp(-(x^2 + y^2));
  fs = [func(ins[i, 1], ins[i, 2]) for i in 1:size(ins, 1)];
  dfs = [dfunc(ins[i, 1], ins[i, 2]) for i in 1:size(ins, 1)];

  m = mwlsKd(hcat(ins, fs), 1, (d, eps) -> exp(-d^2), maxDegree = 5)
  cll = mwlsCll(ins, fs, 1, (d, eps) -> exp(-d^2), maxDegree = 5)

  testStatus = true
  for i in 1:size(ins, 1)
    testStatus &= isapprox(mwlsDiff(m, [ins[i, 1], ins[i, 2]], (1, 1)),
                           dfunc(ins[i, 1], ins[i, 2]),
                           atol = 0.1)
    if !testStatus
      @show ins[i, :]
      @show mwlsDiff(m, [ins[i, 1], ins[i, 2]], (1, 1))
      @show dfunc(ins[i, 1], ins[i, 2])
      break
    end
  end
  @test testStatus == true

  testStatus = true
  for i in 1:size(ins, 1)
    testStatus &= isapprox(mwlsDiff(cll, [ins[i, 1], ins[i, 2]], (1, 1)),
                           dfunc(ins[i, 1], ins[i, 2]),
                           atol = 0.1)
    if !testStatus
      @show ins[i, :]
      @show mwlsDiff(m, [ins[i, 1], ins[i, 2]], (1, 1))
      @show dfunc(ins[i, 1], ins[i, 2])
      break
    end
  end
  @test testStatus == true
end

