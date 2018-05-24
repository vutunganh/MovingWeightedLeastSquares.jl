@testset "1d approximation" begin
  xs = collect(0:0.1:10)
  fs = [sin(x) for x in xs]

  m = mwls_kd(xs, fs, 1, (d, e) -> exp(-d^2))
  cll = mwls_cll(xs, fs, 1, (d, e) -> exp(-d^2))

  teststatus = true
  for p in 0:0.1:10
    teststatus &= isapprox(m(p), sin(p), atol=0.1)
  end
  @test teststatus == true

  teststatus = true
  for p in 0:0.1:10
    teststatus &= isapprox(cll(p), sin(p), atol=0.1)
  end
  @test teststatus == true

  # not enough data for approximation, zero returned by default
  @test isapprox(m(-0.9), 0, atol=0.1)
  @test isapprox(cll(-0.9), 0, atol=0.1)
  @test isapprox(m(10.9), 0, atol=0.1)
  @test isapprox(cll(10.9), 0, atol=0.1)
end

