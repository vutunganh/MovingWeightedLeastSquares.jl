@testset "1d derivative" begin
  xs = collect(0:0.1:10)
  fs = [sin(x) for x in xs]

  m = mwls_kd(xs, fs, 1, (d, e) -> exp(-d^2))
  cll = mwls_cll(xs, fs, 1, (d, e) -> exp(-d^2))

  teststatus = true
  for p in 0:0.1:10
    isapprox(mwls_diff(m, p, 1), cos(p), atol=0.1)
  end
  @test teststatus == true

  teststatus = true
  for p in 0:0.1:10
    teststatus &= isapprox(mwls_diff(cll, p, 1), cos(p), rtol=0.1)
  end
  @test teststatus == true

  # not enough data for approximation, zero returned by default
  @test isapprox(mwls_diff(m, -0.9, 1), 0, atol=0.1)
  @test isapprox(mwls_diff(cll, -0.9, 1), 0, atol=0.1)
  @test isapprox(mwls_diff(m, 10.9, 1), 0, atol=0.1)
  @test isapprox(mwls_diff(cll, 10.9, 1), 0, atol=0.1)
end

