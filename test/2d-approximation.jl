@testset "2d approximation" begin
  l = collect(-2:0.2:2);
  ins = transpose(hcat(collect.(collect(product(l, l))[:])...));
  func = (x, y) -> sin(x) + exp(-y^2)
  fs = [func(ins[i, 1], ins[i, 2]) for i in 1:size(ins, 1)]

  m = mwls_kd(ins, fs, 2, (d, eps) -> exp(-d^2), maxDegree = 3)
  cll = mwls_cll(hcat(ins, fs), 2, (d, eps) -> exp(-d^2), maxDegree = 3)
  teststatus = true
  for i in size(ins, 1)
    teststatus &= isapprox(m([ins[i, 1], ins[i, 2]]),
                           func(ins[i, 1], ins[i, 2]),
                           atol = 0.2)
  end
  @test teststatus == true

  for i in size(ins, 1)
    teststatus &= isapprox(cll([ins[i, 1], ins[i, 2]]),
                           func(ins[i, 1], ins[i, 2]),
                           atol = 0.3)
  end
  @test teststatus == true
end


