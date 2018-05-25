@testset "2d derivative" begin
  l = collect(-1.5:0.2:1.5);
  ins = transpose(hcat(collect.(collect(product(l, l))[:])...));
  func = (x, y) -> exp(-(x^2 + y^2));
  dfunc = (x, y) -> 4 * x * y * exp(-(x^2 + y^2));
  fs = [func(ins[i, 1], ins[i, 2]) for i in 1:size(ins, 1)];
  dfs = [dfunc(ins[i, 1], ins[i, 2]) for i in 1:size(ins, 1)];

  m = mwls_kd(hcat(ins, fs), 1, (d, eps) -> exp(-d^2), maxDegree = 5)
  cll = mwls_cll(ins, fs, 1, (d, eps) -> exp(-d^2), maxDegree = 5)

  teststatus = true
  for i in 1:size(ins, 1)
    teststatus &= isapprox(mwls_diff(m, [ins[i, 1], ins[i, 2]], (1, 1)),
                           dfunc(ins[i, 1], ins[i, 2]),
                           atol = 0.1)
    if !teststatus
      @show ins[i, :]
      @show mwls_diff(m, [ins[i, 1], ins[i, 2]], (1, 1))
      @show dfunc(ins[i, 1], ins[i, 2])
      break
    end
  end
  @test teststatus == true

  teststatus = true
  for i in 1:size(ins, 1)
    teststatus &= isapprox(mwls_diff(cll, [ins[i, 1], ins[i, 2]], (1, 1)),
                           dfunc(ins[i, 1], ins[i, 2]),
                           atol = 0.1)
    if !teststatus
      @show ins[i, :]
      @show mwls_diff(m, [ins[i, 1], ins[i, 2]], (1, 1))
      @show dfunc(ins[i, 1], ins[i, 2])
      break
    end
  end
  @test teststatus == true
end

