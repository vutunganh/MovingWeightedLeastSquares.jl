@testset "cell linked list range test" begin
  range = 1000
  dataset = range * rand(2, 100) - range / 2
  cll = CellLinkedList(dataset, 3.)

  teststatus = true
  for c in 1:100
    pt = 1.1 * range * rand(2) - range * 1.1 / 2
    r = range / 3
    cllres = cll_inrange(cll, pt, r)
    naiveres = [p for p in 1:size(dataset, 2) if norm(pt - dataset[:, p]) < r + 1e-9]
    teststatus &= sort(cllres) == sort(naiveres)
    if !teststatus
      @show dataset
      @show pt
      @show r
      break
    end
  end
  @test teststatus == true
end
