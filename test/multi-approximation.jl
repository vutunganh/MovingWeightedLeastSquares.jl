using MovingWeightedLeastSquares

ins = [1. 1.; 1. -1.0; -1.0 1.0; -1.0 -1.0; 0. 0.; 1. 0.; -1. 0.; 0. 1.; 0. -1]
outs = [1., -0.5, 1., 1., -1., 0., 0., 0., 0.]

m = mwls(ins, outs, 1., (d, e) -> exp(-d^2))

println(m(0, 10))
println(m(0, 1))
