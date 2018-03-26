using MovingWeightedLeastSquares

ins = [1. 1.; 1. -1.0; -1.0 1.0; -1.0 -1.0; 0. 0.; 1. 0.; -1. 0.; 0. 1.; 0. -1]
outs = [1., -0.5, 1., 1., -1., 0., 0., 0., 0.]

m = mwls(ins, outs, 2, 1., (d, e) -> exp(-d^2))

println(MovingWeightedLeastSquares.diff(m, [0., 0.], (1, 2); dist = 10.))
