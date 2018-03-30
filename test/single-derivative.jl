using MovingWeightedLeastSquares

xs = [j * 0.1 for j in 0:100]
fs = [sin(x) for x in xs]

m = mwls(xs, fs, 1., (d, e) -> exp(-d^2))

println(MovingWeightedLeastSquares.diff(m, [-1.0], (1)))
