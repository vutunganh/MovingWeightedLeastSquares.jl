using MovingWeightedLeastSquares

xs = [0.1 * i for i in 0:100]
fs = hcat([sin(x) for x in xs], [cos(x) for x in xs])
println(fs)

m = mwls(xs, fs, 0.5, (d, e) -> exp(-d^2))

println(m(1))
