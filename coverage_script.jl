Pkg.add("Coverage")
cd(Pkg.dir("MovingWeightedLeastSquares"))
using Coverage
covered_lines, total_lines = get_summary(process_folder())
percentage = covered_lines / total_lines * 100
println("($(percentage)%) covered")
