var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#MovingWeightedLeastSquares.jl-1",
    "page": "Home",
    "title": "MovingWeightedLeastSquares.jl",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Overview-1",
    "page": "Home",
    "title": "Overview",
    "category": "section",
    "text": "MovingWeightedLeastSquares.jl is a package that provides an implementation of the moving weighted least squares method. A short description of this method by Andy Nealen can be found here.Let theta(d) mathbbR^+ rightarrow mathbbR^+ be a weighting function of the method. Very often the there will be an varepsilon in mathbbR, such that forall delta  varepsilon theta(delta) = 0. Whenever we use varepsilon or EPS in this document, we mean the cutoff distance for the weighting function. An example of a good weighting function is theta(d) = exp(d^2  a^2), where a is the average distance between sample input data.Point data type is an alias to Vector{T} where {T <: Real}.Interaction with this package is done mostly via structure MwlsObject and its subclasses MwlsKdObject, MwlsCllObject and MwlsNaiveObject. This interface is similar to the interface of \"interpolations objects\" from Interpolations.jl.The difference between the subclasses is the solution of the range search problem. MwlsKdObject solves the range search problem by using a k-d tree created by Kristoffer Carlsson, see NearestNeighbors.jl. MwlsCllObject solves the range search problem by using a cell linked list, which is implemented in this package. MwlsNaiveObject solves the range search problem naively.TL;DR: use anything but MwlsNaiveObject."
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "At the moment this package can be installed by manually cloning the packagePkg.clone(\"https://github.com/vutunganh/MovingWeightedLeastSquares.jl\")The package has been merged to Julia\'s METADATA, but I haven\'t tried installing it that way (using Pkg.add) yet."
},

{
    "location": "constructors.html#",
    "page": "Constructors",
    "title": "Constructors",
    "category": "page",
    "text": ""
},

{
    "location": "constructors.html#MovingWeightedLeastSquares.mwls_cll",
    "page": "Constructors",
    "title": "MovingWeightedLeastSquares.mwls_cll",
    "category": "function",
    "text": "mwls_cll(inputs::Array{T, N}, outputs::Array{U}, EPS::Real, weightfunc::Function) where {T <: Real, U <: Real, N}\nmwls_cll(inputs::Array{T, N}, outputs::Array{U}, EPS::Real, weightfunc::Function; maxDegree::Int = 2) where {T <: Real, U <: Real, N}\n\nCreates MwlsCllObject from sample input and sample output data, the cutoff distance ε and a weighting function θ.\n\nArguments\n\ninputs: a 2d array of input points where each point is on a single row,\noutputs: a 2d array or a vector of output scalars where each output is on a single row,\nEPS::Real: ε of the method (cell edge length and the default distance for range search),\nweightfunc::Function: weighting function θ of the method. It should be in form (distance between two vectors, EPS) -> Float64.\n\nKeyword arguments\n\nmaxDegree::Int: the maximal degree of polynomials used for approximation, 2 by default.\n\n\n\nmwls_cll(input::Array{T, 2}, EPS::Real, weightfunc::Function) where {T <: Real}\nmwls_cll(input::Array{T, 2}, EPS::Real, weightfunc::Function; outputDim::Int = 1, maxDegree::Int = 2) where {T <: Real}\n\nIn this mwls_cll function, the sample input and sample output data are passed in a single array. It is assumed that each pair of input and output is on a single row. Dimension of the output is specified with kwarg outputDim.\n\n\n\n"
},

{
    "location": "constructors.html#MovingWeightedLeastSquares.mwls_kd",
    "page": "Constructors",
    "title": "MovingWeightedLeastSquares.mwls_kd",
    "category": "function",
    "text": "mwls_kd(inputs::Array{T, N}, outputs::Array{U}, EPS::Real, weightfunc::Function) where {T <: Real, U <: Real, N}\nmwls_kd(inputs::Array{T, N}, outputs::Array{U}, EPS::Real, weightfunc::Function; leafsize::Int = 10, maxDegree::Int = 2) where {T <: Real, U <: Real, N}\n\nCreates MwlsKdObject from sample input and sample output data, the cutoff distance ε and a weighting function θ.\n\nArguments\n\ninputs: a 2d array or a vector of input points where each point is on a single row,\noutputs: a 2d array or a vector of output scalars where each output is on a single row,\nEPS::Real: ε of the method (default distance threshold for neighbor search),\nweightfunc::Function: weighting function of the method. It should be in form (distance, EPS) -> Float64.\n\nKeyword arguments\n\nleafSize::Int: the size of the leaves in the k-d-tree, 10 by default.\nmaxDegree::Int: the maximal degree of polynomials used for approximation, 2 by default.\n\n\n\nmwls_kd(input::Array{T, 2}, EPS::Real, weightfunc::Function) where {T <: Real}\nmwls_kd(input::Array{T, 2}, EPS::Real, weightfunc::Function; outputDim::Int = 1, leafSize::Int = 10, maxDegree::Int = 2) where {T <: Real}\n\nIn this mwls_kd function, the sample input and sample output data are passed in a single array. It is assumed that each pair of input and output is on a single row. Dimension of the output is specified with kwarg outputDim.\n\n\n\n"
},

{
    "location": "constructors.html#Constructors-1",
    "page": "Constructors",
    "title": "Constructors",
    "category": "section",
    "text": "Constructing MwlsObject directly is not recommended and helper functions described below should be used instead.mwls_cll\nmwls_kd"
},

{
    "location": "approximation.html#",
    "page": "Approximation",
    "title": "Approximation",
    "category": "page",
    "text": ""
},

{
    "location": "approximation.html#Approximation-1",
    "page": "Approximation",
    "title": "Approximation",
    "category": "section",
    "text": "Let obj be a subclass of MwlsObject. The function call obj(pt::Point) returns the approximated function value at pt. This calls the approximate function."
},

{
    "location": "approximation.html#Example-1",
    "page": "Approximation",
    "title": "Example",
    "category": "section",
    "text": "Let\'s initialize a dataset with input data xs and output data fs.using MovingWeightedLeastSquares # hide\nxs = collect(-2:0.1:2);\nfs = [sin(x) for x in xs];Now let\'s construct an approximation object obj. Let\'s choose a weighting function theta(d) = exp(d^2).obj = mwls_kd(xs, fs, 0.5, (d, e) -> (exp(-d^2)));The approximation at 1 can be obtained by usingobj(1)If a different range of neighbor data is needed, then useobj(1, 1)"
},

{
    "location": "approximation.html#MovingWeightedLeastSquares.approximate",
    "page": "Approximation",
    "title": "MovingWeightedLeastSquares.approximate",
    "category": "function",
    "text": "approximate(obj::MwlsObject, pt::Point)\napproximate(obj::MwlsObject, pt::Point; dist::Real = obj.EPS)\n\nThis calculates the approximated value at pt for each dimension of output data. The approximated value is returned.\n\n\n\n"
},

{
    "location": "approximation.html#MovingWeightedLeastSquares.mwls_coefficients",
    "page": "Approximation",
    "title": "MovingWeightedLeastSquares.mwls_coefficients",
    "category": "function",
    "text": "mwls_coefficients(obj::MwlsObject, pt::Point, dist::Real)\n\nCalculates the coefficients of the linear combination of polynomials used for approximation. This is done for each dimension of output data.\n\nnote: Note\nIf the matrix in the system of linear equations used to find the coefficients is singular, then zero coefficients are returned!\n\n\n\n"
},

{
    "location": "approximation.html#Relevant-documentation-1",
    "page": "Approximation",
    "title": "Relevant documentation",
    "category": "section",
    "text": "approximate\nmwls_coefficients"
},

{
    "location": "approximation.html#Approximation-of-derivative-1",
    "page": "Approximation",
    "title": "Approximation of derivative",
    "category": "section",
    "text": ""
},

{
    "location": "approximation.html#Example-2",
    "page": "Approximation",
    "title": "Example",
    "category": "section",
    "text": "Approximation of first derivative on the dataset from the previous example can be done by usingmwls_diff(obj, 1, 1)Formally a tuple is required for specifying the orders of derivates for each variable. Let f be the approximated function. For an example a tuple (1, 2) calculatesfracpartial^3partial x_1 partial x_2^2f"
},

{
    "location": "approximation.html#MovingWeightedLeastSquares.mwls_diff",
    "page": "Approximation",
    "title": "MovingWeightedLeastSquares.mwls_diff",
    "category": "function",
    "text": "mwls_diff(obj::MwlsNaiveObject, inPt::Real, dirs::Integer; dist::Real = obj.EPS)\n\n\n\nmwls_diff(obj::MwlsNaiveObject, inPt::Point, dirs::Integer; dist::Real = obj.EPS)\n\n\n\nmwls_diff(obj::MwlsNaiveObject, inPt::Real, dirs::NTuple{N, Integer}; dist::Real = obj.EPS) where {N}\n\n\n\nmwls_diff(obj::MwlsObject, inPt::Real, dirs::Integer; dist::Real = obj.EPS)\n\n\n\nmwls_diff(obj::MwlsObject, inPt::Real, dirs::NTuple{N, Integer}; dist::Real = obj.EPS) where {N}\n\n\n\nmwls_diff(obj::MwlsObject, inPt::Point, dirs::Integer; dist::Real = obj.EPS)\n\n\n\nmwls_diff(obj::MwlsObject, inPt::Point, dirs::NTuple{N, Integer}; dist::Real = obj.EPS)\n\nCalculates the approximated derivative at inPt, where x[i] is differentiated dirs[i] times.\n\n\n\n"
},

{
    "location": "approximation.html#MovingWeightedLeastSquares.mwls_diff_polys",
    "page": "Approximation",
    "title": "MovingWeightedLeastSquares.mwls_diff_polys",
    "category": "function",
    "text": "mwls_diff_polys(obj::MwlsObject, inPt::Point, dirs::NTuple{N, Integer}; dist::Real = obj.EPS) where {N}\n\nPolynomials created by mwls_coefficients are differentiated according to dirs. For an example if f is the polynomial used for approximation, then dirs = (12) returns fracpartial^3partial x_1 partial x_2^2f.\n\n\n\n"
},

{
    "location": "approximation.html#Relevant-documentation-2",
    "page": "Approximation",
    "title": "Relevant documentation",
    "category": "section",
    "text": "mwls_diff\nmwls_diff_polys"
},

]}
