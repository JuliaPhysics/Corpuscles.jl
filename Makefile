doc:
	julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate(); include("docs/make.jl");'
tests: 
	julia --project=. test/runtests.jl
format:
	julia -e 'using Pkg; Pkg.activate(temp=true); Pkg.add("JuliaFormatter"); using JuliaFormatter; format(".")'
