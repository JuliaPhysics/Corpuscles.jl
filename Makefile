doc:
	julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate(); include("docs/make.jl");'
tests: 
	julia --project=. test/runtests.jl
format:
	julia --project=@runic --startup-file=no -e 'using Pkg; Pkg.add("Runic")'
	git ls-files -z -- '*.jl' | xargs -0 --no-run-if-empty julia --project=@runic --startup-file=no -e 'using Runic; exit(Runic.main(ARGS))' -- --inplace
