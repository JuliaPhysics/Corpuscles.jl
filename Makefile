doc:
	julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate(); include("docs/make.jl");'
tests: 
	julia --project=. test/runtests.jl
format:
	julia --project=@runic -e 'using Pkg; Pkg.add("Runic")'
	julia --project=@runic -e 'using Runic; files = readlines(`git ls-files '\''*.jl'\''`); isempty(files) || exit(Runic.main(vcat(["--inplace"], files)))'
