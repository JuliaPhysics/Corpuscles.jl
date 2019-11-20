using Particles
using Test

@testset "Particles.jl" begin

    catalog_files = map(basename, Particles.available_catalog_files())
    @test "particle2019.csv" in catalog_files

end
