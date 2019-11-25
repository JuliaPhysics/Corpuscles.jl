using Corpuscles
using Test

@testset "Corpuscles.jl" begin

    catalog_files = map(basename, Corpuscles.available_catalog_files())
    @test "particle2019.csv" in catalog_files

    # Particle Identites
    # PDGID <-> GeantID 
    geant_id = GeantID(1)
    pdg_id = convert(PDGID, geant_id)
    @test isequal(22, pdg_id.value)
    geant_id = GeantID(3)
    pdg_id = convert(PDGID, geant_id)
    @test isequal(11, pdg_id.value)
    geant_id = GeantID(3)
    pdg_id = convert(PDGID, geant_id)
    @test !isequal(13, pdg_id.value)
    pdg_id = PDGID(22)
    geant_id = convert(GeantID, pdg_id)
    @test isequal(1, geant_id.value)
    pdg_id = PDGID(11)
    geant_id = convert(GeantID, pdg_id)
    @test isequal(3, geant_id.value)
    pdg_id = PDGID(13)
    geant_id = convert(GeantID, pdg_id)
    @test !isequal(1, geant_id.value)
end
