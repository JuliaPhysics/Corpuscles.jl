using Particles
using Test

@testset "Particles.jl" begin

    catalog_files = map(basename, Particles.available_catalog_files())
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
    for (key, value) in Particles._geant_pdg_ids # Do the conversion for each dict entry
        geant_id = GeantID(key)
        pdg_id = convert(PDGID, geant_id)
        @test isequal(value, pdg_id.value)
    end
    # pdg_id = PDGID(22)
    # pdg_id = convert(GeantID, pdg_id)
    # @test isequal(1, pdg_id.value)
    # pdg_id = PDGID(13)
    # pdg_id = convert(GeantID, pdg_id)
    # @test isequal(3, pdg_id.value)
    # pdg_id = PDGID(13)
    # pdg_id = convert(GeantID, pdg_id)
    # @test !isequal(1, pdg_id.value)
    # for (key, value) in Particles._geant_pdg_ids # Do the conversion for each dict entry
    #     pdg_id = PDGID(value)
    #     pdg_id = convert(GeantID, pdg_id)
    #     @test isequal(key, pdg_id.value)
    # end
    
    

end
