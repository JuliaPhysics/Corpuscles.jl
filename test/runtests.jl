using Test
using Corpuscles
using Unitful

const DATA_DIR = joinpath(@__DIR__, "../data")

@testset "conversions" begin

    # Particle Identites
    # PDGID <-> Geant3ID 
    geant_id = Geant3ID(1)
    pdg_id = convert(PDGID, geant_id)
    @test isequal(22, pdg_id.value)
    geant_id = Geant3ID(3)
    pdg_id = convert(PDGID, geant_id)
    @test isequal(11, pdg_id.value)
    geant_id = Geant3ID(3)
    pdg_id = convert(PDGID, geant_id)
    @test !isequal(13, pdg_id.value)
    pdg_id = PDGID(22)
    geant_id = convert(Geant3ID, pdg_id)
    @test isequal(1, geant_id.value)
    pdg_id = PDGID(11)
    geant_id = convert(Geant3ID, pdg_id)
    @test isequal(3, geant_id.value)
    pdg_id = PDGID(13)
    geant_id = convert(Geant3ID, pdg_id)
    @test !isequal(1, geant_id.value)
    
    try
        convert(PDGID, Geant3ID(5000))
        @test false
    catch e
        @test e isa Corpuscles.IDException
        msg = sprint(showerror, e)
        @test msg == "ParticleID Error: No corresponding PDGID for Geant3ID(5000) found!"
    end
    
    try
        convert(Geant3ID, PDGID(2222212))
        @test false
    catch e
        @test e isa Corpuscles.IDException
        msg = sprint(showerror, e)
        @test msg == "ParticleID Error: No corresponding Geant3ID for PDGID(2222212) found!"
    end

end

@testset "catalogs" begin
    catalog_files = map(basename, Corpuscles.available_catalog_files())
    @test "particle2008.csv" in catalog_files
    @test "particle2018.csv" in catalog_files
    @test "particle2019.csv" in catalog_files

    Corpuscles.use_catalog_file(joinpath(DATA_DIR, "particle2008.csv"))
    @test 0.054u"MeV" == Particle(553).width.value
    Corpuscles.use_catalog_file(joinpath(DATA_DIR, "particle2019.csv"))
    @test 0.05402u"MeV" == Particle(553).width.value

end

@testset "comparison" begin
    a = Corpuscles.MeasuredValue(1u"m", 10u"cm", 10u"cm")
    b = Corpuscles.MeasuredValue(1u"m", 10u"cm", 10u"cm")
    c = Corpuscles.MeasuredValue(1u"m", 5u"cm", 5u"cm")
    d = Corpuscles.MeasuredValue(0.8u"m", 5u"cm", 5u"cm")
    e = Corpuscles.MeasuredValue(0.9u"m", 5u"cm", 5u"cm")
    f = Corpuscles.MeasuredValue(90u"s", 5u"s", 5u"s")
    g = 90u"cm"

    @test a == b
    @test !(a == c)
    @test d < a
    @test !(e < a)
    @test !(e < c)
    @test g < c
    @test g == e

    try
        e < f
        @test false
    catch e
        @test e isa Unitful.DimensionError
    end
end
