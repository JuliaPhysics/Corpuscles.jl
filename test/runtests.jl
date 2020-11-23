using Test
using InteractiveUtils
using Corpuscles
using Unitful

const DATA_DIR = joinpath(@__DIR__, "../data")

include("particles.jl")


@testset "inits" begin
    for T in vcat(map(subtypes, [Signed, Unsigned])...)
        p = Particle(T(11))
        @test PDGID(11) == p.pdgid
    end
end


@testset "pdg digits" begin
    p = Particle(9020113)
    @test p.pdgid.Nj == 3
    @test p.pdgid.Nq3 == 1
    @test p.pdgid.Nq2 == 1
    @test p.pdgid.Nq1 == 0
    @test p.pdgid.Nl == 2
    @test p.pdgid.Nr == 0
    @test p.pdgid.N == 9
    @test p.pdgid.N8 == 0
    @test p.pdgid.N9 == 0
    @test p.pdgid.N10 == 0

    # digits should be positive
    p = Particle(-321)
    @test p.pdgid.Nj == 1
    @test p.pdgid.Nq3 == 2
    @test p.pdgid.Nq2 == 3
end

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

    @test PDGID(13) === convert(PDGID, PDGID(13))
    
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
        @test occursin("2222212", msg)
    end

    for particle in particles()
        if particle.pdgid in [PDGID(-9000321), PDGID(9000321), PDGID(9020113),
                              PDGID(-9000311), PDGID(9020213), PDGID(9000311),
                              PDGID(-9020213)]
            continue
        end
        @test convert(Geant3ID, particle.pdgid) !== nothing
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
    g = 91u"cm"
    h = 91u"s"

    @test a == b
    @test !(a == c)
    @test d < a
    @test !(e < a)
    @test !(e < c)
    @test g < c
    @test !(f < h)
    @test isapprox(g, e)
    @test isapprox(g, a)
    @test isapprox(a, b)
    @test isapprox(a, c)
    @test isapprox(a, e)
    @test !isapprox(a, d)

    try
        e < f
        @test false
    catch e
        @test e isa Unitful.DimensionError
    end
end

@testset "show and print" begin
    io = IOBuffer(write=true)
    show(io, Particle(1))
    seekstart(io)
    @test "Particle(1) d" == String(read(io))

    io = IOBuffer(write=true)
    print(io, Particle(1))
    seekstart(io)
    output = String(read(io))
    @test occursin(r"PDG ID:\s*1", output)
    @test occursin(r"Name:\s*d\n", output)
end


@testset "helpers" begin
    Corpuscles.use_catalog_file(joinpath(DATA_DIR, "particle2020.csv"))
    for _particles ∈ [particles(), [p.pdgid for p in particles()], [p.pdgid.value for p in particles()]]
        @test sum(map(isstandard, _particles)) > 0
        @test sum(map(isfundamental, _particles)) > 0
        @test sum(map(Corpuscles.isstandard, _particles)) == 610
        @test sum(map(Corpuscles.fundamentalid, _particles)) == 413

        @test sum(map(isquark, _particles)) == 12
        @test sum(map(islepton, _particles)) == 16
        @test sum(map(ismeson, _particles)) == 234
        @test sum(map(isbaryon, _particles)) == 292
        @test sum(map(ishadron, _particles)) == 526
        @test sum(map(isRhadron, _particles)) == 1
        @test sum(map(isSUSY, _particles)) == 0
        @test sum(map(ispentaquark, _particles)) == 0
        @test sum(map(isgaugebosonorhiggs, _particles)) == 6
        @test sum(map(issmgaugebosonorhiggs, _particles)) == 6
        @test sum(map(isdyon, _particles)) == 0
        @test sum(map(isnucleus, _particles)) == 4
        @test sum(map(isdiquark, _particles)) == 50
        @test sum(map(istechnicolor, _particles)) == 0
        @test sum(map(iscompositequarkorlepton, _particles)) == 0
        @test sum(map(isgeneratorspecific, _particles)) == 0
        @test sum(map(isspecial, _particles)) == 0
        @test sum(map(isQball, _particles)) == 0
        @test sum(map(hasfundamentalanti, _particles)) == 30

        @test sum(map(hasdown, _particles)) == 328
        @test sum(map(hasup, _particles)) == 346
        @test sum(map(hascharm, _particles)) == 107
        @test sum(map(hasstrange, _particles)) == 257
        @test sum(map(hasbottom, _particles)) == 68
        @test sum(map(hastop, _particles)) == 0
    end

    @testset "isquark" begin
        quarks = [DQuark, UQuark, SQuark, CQuark, BQuark, TQuark, BPrimeQuark, TPrimeQuark]
        nonquarks = setdiff(Set(instances(PDGIDS)), Set(quarks))
        @test all(isquark(p) for p in quarks)
        @test all(!isquark(p) for p in nonquarks)
    end

    @testset "islepton" begin
        f = islepton
        candidates = [Electron, Positron, Muon, AntiMuon, Tau, TauPrime, Nu_e, NuBar_tau, AntiElectronStar]
        noncandidates = setdiff(Set(instances(PDGIDS)), Set(candidates))
        @test all(f(p) for p in candidates)
        @test all(!f(p) for p in noncandidates)
    end

    @testset "ismeson" begin
        f = ismeson
        candidates = [jpsi, psi_2S, Upsilon_1S, Upsilon_4S, Pi0, PiPlus, eta,
                  eta_prime, a_0_1450_plus, KL, KS, KMinus, phi, omega,
                  rho_770_minus, K1_1270_0, K1_1400_0, rho_1700_0, a2_1320_minus,
                  omega_3_1670, f_4_2300, D0, DPlus, DsPlus, B0, BPlus, Bs,
                  BcPlus, Pi0TC, PiMinusTC, T0, Reggeon, Pomeron, Odderon,
                  RPlus_TTildeDbar, R0_GTildeG]
        noncandidates = setdiff(Set(instances(PDGIDS)), Set(candidates))
        @test all(f(p) for p in candidates)
        @test all(!f(p) for p in noncandidates)
    end
end
