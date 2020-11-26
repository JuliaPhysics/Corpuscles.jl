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


    """
    Runs the given is/has-function on a set of candidates and also checks every
    other noncandidate in PDGIDs for the opposite outcome.
    """
    function check_candidates(f, candidates)
        noncandidates = setdiff(Set(instances(PDGIDS)), Set(candidates))
        for candidate ∈ candidates
            @test f(candidate)
        end
        for noncandidate ∈ noncandidates
            @test !f(noncandidate)
        end
    end


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
        candidates = [DQuark, UQuark, SQuark, CQuark, BQuark, TQuark, BPrimeQuark, TPrimeQuark]
        check_candidates(isquark, candidates)
    end

    @testset "islepton" begin
        candidates = [Electron, Positron, Muon, AntiMuon, Tau, TauPrime, Nu_e, NuBar_tau, AntiElectronStar]
        check_candidates(islepton, candidates)
    end

    @testset "ismeson" begin
        candidates = [jpsi, psi_2S, Upsilon_1S, Upsilon_4S, Pi0, PiPlus, eta,
                  eta_prime, a_0_1450_plus, KL, KS, KMinus, phi, omega,
                  rho_770_minus, K1_1270_0, K1_1400_0, rho_1700_0, a2_1320_minus,
                  omega_3_1670, f_4_2300, D0, DPlus, DsPlus, B0, BPlus, Bs,
                  BcPlus, Pi0TC, PiMinusTC, T0, Reggeon, Pomeron, Odderon,
                  RPlus_TTildeDbar, R0_GTildeG]
        check_candidates(ismeson, candidates)
    end

    @testset "isbaryon" begin
        candidates = [Proton, AntiNeutron, HydrogenNucleus, Lambda, Sigma0,
                      SigmaPlus, SigmaMinus, Xi0, AntiXiMinus, OmegaMinus,
                      LcPlus, Lb, LtPlus, RPlusPlus_GTildeUUU,
                      UCbarCUDPentaquark, AntiUCbarCUDPentaquark]
        check_candidates(isbaryon, candidates)
    end

    @testset "ishadron" begin
        @test all(ishadron(p) == (ismeson(p) || isbaryon(p)) for p in instances(PDGIDS))
    end

    @testset "ispentaquark" begin
        candidates = [UCbarCUDPentaquark, AntiUCbarCUDPentaquark]
        check_candidates(ispentaquark, candidates)
    end

    @testset "isgaugebosonorhiggs" begin
        candidates = [Gluon, Photon, Z0, WMinus, HiggsBoson, ZPrime, Graviton]
        check_candidates(isgaugebosonorhiggs, candidates)
    end

    @testset "issmgaugebosonorhiggs" begin
        candidates = [Gluon, Photon, Z0, WMinus, HiggsBoson]
        check_candidates(issmgaugebosonorhiggs, candidates)
    end

    @testset "isgeneratorspecific" begin
        candidates = [AntiCHadron]
        check_candidates(isgeneratorspecific, candidates)
    end

    @testset "isspecial" begin
        candidates = [Graviton, Reggeon, Pomeron, Odderon, AntiCHadron]
        check_candidates(isspecial, candidates)
    end

    @testset "isnucleus" begin
        candidates = [Proton, AntiNeutron, HydrogenNucleus, Carbon12]
        check_candidates(isnucleus, candidates)
    end

    @testset "isdiquark" begin
        candidates = [DD1, SD0]
        check_candidates(isdiquark, candidates)
    end

    @testset "isRhadron" begin
        candidates = [RPlus_TTildeDbar, R0_GTildeG, RPlusPlus_GTildeUUU]
        check_candidates(isRhadron, candidates)
    end

    @testset "isQball" begin
        candidates = [QBall1, QBall2]
        check_candidates(isQball, candidates)
    end

    @testset "isdyon" begin
        candidates = [DyonSameMagElecChargeSign, DyonOppositeMagElecChargeSign]
        check_candidates(isdyon, candidates)
    end

    @testset "isSUSY" begin
        candidates = [Gluino, Gravitino, STildeL, CTildeR]
        check_candidates(isSUSY, candidates)
    end

    @testset "istechnicolor" begin
        candidates = [Pi0TC, PiMinusTC]
        check_candidates(istechnicolor, candidates)
    end

    @testset "iscompositequarkorlepton" begin
        candidates = [UQuarkStar, AntiElectronStar]
        check_candidates(iscompositequarkorlepton, candidates)
    end

    @testset "hasdown" begin
        f = hasdown
        @test f(Photon) == false
        @test f(Gluon) == false
        @test f(Electron) == false
        @test f(AntiMuon) == false
        @test f(jpsi) == false
        @test f(Upsilon_1S) == false
        @test f(PiPlus) == true
        @test f(KMinus) == false
        @test f(D0) == false
        @test f(DPlus) == true
        @test f(DsPlus) == false
        @test f(B0) == true
        @test f(Bs) == false
        @test f(BcPlus) == false
        @test f(Proton) == true
        @test f(LcPlus) == true
        @test f(Lb) == true
        @test f(DD1) == true
        @test f(SD0) == true
        @test f(Invalid1) == false
        @test f(Invalid2) == false
    end

    @testset "hasup" begin
        f = hasup
        @test f(Photon) == false
        @test f(Gluon) == false
        @test f(Electron) == false
        @test f(AntiMuon) == false
        @test f(jpsi) == false
        @test f(Upsilon_1S) == false
        @test f(PiPlus) == true
        @test f(KMinus) == true
        @test f(D0) == true
        @test f(DPlus) == false
        @test f(DsPlus) == false
        @test f(B0) == false
        @test f(Bs) == false
        @test f(BcPlus) == false
        @test f(Proton) == true
        @test f(LcPlus) == true
        @test f(Lb) == true
        @test f(DD1) == false
        @test f(SD0) == false
        @test f(Invalid1) == false
        @test f(Invalid2) == false
    end

    @testset "hasstrange" begin
        f = hasstrange
        @test f(Photon) == false
        @test f(Gluon) == false
        @test f(Electron) == false
        @test f(AntiMuon) == false
        @test f(jpsi) == false
        @test f(Upsilon_1S) == false
        @test f(PiPlus) == false
        @test f(KMinus) == true
        @test f(D0) == false
        @test f(DPlus) == false
        @test f(DsPlus) == true
        @test f(B0) == false
        @test f(Bs) == true
        @test f(BcPlus) == false
        @test f(Proton) == false
        @test f(LcPlus) == false
        @test f(Lb) == false
        @test f(DD1) == false
        @test f(SD0) == true
        @test f(Invalid1) == false
        @test f(Invalid2) == false
    end

    @testset "hascharm" begin
        candidates = [jpsi, psi_2S, D0, DPlus, DsPlus, BcPlus, LcPlus, UCbarCUDPentaquark, AntiUCbarCUDPentaquark]
        check_candidates(hascharm, candidates)
    end

    @testset "hasbottom" begin
        candidates = [Upsilon_1S, Upsilon_4S, B0, BPlus, Bs, BcPlus, Lb]
        check_candidates(hasbottom, candidates)
    end

    @testset "hastop" begin
        candidates = [T0, LtPlus]
        check_candidates(hastop, candidates)
    end

    @testset "hasfundamentalanti" begin
        candidates = [WMinus, Electron, Positron, Muon, AntiMuon, Tau, TauPrime,
                      Nu_e, NuBar_tau, DQuark, UQuark, SQuark, CQuark, BQuark,
                      TQuark, BPrimeQuark, TPrimeQuark, UQuarkStar,
                      AntiElectronStar, STildeL, CTildeR, AntiCHadron]
        check_candidates(hasfundamentalanti, candidates)
    end

    @testset "A" begin
        candidates = Dict(Proton => 1, AntiNeutron => 1, HydrogenNucleus => 1, Carbon12 => 12)
        for (candidate, value) ∈ candidates
            @test A(candidate) == value
        end
        noncandidates = setdiff(Set(keys(candidates)), Set(instances(PDGIDS)))
        for noncandidate in noncandidates
            @test isnothiung(A(candidate))
        end
    end

    @testset "Z" begin
        candidates = Dict(Proton => 1, AntiNeutron => 0, HydrogenNucleus => 1, Carbon12 => 6)
        for (candidate, value) ∈ candidates
            @test Z(candidate) == value
        end
        noncandidates = setdiff(Set(keys(candidates)), Set(instances(PDGIDS)))
        for noncandidate in noncandidates
            @test isnothiung(Z(candidate))
        end
    end

    @testset "isvalid" begin
        f = isvalid
        candidates = [Photon, Gluon, Electron, AntiMuon, jpsi, Upsilon_1S,
                      PiPlus, KMinus, D0, DPlus, DsPlus, B0, Bs, BcPlus, Proton,
                      LcPlus, Lb, DD1, SD0]
        noncandidates = [Invalid1, Invalid2]
        for candidate ∈ candidates
            @test f(Corpuscles.pdgid(candidate))
        end
        for noncandidate in noncandidates
            @test !f(Corpuscles.pdgid(noncandidate))
        end
    end

    @testset "charge" begin
        @test charge(Photon) == 0
        @test charge(Gluon) == 0
        @test charge(Electron) == -1
        @test charge(AntiMuon) == +1
        @test charge(jpsi) == 0
        @test charge(Upsilon_1S) == 0
        @test charge(PiPlus) == +1
        @test charge(KMinus) == -1
        @test charge(D0) == 0
        @test charge(DPlus) == +1
        @test charge(DsPlus) == +1
        @test charge(B0) == 0
        @test charge(Bs) == 0
        @test charge(BcPlus) == +1
        @test charge(Proton) == +1
        @test charge(LcPlus) == +1
        @test charge(Lb) == 0
        @test charge(DD1) == -2 / 3
        @test charge(SD0) == -2 / 3
        @test Corpuscles.isnothing(charge(Invalid1))
        @test Corpuscles.isnothing(charge(Invalid2))
    end

    @testset "threecharge" begin
        @test threecharge(Photon) == 0
        @test threecharge(Electron) == -3
        @test threecharge(jpsi) == 0
        @test threecharge(Upsilon_1S) == 0
        @test threecharge(KMinus) == -3
        @test threecharge(D0) == 0
        @test threecharge(Proton) == +3
        @test threecharge(LcPlus) == +3
        @test threecharge(Lb) == 0
        @test threecharge(DD1) == -2
        @test threecharge(SD0) == -2
        @test Corpuscles.isnothing(threecharge(Invalid1))
        @test Corpuscles.isnothing(threecharge(Invalid2))
    end

    @testset "jspin" begin
        # TODO
    end

    @testset "J" begin
        # TODO
    end
end
