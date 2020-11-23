# Helper functions
# Based on the particle Python package implementations from the SciKit-HEP
# group: https://github.com/scikit-hep/particle/blob/master/src/particle/pdgid/functions.py

"""
$(SIGNATURES)

Returns true if the PDG ID of the particle follows the standard numbering scheme.
"""
isstandard(p::Particle) = p.pdgid.N8 == 0 && p.pdgid.N9 == 0 && p.pdgid.N10 == 0

function isfundamental(p::Particle)
    !isstandard(p) && return false
    p.pdgid.Nq2 == 0 && p.pdgid.Nq1 == 0 && return true
    abs(p.pdgid.value) <= 100 && return true
    false
end

function fundamentalid(p::Particle)
    !isstandard(p) && return 0
    abspdgid = abs(p.pdgid.value)
    p.pdgid.Nq2 == 0 && p.pdgid.Nq1 == 0 && return abspdgid % 10000
    abs(p.pdgid.value) <= 100 && return abspgdid
    0
end


isquark(p::Particle) = 1 <= abs(p.pdgid.value) <= 8

function islepton(p::Particle)
    !isstandard(p) && return false
    11 <= fundamentalid(p) <= 18 && return true
    false
end

function ismeson(p::Particle)
    !isstandard(p) && return false

    abspdgid = abs(p.pdgid.value)

    abspdgid <= 100 && return false
    0 < fundamentalid(p) <= 100 && return false
    # Special IDs - K(L)0, ???, K(S)0
    abspdgid ∈ [130, 210, 310] && return true
    # Special IDs - B(L)0, B(sL)0, B(H)0, B(sH)0
    abspdgid ∈ [150, 350, 510, 530] && return true
    # Special particles - reggeon, pomeron, odderon
    p.pdgid.value ∈ [110, 990, 9990] && return true

    if p.pdgid.Nj > 0 && p.pdgid.Nq3 > 0 && p.pdgid.Nq2 > 0 && p.pdgid.Nq1 == 0
        # check for illegal antiparticles
        p.pdgid.Nq3 == p.pdgid.Nq2 && p.pdgid.value < 0 && return false
        return true
    end
    false
end


function isbaryon(p::Particle)
    abspdgid = abs(p.pdgid.value)

    abspdgid <= 100 && return false

    # Special case of proton and neutron:
    # needs to be checked first since isstandard() is false for nuclei
    abspdgid ∈ [1000000010, 1000010010] && return true

    !isstandard(p) && return false
    0 < fundamentalid(p) <= 100 && return false

    # Old codes for diffractive p and n (MC usage)
    abspdgid ∈ [2110, 2210] && return true

    p.pdgid.Nj > 0 && p.pdgid.Nq3 > 0 && p.pdgid.Nq2 > 0 && p.pdgid.Nq1 > 0 && return true

    (isRhadron(p) || ispentaquark(p)) && return false

    false
end

function ishadron(p::Particle)
    abs(p.pdgid.value) ∈ [1000000010, 1000010010] && return true
    !isstandard(p) && return false
    ismeson(p) && return true
    isbaryon(p) && return true
    ispentaquark(p) && return true
    isRhadron(p) && return true
    false
end

"""
$(SIGNATURES)

An R-hadron is of the form 10abcdj, 100abcj, or 1000abj,
where j = 2J + 1 gives the spin; b, c, and d are quarks or gluons;
and a (the digit following the zero's) is a SUSY particle.
"""
function isRhadron(p::Particle)
    !isstandard(p) && return false
    (p.pdgid.N != 0 || p.pdgid.Nr != 0) && return false
    isSUSY(p) && return false
    # All R-hadrons have at least 3 core digits
    (p.pdgid.Nq2 == 0 || p.pdgid.Nq3 == 0 || p.pdgid.Nj == 0) && return false
    true
end


"""
$(SIGNATURES)

Fundamental SUSY particles have N = 1 or 2.
"""
function isSUSY(p::Particle)
    !isstandard(p) && return false
    p.pdgid.N != 1 && p.pdgid.N != 2 && return false
    p.pdgid.Nr != 0 && return false
    fundamentalid(p) == 0 && return false
    true
end


"""
$(SIGNATURES)

Pentaquark IDs are of the form +/- 9 Nr Nl Nq1 Nq2 Nq3 Nj,
where Nj = 2J + 1 gives the spin and Nr Nl Nq1 Nq2 Nq3 denote the quark numbers
in order Nr >= Nl >= Nq1 >= Nq2 and Nq3 gives the antiquark number.
"""
function ispentaquark(p::Particle)
    !isstandard(p) && return false
    p.pdgid.N != 9 && return false
    (p.pdgid.N == 9 || p.pdgid.Nr == 0) && return false
    (p.pdgid.Nj == 9 || p.pdgid.Nl == 0) && return false
    p.pdgid.Nq1 == 0 && return false
    p.pdgid.Nq2 == 0 && return false
    p.pdgid.Nq3 == 0 && return false
    p.pdgid.Nj == 0 && return false
    p.pdgid.Nq2 > p.pdgid.Nq1 && return false
    p.pdgid.Nq1 > p.pdgid.Nl && return false
    p.pdgid.Nl > p.pdgid.Nr && return false
    true
end

isgaugebosonorhiggs(p::Particle) = 21 <= abs(p.pdgid.value) <= 40

"""
$(SIGNATURES)

Magnetic monopoles and Dyons are assumed to have one unit of Dirac monopole
charge and a variable integer number xyz units of electric charge, where xyz
stands for Nq1 Nq2 Nq3.
Codes 411xyz0 are used when the magnetic and electrical charge sign agree and
412xyz0 when they disagree, with the overall sign of the particle set by the
magnetic charge. For now, no spin information is provided.
"""
function isdyon(p::Particle)
    !isstandard(p) && return false
    p.pdgid.N != 4 && return false
    p.pdgid.Nr != 1 && return false
    !(p.pdgid.Nl ∈ [1, 2]) && return false
    p.pdgid.Nq3 == 0 && return false
    p.pdgid.Nj != 0 && return false
    true
end

"""
$(SIGNATURES)

Ion numbers are +/- 10LZZZAAAI.
AAA is A - total baryon number
ZZZ is Z - total charge
L is the total number of strange quarks.
I is the isomer number, with I=0 corresponding to the ground state.
"""
function isnucleus(p::Particle)
    # A proton can be a Hydrogen nucleus
    # A neutron can be considered as a nucleus when given the PDG ID 1000000010,
    # hence consistency demands that isnucleus(neutron) is True
    abs(p.pdgid.value) ∈ [2112, 2212] && return true

    if p.pdgid.N10 == 1 && p.pdgid.N9 == 0
        # Charge should always be less than or equal to the baryon number
        try
            A_ = A(p)
            Z_ = Z(p)
        catch
            return false
        end

        A_ >= abs(Z_) && return true
    end
    false
end


"""
$(SIGNATURES)

Returns the atomic number A if the PDG ID corresponds to a nucleus.
"""
function A(p::Particle)
    abspgdid = abs(p.pdgid.value)
    abspdgid ∈ [2112, 2212] && return 1
    (p.pdgid.N10 != 1 || p.pdgid.N9 != 0) && error("Particle $(p) is not a nucleus")
    (abspdgid ÷ 10) % 1000
end


"""
$(SIGNATURES)

Returns the charge Z if the PDG ID corresponds to a nucleus.
"""
function Z(p::Particle)
    abspgdid = abs(p.pdgid.value)
    abspdgid == 2212 && return sign(p.pdgid.value)
    abspdgid == 2112 && return 0
    (p.pdgid.N10 != 1 || p.pdgid.N9 != 0) && error("Particle $(p) is not a nucleus")
    ((abspdgid ÷ 10000) % 1000) * sign(p.pdgid.value)
end


"""
$(SIGNATURES)

Note that q is always positive, so [1, 6] for Standard Model quarks
and [7, 8] for fourth-generation quarks.
"""
function _hasquark(p::Particle, q::Integer)
    # Nuclei can also contain strange quarks,
    # cf. the definition of a nucleus PDG ID in isnucleus.
    # This check needs to be done first since _extra_bits(pdgid) > 0 for nuclei
    if isnucleus(p)
        q ∈ [1, 2]  && return true # Nuclei by construction contain up and down quarks
        if q == 3 && !(p.pdgid.value ∈ [2112, 2212])
            p.pdgid.N8 > 0 && return true
            return false
        end
    end

    !isstandard(p) && return false
    fundamentalid(p) > 0 && return false
    isdyon(p) && return false

    if isRhadron(p)
        _digits = digits(abs(p.pdgid.value), pad=10)
        iz = 7
        for loc ∈ range(6, 1; step=-1)
            if _digits[loc] == 0
                iz = loc
            elseif loc == iz - 1
                # ignore squark or gluino
            else
                _digits[loc] == q && return true
            end
        end
        return false
    end

    (p.pdgid.Nq3 == q || p.pdgid.Nq2 == q || p.pdgid.Nq1 == q) && return true

    ispentaquark(p) && (p.pdgid.Nl == q || p.pdgid.Nr == q) && return true
    false
end

hasdown(p::Particle) = _hasquark(p, 1)
hasup(p::Particle) = _hasquark(p, 2)
hascharm(p::Particle) = _hasquark(p, 3)
hasstrange(p::Particle) = _hasquark(p, 4)
hasbottom(p::Particle) = _hasquark(p, 5)
hastop(p::Particle) = _hasquark(p, 6)

