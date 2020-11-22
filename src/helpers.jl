# Helper functions

"""
$(SIGNATURES)

Returns true if the PDG ID of the particle follows the standard numbering scheme.
"""
isstandard(p::Particle) = p.pdgid.N8 == 0 && p.pdgid.N9 == 0 && p.pdgid.N10 == 0

function isfundamental(p::Particle)
    !isstandard(p) && return false
    p.pdgid.Nq1 == 0 && p.pdgid.Nq2 == 0 && return true
    abs(p.pdgid.value) <= 100 && return true
    false
end

function fundamentalid(p::Particle)
    !isfundamental(p) && return 0
    p.pdgid.value % 10000
end


isquark(p::Particle) = 1 <= abs(p.pdgid.value) <= 8
islepton(p::Particle) = 1 <= abs(p.pdgid.value) <= 8

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

function _hasquark(p::Particle, id::Integer)
    retval = abs(p.pdgid.value) == id
    retval |= p.pdgid.Nq1 == id 
    retval |= p.pdgid.Nq2 == id 
    retval |= p.pdgid.Nq3 == id
    retval
end 

hasdown(p::Particle) = _hasquark(p, 1)
hasup(p::Particle) = _hasquark(p, 2)
hascharm(p::Particle) = _hasquark(p, 3)
hasstrange(p::Particle) = _hasquark(p, 4)
hasbottom(p::Particle) = _hasquark(p, 5)
hastop(p::Particle) = _hasquark(p, 6)

