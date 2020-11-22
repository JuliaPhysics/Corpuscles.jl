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

