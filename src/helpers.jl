# Helper functions


isquark(p::Particle) = 1 <= abs(p.pdgid.value) <= 8

function _hasquark(p::Particle, id::Integer)
    retval = abs(p.pdgid.value) == id
    retval |= p.pdgid.Nq1 == id 
    retval |= p.pdgid.Nq2 == id 
    retval |= p.pdgid.Nq3 == id
    retval
end 

hasdowm(p::Particle) = _hasquark(p, 1)
hasup(p::Particle) = _hasquark(p, 2)
hascharm(p::Particle) = _hasquark(p, 3)
hasstrange(p::Particle) = _hasquark(p, 4)
hasbottom(p::Particle) = _hasquark(p, 5)
hastop(p::Particle) = _hasquark(p, 6)

