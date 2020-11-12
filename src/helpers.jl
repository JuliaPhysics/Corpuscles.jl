# Helper functions


isquark(p::Particle) = 1 <= abs(p.pdgid.value) <= 8

function hastop(p::Particle)
    retval = abs(p.pdgid.value) == 6
    retval |= abs(p.pdgid.Nq1) == 6
    retval |= abs(p.pdgid.Nq2) == 6
    retval |= abs(p.pdgid.Nq3) == 6
    retval
end 
