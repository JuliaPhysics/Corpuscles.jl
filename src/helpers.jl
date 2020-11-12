# Helper functions


isquark(p::Particle) = 1 <= abs(p.pdgid.value) <= 8
