# Helper Functions

There are many useful functions available to check additional properties of
particles, like their types, composition, charge etc.

All of them take either a `Particle`, a `PDGID` or a simple `Integer` (or
anything which can be converted to an `Integer`) which represents a PDG ID, as
input, so that the general API is the following:

```julia
helperfunction(p::Union{Particle, PDGID, Integer})
```

Here is a list of the currently available helper functions:

- `hasup(p)`
- `hasdown(p)`
- `hasstange(p)`
- `hascharm(p)`
- `hasbottom(p)`
- `hastop(p)`
- `isquark(p)`
- `isstandard(p)`
- `isfundamental(p)`
- `fundamentalid(p)`
- `islepton(p)`
- `ismeson(p)`
- `isbaryon(p)`
- `ishadron(p)`
- `isRhadron(p)`
- `isSUSY(p)`
- `ispentaquark(p)`
- `isgaugebosonorhiggs(p)`
- `issmgaugebosonorhiggs(p)`
- `istechnicolor(p)`
- `iscompositequarkorlepton(p)`
- `isdyon(p)`
- `isdiquark(p)`
- `isgeneratorspecific(p)`
- `isspecial(p)`
- `isQball(p)`
- `hasfundamentalanti(p)`
- `isnucleus(p)`
- `A(p)`
- `Z(p)`
