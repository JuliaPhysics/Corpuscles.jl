# Usage

The `Particle` struct can be used to create a particle. If an integer value is
passed, it will be interpreted as PDG ID, which is the primary particle
encoding in **Corpuscles.jl**:

```julia
julia> using Corpuscles

julia> p = Particle(-319)
Particle(-319)
```

To get an overview of the available particle information, use `print()`:

```julia
julia> print(p)
Name:    K(4)*(2045)
PDG ID:  -319
LaTeX:   $\bar{K}_{4}^{*}(2045)^{0}$
Status:  Common
Width = 198.0 MeV ± 30.0 MeV
Q (charge) = 0//1 e
Composition = Ds
Isospin = 1//2
Mass = 2045.0 MeV ± 9.0 MeV
P (space parity) = 1
```

The properties are accessible via attributes:

```julia
julia> fieldnames(Particle)
(:pdgid, :mass, :width, :charge, :isospin, :parity, :gparity, :cparity, :antiprop, :rank, :status, :name, :quarks, :latex)

julia> p.quarks
"Ds"

julia> p.isospin
1//2
```

For the mass, we use the [Unitful](https://github.com/PainterQubits/Unitful.jl)
package, which makes it easy to combine values with physical units:

```julia
julia> typeof(p.mass)
Corpuscles.MeasuredValue{𝐋^2 𝐌 𝐓^-2}

julia> p.mass
2045.0 MeV ± 9.0 MeV

julia> p.mass.value
2045.0 MeV

julia> p.mass.lower_limit
9.0 MeV

julia> p.mass.upper_limit
9.0 MeV
```

and also `Base.isless` and `Base.isapprox` are implemented so that the
lower and upper limits are taken into account, as seen here:

```julia
julia> using Unitful

julia> p.mass
2045.0 MeV ± 9.0 MeV

julia> p.mass > 2034u"MeV"
true
```
