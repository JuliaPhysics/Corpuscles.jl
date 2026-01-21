module Corpuscles

using DocStringExtensions
using DelimitedFiles
using Unitful
using UnitfulAtomic
using Printf

import Base

export Particle, PDGID, PythiaID, Geant3ID, particles

# helpers.jl
export isfundamental, isstandard
export isquark, islepton, ismeson, isbaryon, ishadron
export isRhadron, isSUSY, ispentaquark, isdyon, isnucleus, isdiquark
export istechnicolor, iscompositequarkorlepton
export isgaugebosonorhiggs, issmgaugebosonorhiggs
export isgeneratorspecific, isspecial, isQball, hasfundamentalanti
export hasdown, hasup, hascharm, hasstrange, hasbottom, hastop
export A, Z, charge, threecharge, J, S, L, jspin, sspin, lspin
export value_GeV, value_MeV, uncertainty_GeV, uncertainty_MeV

# Julia 1.0 compatibility
eachrow_(x) = (x[i, :] for i = 1:size(x)[1])
isnothing(::Any) = false
isnothing(::Nothing) = true


const _data_dir = abspath(joinpath(@__DIR__, "..", "data"))
include("pidNames.jl")
include("similarity.jl")

function parserational(::Type{Rational{T}}, val::AbstractString) where {T<:Integer}
    !('/' in val) && return parse(T, val) // 1
    nums, denoms = split(val, '/', keepempty = false)
    num = parse(T, nums)
    denom = parse(T, denoms)
    return num//denom
end

abstract type ParticleID end

"""
$(SIGNATURES)

PDG IDs consist of 7 digits prefixed by a sign, following the scheme:

    +/- N Nr Nl Nq1 Nq2 Nq3 Nj

Those are accessible as fields. There are PDG IDs with more than 7 digits for
non-standard particles such as Q-balls. To support those, we follow the
implementation in the SciKit-HEP `particle` Python package, which
introduced N8, N9 and N10.

"""
struct PDGID <: ParticleID
    value::Int32
    Nj::Int8
    Nq3::Int8
    Nq2::Int8
    Nq1::Int8
    Nl::Int8
    Nr::Int8
    N::Int8
    N8::Int8
    N9::Int8
    N10::Int8

    function PDGID(value)
        d = digits(abs(value), pad = 10)
        # splatting in `new()` is only supported in Julia 1.2+
        return new(value, d[1], d[2], d[3], d[4], d[5], d[6], d[7], d[8], d[9], d[10])
    end
end

struct Geant3ID <: ParticleID
    value::Any
end

@deprecate GeantID Geant3ID true

struct PythiaID <: ParticleID
    value::Any
end

@enum PDGStatus begin
    Common = 0
    Rare = 1
    Unsure = 2
    Further = 3
    Nonexistent = 4
end

@enum InvProperty begin
    Same = 0
    Barred = 1
    ChargeInv = 2
end

struct MeasuredValue{D}
    value::Quantity{T1,D,U1} where {T1<:Real,U1<:Unitful.Units}
    lower_limit::Quantity{T2,D,U2} where {T2<:Real,U2<:Unitful.Units}
    upper_limit::Quantity{T3,D,U3} where {T3<:Real,U3<:Unitful.Units}
end

function Base.isless(x::MeasuredValue, y::Quantity)
    x.value + x.upper_limit < y
end

function Base.isless(x::Quantity, y::MeasuredValue)
    x < y.value - y.lower_limit
end

function Base.isless(x::MeasuredValue, y::MeasuredValue)
    x.value + x.upper_limit < y.value - y.lower_limit
end

function Base.isapprox(x::Quantity, y::MeasuredValue)
    (x < y.value + y.upper_limit) & (x > y.value - y.lower_limit)
end

function Base.isapprox(y::MeasuredValue, x::Quantity)
    (x < y.value + y.upper_limit) & (x > y.value - y.lower_limit)
end

function Base.isapprox(x::MeasuredValue, y::MeasuredValue)
    # TODO: I think this can be done way more elegant, but functionality is now given
    retval = isapprox(x, y.value + y.upper_limit)
    retval |= isapprox(x, y.value - y.lower_limit)
    retval |= isapprox(x, y.value)
    retval |= isapprox(y, x.value + x.upper_limit)
    retval |= isapprox(y, x.value - x.upper_limit)
    retval |= isapprox(y, x.value)
    retval
end

const _energy_dim = Unitful.dimension(u"J")
const _mass_dim = Unitful.dimension(u"kg")
const _charge_dim = Unitful.dimension(u"C")

struct Particle
    pdgid::PDGID
    mass::MeasuredValue{_mass_dim}
    width::Union{Missing,MeasuredValue{_mass_dim}}
    charge::Quantity{T,_charge_dim,U} where {T<:Real,U<:Unitful.Units}
    isospin::Union{Missing,Rational{Int8}}
    parity::Union{Missing,Int8}
    gparity::Union{Missing,Int8}
    cparity::Union{Missing,Int8}
    antiprop::InvProperty
    rank::Int8
    status::PDGStatus
    name::String
    quarks::String
    latex::String
    glyph::String
end

pdgid(p::PDGID) = p
pdgid(p::Particle) = p.pdgid
pdgid(x) = PDGID(Integer(x))

import Base: -
# anti-particle
function (-)(p::Particle)
    try
        Particle(-p.pdgid.value)
    catch # anti-particle is itself
        p
    end
end

function read_conversion_csv(filepath::AbstractString)
    file_content = readdlm(filepath, ',', AbstractString, skipstart = 2, comments = true)
    conversions = parse.(Int, file_content[:, 1:3])
end

const _id_conversion_cols = Dict(PDGID => 1, Geant3ID => 3, PythiaID => 2)
const _id_conversion_tbl = read_conversion_csv(joinpath(_data_dir, "conversions.csv"))

Particle(id::ParticleID) = catalog.particle_dict[Base.convert(PDGID, id)]
Particle(id::Integer) = Particle(PDGID(id))

struct IDException <: Exception
    var::AbstractString
end

Base.showerror(io::IO, e::IDException) = Printf.@printf(io, "ParticleID Error: %s", e.var)

function Base.convert(t::Type{X}, id::Y) where {X<:ParticleID,Y<:ParticleID}
    if isequal(X, Y)
        return id
    end
    val_col = _id_conversion_cols[t]
    key_col = _id_conversion_cols[Y]
    row = findfirst(x->isequal(x, id.value), _id_conversion_tbl[:, key_col])
    if iszero(id.value) || isequal(row, nothing)
        throw(IDException("No corresponding $t for $id found!"))
    end
    X(_id_conversion_tbl[row, val_col])
end


function read_parity(val::AbstractString)
    tmp = parse(Int8, val)
    if tmp == 5
        return missing
    else
        return tmp
    end
end

const ParticleDict = Dict{PDGID,Particle}

function read_particle_csv(filepath::AbstractString)
    file_content = readdlm(filepath, ',', AbstractString, comments = true)
    header = string.(file_content[1, :])
    dct_particles = ParticleDict()
    for row in eachrow_(file_content[2:end, :])
        pdgid = PDGID(parse(Int64, row[1]))
        mass_value = parse(Float64, row[2]) * u"MeV/c^2"
        mass_lower = parse(Float64, row[3]) * u"MeV/c^2"
        mass_upper = parse(Float64, row[4]) * u"MeV/c^2"
        mass = MeasuredValue{_mass_dim}(mass_value, mass_lower, mass_upper)
        width_value = parse(Float64, row[5]) * u"MeV/c^2"
        width_lower = parse(Float64, row[6]) * u"MeV/c^2"
        width_upper = parse(Float64, row[7]) * u"MeV/c^2"
        width = MeasuredValue{}(width_value, width_lower, width_upper)
        isospin = if (row[8] in ["", "?"])
            missing
        else
            parserational(Rational{Int8}, row[8])
        end
        gparity = read_parity(row[9])
        parity = read_parity(row[10])
        cparity = read_parity(row[11])
        antiprop = InvProperty(parse(Int8, row[12]))
        charge = parse(Int8, row[13]) // 3 * u"e_au"
        rank = parse(Int8, row[14])
        status = PDGStatus(parse(Int8, row[15]))
        name = _generatename(row[16], pdgid.value, antiprop, charge)
        glyph = get(PIDNames, pdgid.value, row[16])
        quarks = row[17]
        latex = row[18]
        dct_particles[pdgid] = Particle(
            pdgid,
            mass,
            width,
            charge,
            isospin,
            parity,
            gparity,
            cparity,
            antiprop,
            rank,
            status,
            name,
            quarks,
            latex,
            glyph,
        )
    end
    dct_particles
end

"""
$(SIGNATURES)

Generate a name for the particle, adding tilde if anti and -/+ according to the
charge.
"""
function _generatename(name, pdgid, antiprop, charge)
    prefix_tilde = antiprop == Barred && pdgid < 0 ? "~" : ""
    suffix_charge = ""
    if isinteger(charge.val)
        chrg = floor(Int, charge.val)
        if chrg == 0
            suffix_charge = "0"
        elseif chrg > 0
            suffix_charge = '+'^abs(chrg)
        else
            suffix_charge = '-'^abs(chrg)
        end
    end
    return prefix_tilde * name * suffix_charge
end



"""
$(SIGNATURES)

Function to get the available catalog files which are available within
the package and returns a list with the absolute filepaths.

# Examples
```julia-repl
julia> Corpuscles.available_catalog_files()
["/home/foobar/dev/Corpuscles.jl/data/particle2019.csv"]
```
"""
function available_catalog_files()
    dir_content = readdir(_data_dir)
    filter!(s->occursin(".csv", s), dir_content)
    filter!(s->!occursin("conversions", s), dir_content)
    joinpath.(_data_dir, dir_content)
end

const _catalogs = available_catalog_files()

const _default_year = "2025"
const _default_catalog = filter(s->occursin(_default_year, s), _catalogs)[end]

mutable struct Catalog
    particle_dict::ParticleDict
    particles::Vector{Particle}

    Catalog(d::ParticleDict) = new(d, collect(values(d)))
end

"""
$(SIGNATURES)

Returns the full list of particles from the currently selected catalog.
"""
particles() = catalog.particles

const catalog = Catalog(read_particle_csv(_default_catalog))

# inverse catelog for `Particle(name)` construction
const inv_catalog =
    Dict{String,PDGID}(p.name => p.pdgid for p in catalog.particles if p.pdgid.value>0)

"""
$(SIGNATURES)

This function reads a given catalog file and sets it as reference

# Arguments
- `filepath::AbstractString`: filepath to the catalog file

# Examples
```julia-repl
julia> Corpuscles.use_catalog_file("/home/foobar/dev/Corpuscles.jl/data/particle2019.csv")
```
"""
function use_catalog_file(filepath::AbstractString)
    # TODO: there is surely a better way to manage this "closure"-like thing
    catalog.particle_dict = read_particle_csv(filepath)
    catalog.particles = collect(values(catalog.particle_dict))
    return
end

function Base.show(io::IO, m::MeasuredValue)
    if isapprox(m.upper_limit, m.lower_limit)
        print(io, "$(m.value) ± $(m.lower_limit)")
    else
        print(io, "$(m.value) + $(m.upper_limit) - $(m.lower_limit)")
    end
    return
end

function Base.show(io::IO, p::PDGID)
    Printf.@printf(io, "PDGID(%s)", p.value)
end

function Base.show(io::IO, p::Particle)
    Printf.@printf(io, "Particle(%s) %s", p.pdgid.value, p.name)
end

function Base.print(io::IO, p::Particle)
    Printf.@printf(io, "%-8s %s\n", "Name:", p.name)
    Printf.@printf(io, "%-8s %s\n", "PDG ID:", p.pdgid.value)
    Printf.@printf(io, "%-8s %s\n", "LaTeX:", "\$$(p.latex)\$")
    Printf.@printf(io, "%-8s %s\n", "Status:", p.status)
    fields = Dict(
        "Mass" => p.mass,
        "Width" => p.width,
        "Q (charge)" => p.charge,
        "C (charge parity)" => p.cparity,
        "P (space parity)" => p.parity,
        "G (G-parity)" => p.gparity,
        "Isospin" => p.isospin,
        "Composition" => p.quarks,
    )
    for (key, value) in fields
        if value isa MeasuredValue || !ismissing(value) && !isempty(value)
            Printf.@printf(io, "%s = %s\n", key, value)
        end
    end
end


include("helpers.jl")
include("literals.jl")

end # module
