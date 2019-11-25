module Corpuscles

using DelimitedFiles
using Unitful
using UnitfulAtomic
using Printf

import Base: show, convert, showerror

export Particle, PDGID, PythiaID, GeantID

# Julia 1.0 compatibility
eachrow_(x) = (x[i, :] for i in 1:size(x)[1])
eachcol_(x) = (x[:, j] for j in 1:size(x)[2])

const _data_dir = abspath(joinpath(@__DIR__, "..", "data"))

function Base.parse(::Type{Rational{T}}, val::AbstractString) where {T <: Integer}
    !('/' in val) && return parse(T, val) // 1
    nums, denoms = split(val, '/', keepempty=false)
    num = parse(T, nums)
    denom = parse(T, denoms)
    return num//denom
end

abstract type ParticleID end
struct PDGID <: ParticleID
    value
end
struct GeantID <: ParticleID
    value
end
struct PythiaID <: ParticleID
    value
end

@enum PDGStatus begin
    Common      = 0
    Rare        = 1
    Unsure      = 2
    Further     = 3
    Nonexistent = 4
end

@enum InvProperty begin
    Same = 0
    Barred = 1
    ChargeInv = 2
end

struct MeasuredValue{D}
    value::Quantity{T1,D,U1} where {T1 <: Real, U1 <: Unitful.Units}
    lower_limit::Quantity{T2,D,U2} where {T2 <: Real, U2 <: Unitful.Units}
    upper_limit::Quantity{T3,D,U3} where {T3 <: Real, U3 <: Unitful.Units}
end

const _energy_dim = Unitful.dimension(u"J")
const _charge_dim = Unitful.dimension(u"C")

struct Particle
    pdgid::PDGID
    mass::MeasuredValue{_energy_dim}
    width::Union{Missing, MeasuredValue{_energy_dim}}
    charge::Quantity{T,_charge_dim,U} where {T<:Real, U<: Unitful.Units}
    isospin::Union{Missing, Rational{Int8}}
    parity::Union{Missing, Int8}
    gparity::Union{Missing, Int8}
    cparity::Union{Missing, Int8}
    antiprop::InvProperty
    rank::Int8
    status::PDGStatus
    name::String
    quarks::String
    latex::String
end

function read_conversion_csv(filepath::AbstractString)
    file_content = readdlm(filepath, ',', AbstractString, skipstart=2, comments=true)
    conversions = parse.(Int, file_content[:,1:3])
end

const _id_conversion_cols = Dict(PDGID => 1, GeantID => 3, PythiaID => 2)
const _id_conversion_tbl = read_conversion_csv(joinpath(_data_dir, "conversions.csv"))

Particle(id::ParticleID) = _current_particle_dct[convert(PDGID, id)]
Particle(id::Integer) = Particle(PDGID(id))

struct IDException <: Exception 
    var::AbstractString
end

Base.showerror(io::IO, e::IDException) = Printf.@printf(io, "ParticleID Error: %s", e.var)

function Base.convert(t::Type{X}, id::Y) where {X <: ParticleID, Y <: ParticleID}
    if isequal(X, Y)
        return id
    end
    val_col = _id_conversion_cols[t]
    key_col = _id_conversion_cols[Y]
    row = findfirst(x->isequal(x, id.value), _id_conversion_tbl[:,key_col])
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

const ParticleDict = Dict{PDGID, Particle}

function read_particle_csv(filepath::AbstractString)
    file_content = readdlm(filepath, ',', AbstractString)
    header = string.(file_content[1,:])
    dct_particles = ParticleDict()
    for row in eachrow_(file_content[2:end,:])
        pdgid       = PDGID(parse(Int64, row[1]))
        mass_value  = parse(Float64, row[2]) * u"MeV"
        mass_lower  = parse(Float64, row[3]) * u"MeV"
        mass_upper  = parse(Float64, row[4]) * u"MeV"
        mass = MeasuredValue{_energy_dim}(mass_value, mass_lower, mass_upper)
        width_value  = parse(Float64, row[5]) * u"MeV"
        width_lower  = parse(Float64, row[6]) * u"MeV"
        width_upper  = parse(Float64, row[7]) * u"MeV"
        width = MeasuredValue{}(width_value, width_lower, width_upper)
        isospin = if (row[8] in ["", "?"]) missing else parse(Rational{Int8}, row[8]) end
        gparity = read_parity(row[9])
        parity = read_parity(row[10])
        cparity = read_parity(row[11])
        antiprop = InvProperty(parse(Int8, row[12]))
        charge = parse(Int8, row[13]) // 3 * u"e_au"
        rank = parse(Int8, row[14])
        status = PDGStatus(parse(Int8, row[15]))
        name = row[16]
        quarks = row[17]
        latex = row[18]
        dct_particles[pdgid] = Particle(pdgid,
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
                                        latex)
    end
    dct_particles
end

"""
    available_catalog_files()

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
    filter!(s->occursin(".csv",s), dir_content)
    filter!(s->!occursin("conversions", s), dir_content)
    joinpath.(_data_dir, dir_content)
end

const _catalogs = available_catalog_files()

const _default_year = "2019"
const _default_catalog = filter(s->occursin(_default_year,s), _catalogs)[end]

const _current_particle_dct = read_particle_csv(_default_catalog)

"""
    use_catalog_file(filepath::AbstractString)

This function reads a given catalog file and sets it as reference

# Arguments
- `filepath::AbstractString`: filepath to the catalog file

# Examples
```julia-repl
julia> Corpuscles.use_catalog_file("/home/foobar/dev/Corpuscles.jl/data/particle2019.csv")
```
"""
function use_catalog_file(filepath::AbstractString)
    _current_particle_dct = read_particle_csv(filepath)
    return
end

function show(io::IO, m::MeasuredValue)
    if isapprox(m.upper_limit, m.lower_limit)
        print(io, "$(m.value) ± $(m.lower_limit)")
    else
        print(io, "$(m.value) + $(m.lower_limit) - $(m.lower_limit)")
    end
    return
end

function show(io::IO, p::Particle)
    Printf.@printf(io, "\n%-8s %-12s", "Name:", p.name)
    Printf.@printf(io, "%-7s %-10s", "PDGid:", p.pdgid)
    Printf.@printf(io, " %-7s %s", "LaTex:", "\$$(p.latex)\$\n\n")
    Printf.@printf(io, "%-8s %s\n", "Status:", p.status)
    println(io, "\nParameters:")
    println(io, "-----------")
    fields = Dict("Mass" => p.mass,
                  "Width" => p.width,
                  "Q (charge)" => p.charge,
                  "C (charge parity)" => p.cparity,
                  "P (space parity)" => p.parity,
                  "G (G-parity)" => p.gparity,
                  "Isospin" => p.isospin,
                  "Composition" => p.quarks)
    for (key, value) in fields
        if value isa MeasuredValue || !ismissing(value) && !isempty(value)
            Printf.@printf(io, "%-19s = %s\n",key, value)
        end
    end
end

end # module
