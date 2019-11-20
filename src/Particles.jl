module Particles

using DelimitedFiles

# Julia 1.0 compatibility
eachrow_(x) = (x[i, :] for i in 1:size(x)[1])


function Base.parse(::Type{Rational{T}}, val::AbstractString) where {T <: Integer}
    !('/' in val) && return parse(T, val) // 1
    nums, denoms = split(val, '/', keepempty=false)
    num = parse(T, nums)
    denom = parse(T, denoms)
    return num//denom 
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

struct MeasuredValue{T}
    value::T
    lower_limit::T
    upper_limit::T
end

struct ParticleInfo
    pdgid::Int64
    mass::MeasuredValue{Float64}
    width::Union{Missing, MeasuredValue{Float64}}
    charge::Rational{Int8}
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

function read_parity(val::AbstractString)
    tmp = parse(Int8, val)
    if tmp == 5
        return missing
    else
        return tmp
    end
end

const ParticleDict = Dict{Int, ParticleInfo}

function read_particle_csv(filepath::AbstractString)
    file_content = readdlm(filepath, ',', AbstractString)
    header = string.(file_content[1,:])
    dct_particles = ParticleDict()
    for row in eachrow_(file_content[2:end,:])
        pdgid       = parse(Int64, row[1])
        mass_value  = parse(Float64, row[2])
        mass_lower  = parse(Float64, row[3])
        mass_upper  = parse(Float64, row[4])
        mass = MeasuredValue{Float64}(mass_value, mass_lower, mass_upper)
        width_value  = parse(Float64, row[5])
        width_lower  = parse(Float64, row[6])
        width_upper  = parse(Float64, row[7])
        width = MeasuredValue{Float64}(width_value, width_lower, width_upper)
        isospin = if (row[8] in ["", "?"]) missing else parse(Rational{Int16}, row[8]) end
        gparity = read_parity(row[9])
        parity = read_parity(row[10])
        cparity = read_parity(row[11])
        antiprop = InvProperty(parse(Int8, row[12]))
        charge = parse(Int8, row[13]) // 3
        rank = parse(Int8, row[14])
        status = PDGStatus(parse(Int8, row[15]))
        name = row[16]
        quarks = row[17]
        latex = row[18]
        dct_particles[pdgid] = ParticleInfo(pdgid,
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

const _data_dir = abspath(joinpath(@__DIR__, "..", "data"))

"""
    available_catalog_files()

Function to get the available catalog files which are available within 
the package and returns a list with the absolute filepaths.

# Examples
```julia-repl
julia> Particles.available_catalog_files()
["/home/foobar/dev/Particles.jl/data/particle2019.csv"]
```
"""
function available_catalog_files()
    dir_content = readdir(_data_dir)
    filter!(s->occursin(".csv",s), dir_content)
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
julia> Particles.use_catalog_file("/home/foobar/dev/Particles.jl/data/particle2019.csv")
```
"""
function use_catalog_file(filepath::AbstractString)
    _current_particle_dct = read_particle_csv(filepath)
    return
end

end # module
