const COMPOSE_DATA_PATH = joinpath(@__DIR__, "..", "..", "data", "CompOSE")

# -------------------------------------------------------------------------- #
# Compose custom errors                                                      #  
# -------------------------------------------------------------------------- #

abstract type CompOSError <: Exception end
abstract type MissingCompOSEFileError <: CompOSError end

struct MissingCompOSEThermoFileError <: MissingCompOSEFileError
    file_name::String
end

function Base.showerror(io::IO, error::MissingCompOSEThermoFileError)
    print(io, "Missing CompOSE thermodynamic file: $(error.file_name)")
end


struct MissingCompOSETFileError <: MissingCompOSEFileError
    file_name::String
end

function Base.showerror(io::IO, error::MissingCompOSETFileError)
    print(io, "Missing CompOSE temperature grid file required by varying thermo column 1: $(error.file_name)")
end


struct MissingCompOSENBFileError <: MissingCompOSEFileError
    file_name::String
end

function Base.showerror(io::IO, error::MissingCompOSENBFileError)
    print(io, "Missing CompOSE baryon-density grid file required by varying thermo column 2: $(error.file_name)")
end


struct MissingCompOSEYQFileError <: MissingCompOSEFileError
    file_name::String
end

function Base.showerror(io::IO, error::MissingCompOSEYQFileError)
    print(io, "Missing CompOSE charge-fraction grid file required by varying thermo column 3: $(error.file_name)")
end


struct UnsupportedCompOSETableError <: CompOSError
    path::String
    varying::Tuple{Bool, Bool, Bool}
end

function Base.showerror(io::IO, error::UnsupportedCompOSETableError)
    print(io, "Unsupported CompOSE table at $(error.path): expected a cold NS EOS 1D table with varying columns (false, true, false), got $(error.varying)")
end

# -------------------------------------------------------------------------- #
# Compose consistency checks                                                 #  
# -------------------------------------------------------------------------- #

function resolveEOSDataPath(path_or_name::AbstractString)
    ispath(path_or_name) && return path_or_name

    compose_path = joinpath(COMPOSE_DATA_PATH, path_or_name)
    ispath(compose_path) && return compose_path

    throw(ArgumentError("EOS data path not found: $(path_or_name)"))
end

function getComposeTableDim(path::AbstractString; eos_name::AbstractString="eos")
    thermo_file = joinpath(path, eos_name * ".thermo")
    isfile(thermo_file) || throw(MissingCompOSEThermoFileError(thermo_file))

    function _firstThreeComposeValues(line::AbstractString)
        tokens = split(strip(split(line, '#'; limit=2)[1]); limit=4)
        length(tokens) < 3 && return nothing

        value1 = tryparse(Float64, tokens[1])
        value2 = tryparse(Float64, tokens[2])
        value3 = tryparse(Float64, tokens[3])
        (value1 === nothing || value2 === nothing || value3 === nothing) && return nothing

        return (value1, value2, value3)
    end

    reference = nothing
    varying = (false, false, false)
    open(thermo_file, "r") do io
        eof(io) && throw(ArgumentError("Empty CompOSE thermo file: $(thermo_file)"))
        readline(io)

        while !eof(io)
            values = _firstThreeComposeValues(readline(io))
            values === nothing && continue

            if reference === nothing
                reference = values
                continue
            end

            varying = (
                varying[1] | (values[1] != reference[1]),
                varying[2] | (values[2] != reference[2]),
                varying[3] | (values[3] != reference[3])
            )

            all(varying) && break
        end
    end

    reference === nothing && throw(ArgumentError("No thermo data rows found in CompOSE file: $(thermo_file)"))

    return count(identity, varying), varying
end

function checkCompOSEFileConsistency(path::AbstractString; eos_name::AbstractString="eos")
    _, varying = getComposeTableDim(path; eos_name=eos_name)

    t_file = joinpath(path, eos_name * ".t")
    varying[1] && !isfile(t_file) && throw(MissingCompOSETFileError(t_file))

    nb_file = joinpath(path, eos_name * ".nb")
    varying[2] && !isfile(nb_file) && throw(MissingCompOSENBFileError(nb_file))

    yq_file = joinpath(path, eos_name * ".yq")
    varying[3] && !isfile(yq_file) && throw(MissingCompOSEYQFileError(yq_file))

    return true
end

# -------------------------------------------------------------------------- #
# load EOS                                                                   #  
# -------------------------------------------------------------------------- #

function readComposeColdNS1dEos(path_or_name::AbstractString, os_name::AbstractString="eos")
    path = resolveEOSDataPath(path_or_name)
    _, varying = getComposeTableDim(path; eos_name=os_name)
    varying == (false, true, false) || throw(UnsupportedCompOSETableError(path, varying))
    checkCompOSEFileConsistency(path; eos_name=os_name)

    nb_file = joinpath(path, os_name * ".nb")
    thermo_file = joinpath(path, os_name * ".thermo")

    # Load grid files
    nb_values = Float64[]
    for line in eachline(nb_file)
        stripped = strip(split(line, '#'; limit=2)[1])
        isempty(stripped) && continue
        push!(nb_values, parse(Float64, stripped))
    end
    number_density_grid = nb_values[3:end] #remove first two values. 
                                           # See README in data/CompOSE/README.md
    # Load thermo quantities
    nb_indices = Int[]
    pressure_over_nb = Float64[]
    energy_per_mass  = Float64[]
    neutron_mass     = 0.0
    open(thermo_file, "r") do io
        header = split(strip(split(readline(io), '#'; limit=2)[1])) # See README in data/CompOSE/README.md
        neutron_mass = parse(Float64, header[1])                    # Explains which entries are neutron mass

        while !eof(io)
            tokens = split(strip(split(readline(io), '#'; limit=2)[1]))
            isempty(tokens) && continue
            push!(nb_indices, parse(Int, tokens[2]))
            push!(pressure_over_nb, parse(Float64, tokens[4])) # See README in data/CompOSE/README.md
            push!(energy_per_mass, parse(Float64, tokens[10])) # Explains which rows have which quantity
        end
    end

    if minimum(nb_indices) == 0
        @inbounds for i in eachindex(nb_indices)
            nb_indices[i] += 1
        end
    end

    n = length(nb_indices)
    number_density = Vector{Float64}(undef, n)
    pressure = Vector{Float64}(undef, n)
    energy_density = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        nb = number_density_grid[nb_indices[i]]
        number_density[i] = nb
        pressure[i] = pressure_over_nb[i] * nb
        energy_density[i] = (energy_per_mass[i] + 1.0) * nb * neutron_mass
    end

    order = sortperm(number_density)
    data = [number_density[order], energy_density[order], pressure[order]]
    labels = Dict(
        "number_density" => 1,
        "energy_density" => 2,
        "pressure" => 3
    )
    symbols = ["n", "e", "p"]
    dimensions = [NumberDensityDim(), EnergyDensityDim(), PressureDim()]

    return BarotropicEOS(basename(normpath(path)), data, labels, symbols, dimensions, NuclearUnits())
end
