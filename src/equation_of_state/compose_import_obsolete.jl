const COMPOSE_NEUTRON_MASS_MEV = 939.56542052



function _readNumericRows(file_name::AbstractString; min_columns::Int=1)
    rows = Vector{Vector{Float64}}()
    for line in eachline(file_name)
        stripped = strip(split(line, '#'; limit=2)[1])
        isempty(stripped) && continue

        values = Float64[]
        parse_failed = false
        for token in split(stripped)
            value = tryparse(Float64, token)
            if value === nothing
                parse_failed = true
                break
            end
            push!(values, value)
        end

        if !parse_failed && length(values) >= min_columns
            push!(rows, values)
        end
    end

    if isempty(rows)
        throw(ArgumentError("No numeric data found in `$(file_name)`."))
    end

    ncols = maximum(length, rows)
    data = fill(NaN, length(rows), ncols)
    for (i, row) in pairs(rows)
        data[i, 1:length(row)] .= row
    end
    return data
end

function _readComposeGridValues(file_name::AbstractString)
    values = _readNumericRows(file_name)

    if size(values, 2) == 1
        grid = vec(values)
        if length(grid) > 2 && isinteger(grid[2]) && Int(grid[2]) <= length(grid) - 2
            return grid[3:end]
        end
        if length(grid) > 1 && isinteger(grid[1]) && Int(grid[1]) == length(grid) - 1
            return grid[2:end]
        end
        return grid
    end

    if all(values[:, 1] .== 1:size(values, 1))
        return values[:, end]
    end

    if size(values, 1) > 1 && isinteger(values[1, 1]) && Int(values[1, 1]) == size(values, 1) - 1
        return values[2:end, end]
    end

    return values[:, end]
end

function _readComposeThermoData(file_name::AbstractString)
    header = _readNumericRows(file_name)[1, :]
    raw = _readNumericRows(file_name; min_columns=10)
    if size(raw, 1) < 1 || size(raw, 2) < 10
        throw(ArgumentError("`$(file_name)` does not look like a CompOSE eos.thermo file."))
    end
    return header, raw
end

"""
    readCompOSEBarotropicEOS(path; name=basename(path), t_index=0, yq_index=1)

Read an unpacked CompOSE EoS directory and return a `BarotropicEOS`.

This importer targets one-dimensional cold/barotropic slices represented by
`eos.nb` and `eos.thermo`. It filters `eos.thermo` by `t_index` and `yq_index`,
then reconstructs pressure and energy density in CompOSE nuclear units:
`n_B [fm^-3]`, `p [MeV/fm^3]`, and `e [MeV/fm^3]`.

INPUTS: 
    path     : 
    name     : Name of te EOS object generated. 
               Default value is created from path variable. 
    t_index  : 
    yq_index :
NOTE:
For standard CompOSE `eos.nb` files, the first two numeric entries are treated
as grid metadata and skipped before reading the density grid.
"""
function readCompOSEBarotropicEOS(
    path::AbstractString; 
    name::AbstractString=basename(normpath(path)), 
    t_index::Int=0, 
    yq_index::Int=1
    )
    
    nb_file = joinpath(path, "eos.nb")
    thermo_file = joinpath(path, "eos.thermo")

    if !isfile(nb_file)
        throw(ArgumentError("Missing CompOSE density grid file: $(nb_file)"))
    elseif !isfile(thermo_file)
        throw(ArgumentError("Missing CompOSE thermodynamic file: $(thermo_file)"))
    end

    number_density_grid   = _readComposeGridValues(nb_file)
    thermo_header, thermo = _readComposeThermoData(thermo_file)
    neutron_mass = thermo_header[1] > 0 ? thermo_header[1] : COMPOSE_NEUTRON_MASS_MEV

    idx = (thermo[:, 1] .== t_index) .& (thermo[:, 3] .== yq_index)
    thermo = thermo[idx, :]
    if isempty(thermo)
        throw(ArgumentError("No CompOSE thermo rows found for t_index=$(t_index), yq_index=$(yq_index)."))
    end

    order = sortperm(thermo[:, 2])
    thermo = thermo[order, :]

    nb_indices = Int.(thermo[:, 2])
    if minimum(nb_indices) == 0
        nb_indices .+= 1
    end
    if minimum(nb_indices) < 1 || maximum(nb_indices) > length(number_density_grid)
        throw(ArgumentError("CompOSE density indices are outside the eos.nb grid range."))
    end
    number_density = number_density_grid[nb_indices]
    pressure = thermo[:, 4] .* number_density
    energy_density = (thermo[:, 10] .+ 1.0) .* number_density .* neutron_mass

    data = [number_density, energy_density, pressure]
    labels = Dict(
        "number_density" => 1,
        "energy_density" => 2,
        "pressure" => 3
    )
    symbols = ["n", "e", "p"]
    dimensions = [NumberDensityDim(), EnergyDensityDim(), PressureDim()]

    return BarotropicEOS(name, data, labels, symbols, dimensions, NuclearUnits())
end
