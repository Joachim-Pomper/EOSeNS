struct BarotropicEOS <: EquationOfState 
    name::String
    data::Vector{Vector{Float64}}
    labels::Dict{String, Int64}
    symbols::Vector{String}
    dimensions::Vector{UniCon.Dimension}
    unit_system::UnitSystem 
end

function Base.string(eos::BarotropicEOS)
    str = "BarotropicEOS: $(eos.name)"*"\n"
    n = length(eos.data[1])
    str = str * "n_data: $(n) \n"
    len_str = 0
    len_unit = 0
    for label in getLabels(eos)
        len_str = max(length(label), len_str)
        len_unit = max(length(string(getUnit(eos, label))), len_unit)
    end
    hline = "-"^(len_str+2)*"|"*"-"^21*"|"*"-"^(1+len_unit)*"\n"  
    str = str*hline
    str = str*"Lables"*" "^(len_str - 4)*":  min_val - max_val "
    str = str*" :  Units \n" 
    str = str*hline
    for label in getLabels(eos)
        data = eos[label]
        min_val = minimum(data)
        max_val = maximum(data)
        str = str * label * " "^(len_str-length(label))
        str = str * "  :  $(@sprintf("%7.4g", min_val)) - $(@sprintf("%-7.4g", max_val))"
        str = str * "  : "*string(getUnit(eos, label))*"\n"
    end 
    
    return str
end

function Base.show(io::IO, ::MIME"text/plain", eos::BarotropicEOS)
     print(io, string(eos))
end

function getLabels(eos::BarotropicEOS)
    # ToDo: Sort them in the right order 
    return collect(keys(eos.labels))
end

function getDimension(eos::BarotropicEOS, i::String)
    return eos.dimensions[eos.labels[i]]
end

function getUnit(eos::BarotropicEOS, i::String)
    return (getDimension(eos, i), eos.unit_system)
end

function getSymbol(eos::BarotropicEOS, i::String)
    return eos.symbols[eos.labels[i]]
end

function isaLabel(eos::BarotropicEOS, label::String)
    if label in keys(eos.labels)
        return true
    else
        return false
    end
end

function renameLabel!(eos::BarotropicEOS, old_name::String, new_name::String)

    if isaLabel(eos::BarotropicEOS, old_name)
        idx = eos.labels[old_name]
        eos.labels[new_name] = idx 
        delete!(eos.labels, old_name)
    else
        throw(KeyError(old_name))
        #Todo: Make custom error
    end
    return eos
end

function Base.getindex(eos::BarotropicEOS, i::String)

    try
        return eos.data[eos.labels[i]]
    
    catch error
        if isa(error, KeyError)
            throw(KeyError(i))
            #Todo: Make custom error
        else 
            throw(error)
        end
    end
    
end

function interpolated(eos::BarotropicEOS, y_name::String, x_name::String)
    
    x_data = eos[x_name]
    y_data = eos[y_name]
    return QuadraticInterpolation(y_data, x_data, :Backward)
end

function writeHDF5BarotropicEOS(eos::BarotropicEOS, file_name)

    eos_name = eos.name
    n_label = length(eos.data)
    h5open(file_name, "w") do file
        create_group(file, eos_name)
        g = file[eos_name]
        base_unit_names, base_unit_si_values = docRep(eos.unit_system)
        attrs(g)["base_unit_names"] = base_unit_names
        attrs(g)["base_unit_si_values"] = base_unit_si_values
        
        # Write datasets
        for label in getLabels(eos)
            g[label] = eos[label]
            dataset = g[label]
            attrs(dataset)["symbol"] = getSymbol(eos, label)
            attrs(dataset)["dimension"] = docRep(getDimension(eos, label))
        end
    end
end

function readHDF5BarotropicEOS(file_name, eos_name)

    data = Vector{Vector{Float64}}(undef,0)
    labels = Dict{String,Int64}()
    dimensions = Vector{UniCon.Dimension}(undef,0)
    symbols = Vector{String}(undef,0)
    
    h5open(file_name, "r") do file
        group_names = keys(file)
        if eos_name in group_names
            group = file[eos_name]

            #reconstruct unit system
            base_unit_names = attrs(group)["base_unit_names"]
            base_unit_si_values = attrs(group)["base_unit_si_values"]
            unit_system = UnitSystem(collect(zip(base_unit_names,base_unit_si_values))...)
        
            #read the data
            idx_label = 0
            for label in keys(group)
                idx_label +=1
                labels[label] = idx_label
                dataset = group[label]
                append!(data, [read(dataset)])
                append!(dimensions, [UniCon.Dimension(attrs(dataset)["dimension"]...)])
                append!(symbols, [attrs(dataset)["symbol"]])
            end

            return BarotropicEOS(eos_name, data, labels, symbols, dimensions, unit_system)
            
        else
            throw(DomainError(eos_name, "`"*eos_name*"` is not a group name in the HDF5 file."))
        end
    end
end

function convertUnits(eos::BarotropicEOS, new_usy::UnitSystem)
    
    new_data = Vector{Vector{Float64}}(undef,0)
    labels = eos.labels
    dimensions = eos.dimensions
    old_usy = eos.unit_system

    for idx in 1:length(eos.data)
        new_values = convertValue(eos.data[idx], eos.dimensions[idx], old_usy, new_usy)
        append!(new_data, [new_values])
    end
    return BarotropicEOS(eos.name, new_data, eos.labels, eos.symbols, eos.dimensions, new_usy)
end


