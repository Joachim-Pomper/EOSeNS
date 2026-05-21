const BBP_TABLE_PATH = joinpath(@__DIR__, "..", "..", "data", "crust_eos_data", "bbp_table.txt")

function loadBBPTable()
    return Matrix{Float64}(readdlm(BBP_TABLE_PATH; comments=true, comment_char='#'))
end

function bbpnvCrustEOS()
    
    c_light = 299792458.0*100.0 # cm/s
    bbp_table = loadBBPTable()
    data = [bbp_table[:,1]*c_light^2,bbp_table[:,2],bbp_table[:,3],bbp_table[:,4],bbp_table[:,5],bbp_table[:,6]]

    labels = Dict(
        "energy_density"       => 1,
        "pressure"             => 2,
        "number_density"       => 3,
        "number_of_protons"    => 4,
        "nucleus_mass_number"  => 5,
        "adiabatic_index"      => 6,
    )

    symbols = ["ρ", "P", "n_b", "Z", "A", "𝛾" ]

    name = "Baym-Bethte-Pethick-Negele-Vautherin EOS"

    dimensions = [
        EnergyDensityDim(),
        PressureDim(),
        NumberDensityDim(),
        DimensionLessDim(),
        DimensionLessDim(),
        DimensionLessDim()
    ]

    unit_system = CGS()

    return BarotropicEOS(
        name, 
        data, 
        labels, 
        symbols,
        dimensions,
        unit_system
    )
end

function addCrustByJumping(eos_inner::BarotropicEOS, n0 = 0.08)

    # get crust EOS
    eos_crust = convertUnits(bbpnvCrustEOS(), eos_inner.unit_system)
    
    # convert n0 into the correct units
    n0 = convertValue(n0, NumberDensityDim(), NuclearUnits(), eos_inner.unit_system)

    if eos_crust["number_density"][end] < n0

        n0_err = convertValue(n0, NumberDensityDim(), eos_inner.unit_system, NuclearUnits())
        bound = eos_crust["number_density"][end] 
        bound = convertValue(bound, NumberDensityDim(), eos_inner.unit_system, NuclearUnits())
        str = "Crust cannot be added to EOS. "
        str = str*"The number_density n0 of the mathching point is $(n0). "
        str = str*"The crust EOS requires the value to be above $(eos_crust["number_density"][end] )."
        throw(ErrorException(str))

    elseif eos_inner["number_density"][1] > n0

        n0_err = convertValue(n0, NumberDensityDim(), eos_inner.unit_system, NuclearUnits())
        bound = eos_inner["number_density"][1]
        bound = convertValue(bound, NumberDensityDim(), eos_inner.unit_system, NuclearUnits())
        str = "Crust cannot be added to EOS. "
        str = str*"The number_density n0 of the mathching point is $(n0) [1/fm^3]. "
        str = str*"The inner EOS requires the value to be below $(eos_inner["number_density"][1]) [1/fm^3]."
        throw(ErrorException(str))

    end

    # get data
    n_inner = eos_inner["number_density"] 
    n_curst = eos_crust["number_density"]
    p_inner = eos_inner["pressure"] 
    p_curst = eos_crust["pressure"]
    e_inner = eos_inner["energy_density"] 
    e_curst =eos_crust["energy_density"]

    # define matched interpolation
    idx_inner = n_inner.>= n0
    idx_crust = n_curst .< n0 

    # combine equations of state
    p_data = vcat(p_curst[idx_crust], p_inner[idx_inner])
    e_data = vcat(e_curst[idx_crust], e_inner[idx_inner])
    n_data = vcat(n_curst[idx_crust], n_inner[idx_inner])
    data = [p_data, e_data, n_data]
    labels = Dict(
        "pressure"=>1,
        "energy_density"=>2,
        "number_density"=>3
    )
    symbols = ["P", "e", "n_b"]
    dimensions = [PressureDim(), EnergyDensityDim(), NumberDensityDim()]
    unit_system = eos_inner.unit_system
    eos_name = eos_inner.name*" with BBPNV crust"

    return BarotropicEOS(
        eos_name, 
        data,
        labels,
        symbols, 
        dimensions, 
        unit_system
    )

    return
end

function addCrustByGlueing(eos_inner::BarotropicEOS)

    # get crust EOS
    eos_crust = convertUnits(bbpnvCrustEOS(), eos_inner.unit_system)

    # get interpolating functions 
    inner_interp = interpolated(eos_inner, "energy_density", "pressure")
    crust_interp = interpolated(eos_crust, "energy_density", "pressure")

    # get matching interval 
    lower_bound = eos_inner["pressure"][1]
    upper_bound = eos_crust["pressure"][end]

    if lower_bound > upper_bound
        str = "Crust cannot be added to EOS. "
        str = str*"The lowest value of the pressure is $(lower_bound). "
        str = str*"The crust EOS requires the pressure to be below $(upper_bound)."
        throw(ErrorException(str))
    end

    # find the matching point
    diff_interp(x) =  inner_interp(x) - crust_interp(x)
    matching_x = find_zero(diff_interp, (lower_bound, upper_bound)) 

    # define matched interpolation
    function match_interp(x::Float64)
        if x >= matching_x
            return inner_interp(x) 
        else 
            return crust_interp(x)
        end
    end
    
    # match the equations 
    x_values = sort(unique(vcat(eos_inner["pressure"], eos_crust["pressure"], [matching_x]))) 
    y_values = match_interp.(x_values)  

    data = [x_values, y_values]
    labels = Dict(
        "pressure"=>1,
        "energy_density"=>2
    )
    symbols = ["P", "e"]
    dimensions = [PressureDim(), EnergyDensityDim()]
    unit_system = eos_inner.unit_system
    eos_name = eos_inner.name*" with BBPNV crust"

    return BarotropicEOS(
        eos_name, 
        data,
        labels,
        symbols, 
        dimensions, 
        unit_system
    )
    
end
