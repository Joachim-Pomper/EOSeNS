const BBL_TABLE_DIR = joinpath(@__DIR__, "..", "..", "data", "bombaci_logoteta")
const BBL_TABLE_PATHS = Dict(
    1 => joinpath(BBL_TABLE_DIR, "bbl_delta_1.txt"),
    2 => joinpath(BBL_TABLE_DIR, "bbl_delta_2.txt")
)

function loadBBLTable(delta::Int64)
    if !haskey(BBL_TABLE_PATHS, delta)
        throw(DomainError(delta, "Input must be either 1 or 2"))
    end
    return Matrix{Float64}(readdlm(BBL_TABLE_PATHS[delta]; comments=true, comment_char='#'))
end

"""
See http://arxiv.org/abs/1805.11846
"""
function bblTableEOS(delta::Int64)
    
    if delta == 1
        name = "BBL1T EOS"
    elseif delta == 2
        name = "BBL2T EOS"
    else
        throw(DomainError(delta, "Input must be either 1 or 2"))
    end
    data = loadBBLTable(delta)
    data = [data[:,1], data[:,2], data[:,3], data[:,4], data[:,5]]


    labels = Dict(
        "number_density"    => 1,
        "energy_density"    => 2,  
        "pressure"          => 3, 
        "proton_fraction"   => 4,
        "electron_fraction" => 5
    )

    symbols = ["n", "ε", "P", "x_p", "x_e"]

    dimensions = [
        NumberDensityDim(),
        EnergyDensityDim(),
        PressureDim(),
        DimensionLessDim(),
        DimensionLessDim()
    ]

    unit_system = NuclearUnits()

    return BarotropicEOS(
        name, 
        data, 
        labels, 
        symbols,
        dimensions,
        unit_system
    )
end

### Bombaci Logoteta 

"""
See http://arxiv.org/abs/1805.11846
"""
function bblPolytropeEOS(delta::Int64; delta_step::Float64 = 0.001, n_dense_min::Float64 = 8e-2, n_dense_max::Float64=1.3)
    
    n = collect(n_dense_min:delta_step:n_dense_max)
    if delta == 1
        a = 945.199
        b = 293.551
        gamma = 2.82302
        name = "BBL1P EOS"
    elseif delta == 2
        a = 942.832
        b = 259.852
        gamma = 2.67041
        name = "BBL2P EOS"
    else
        throw(DomainError(delta, "Input must be either 1 or 2"))
    end

    mb = a
    k = (gamma-1)*b

    return polytropic_eos(n, k, gamma, mb, name = name)
end
