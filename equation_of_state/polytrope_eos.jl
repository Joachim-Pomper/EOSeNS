### Polytropic EOS as constructor 
"""
Create a BarotropicEOS struct from a specified polytrop

INPUT:
    gamma   [-]
    k       [MeV*(fm)^(3(gamma-1))]
    mb      []
    n       [(fm)^(-3)] 
"""
function polytropic_eos(n::Vector{Float64}, k::Float64, gamma::Float64, mb::Float64; name = "polytrop_eos") 

    if gamma == 1.
        throw(ArgumentError("Argument `gamma` cannot be 1."))
    end

    p = k.*n.^gamma               
    e = mb.*n .+ k.*n.^gamma ./ (gamma-1)
    de_dp = (p ./ k).^(1 ./ gamma-1) ./ (gamma * k) .+ 1.0 ./ (gamma - 1.0)
    d2e_dp2 = (p ./ k).^(1 ./ gamma - 2) .* (1.0/gamma - 1) ./ (gamma * k^2)
    d3e_dp3 = (p ./ k).^(1 ./ gamma - 3) .* (1.0/gamma - 1.0) * (1.0/gamma - 2.0) / (gamma*k^3)

    data = [
        n, 
        e, 
        p, 
        de_dp,
        d2e_dp2,
        d3e_dp3
    ]
    symbols = [
        "n",
        "e", 
        "p" , 
        "de_dp", 
        "d2e_dp3", 
        "d3e_dp3"
    ]
    dimensions = [
        NumberDensityDim(),
        EnergyDensityDim(),
        PressureDim(),
        EnergyDensityDim() / PressureDim(),
        EnergyDensityDim() / PressureDim() / PressureDim(),
        EnergyDensityDim() / PressureDim() / PressureDim() 
    ] 
    unit_system = NuclearUnits()

    lables = Dict(
        "number_density" => 1,
        "energy_density" => 2,
        "pressure" => 3,
        "de_dp" => 4,
        "d2e_dp3" => 5,
        "d3e_dp3" => 6
    )

    return BarotropicEOS(name, data, lables, symbols, dimensions, unit_system)

end
