module NEOS

    using HDF5
    using Printf
    using DataInterpolations
    using Roots 
    
    import Plots: plot
    using Plots

    export EquationOfState, BarotropicEOS, polytropic_eos
    export bblTableEOS, bblPolytropeEOS
    export addCrustByJumping, addCrustByGlueing, bbpnvCrustEOS
    export writeHDF5BarotropicEOS, readHDF5BarotropicEOS
    export getSymbol, getUnit, getDimension, getLabels
    export interpolated
    export convertUnits
    export plot

    # unit conversion utilities
    using UniCon

    export UnitSystem, Unit, Unitful, SI, CGS, AstroUnits, NuclearUnits, GeometrizedUnits
    export DimensionLessDim, LengthDim, TimeDim, MassDim, TemperatureDim, VolumeDim, DensityDim
    export EnergyDim, VelocityDim, EnergyDensityDim, PressureDim, NumberDensityDim
    export convertValue
    export docRep

    NumberDensityDim() = DimensionLessDim() / VolumeDim()

    # define abstract types
    abstract type EquationOfState end

    # define constant values
    const eos_fields = Dict(
        :energy_density => "e",
        :pressure => "p",
        :number_density => "n",
        :speed_of_sound => "c_s",
        :de_dp => "de_dp",
        :d2e_pd2 => "d2e_dp2",
        :d3e_pd3 => "d3e_dp3"
    )

    include("polytrope_eos.jl")
    include("barotropic_eos.jl")
    include("bombaci_logoteta.jl")
    include("add_crust.jl")

end
