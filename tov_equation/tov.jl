module TOV

    using HDF5
    using Printf
    using DataInterpolations
    using DifferentialEquations
    using Plots
    
    include("../utils/utils.jl")
    using .Utils
    
    include("../equation_of_state/neos.jl")
    using .NEOS

    export TovProblem, TovSolution
    export setupTOVeq, solveTOVeq, getMassRadiusRelation
    
    # working in cgs units
    c_light = 2.99792458e10  # cm/s
    g_newton = 6.6730831e-8  # cm^3 x s^-2 x g^-1
    m_sol = 1.988416e33 # solar mass in g
    
    TovProblem = Tuple{
        SciMLBase.AbstractODEProblem{Vector{Float64}, Tuple{Float64, Float64}, true},  EquationOfState
    }
    
    struct TovSolution
        mass::Float64            # mass of the star
        radius::Float64          # radius of the star
        r::Vector{Float64}       # radial coordinate in cm
        m::Vector{Float64}       # mass within radius in g
        e::Vector{Float64}       # energy_density in g / cm s^2
        p::Vector{Float64}       # pressure in g / cm s^2
        alpha::Vector{Float64}   # dimensionless
        beta::Vector{Float64}    # dimesnionless
    end
    
    include("tov_eq.jl")


end