function setupTOVeq(eos, p_central, r_max = 30e5)
    
    # convert units, if necessary
    if eos.unit_system != CGS()
        eos = convertUnits(eos, CGS())
    end
    
    # The equation of state is given in a nuclear units, the tov equation should be formulated in the SI system
    e_from_p = interpolated(eos, "energy_density", "pressure")
    
    # Check if central pressure is compatible with the equation of state.
    function inPressureDomain(p)
        indomain = (p >= minimum(eos["pressure"]))
        indomain = indomain & (p <= maximum(eos["pressure"]))
        return indomain
    end
    
    if !inPressureDomain(p_central)
        throw(DomainError(p_central, "Central pressure value not supported by provided equation of state"))
    end
    
    # u[1] = m(r)       mass in g
    # u[2] = p(r)/c^2   rescaled pressure in g/cm^3 
    # u[3] = α(r)       dimensionless gravitational potential
    function tovEquation!(du, u, intern_param, r)

        try
            p = u[2]*c_light^2 # pressure in g/(s^2 cm)
            rho = e_from_p(p)/c_light^2 # mass_density in g/cm^3
            
            #calculate dα/dr
            if r>0. # regular case
                r_minus_m = (r/g_newton*c_light^2) - 2.0*u[1]
                dalpha_dr = (u[1]+4.0*pi*r^3*u[2])/(r_minus_m*r) 
            elseif r==0. # limiting case of zero radius
                dalpha_dr = 0
            else # not valid, r>0 must hold
                throw(DomainError(r, "r>=0 must hold!"))
            end
    
            #set values
            du[1] = 4.0*pi*r^2*rho
            du[2] = -(rho+u[2])*dalpha_dr
            du[3] = dalpha_dr
    
        catch error
            # We are now outside the description the EOS allows, abort the integration
            # Must assume that EOS is suitable such that this means we are also outside the Neutron star       
            if !inPressureDomain(u[2]*c_light^2) 
                du[1] = 0 # Mass will no longer increase
                du[2] = 0
                du[3] = 0
            else
                print(u[2]*c_light^2)
                throw(error)
            end
        end
    end
    
    # Set up initial value problem
    u0 = [0., p_central/c_light^2, 0.]
    
    # Interval of radius values
    r_span = (0.0, r_max)

    return (ODEProblem(tovEquation!, u0, r_span), eos)
end


function solveTOVeq(tov_problem::TovProblem, solver; kwargs...)

    # convert units if necessary   
    eos = tov_problem[2]
    if eos.unit_system != CGS()
        eos = convertUnits(eos, CGS())
    end
    
    # use DifferentialEquations.jl to solve the ODE
    sol = solve(tov_problem[1], solver; kwargs...)
    
    # Mass, radius and pressure
    r = sol.t[:] 
    m = sol[1,:]
    p = sol[2,:]*c_light^2
    
    idx_border = findfirst(diff(sol[1,:]) .== 0)
    radius = r[idx_border]
    mass = m[idx_border]
    p[idx_border:end] .= NaN
    
    # get energy_density 
    e_from_p = interpolated(eos, "energy_density", "pressure")
    e = e_from_p(p)

    #alpha and beta must be calculated 
    beta = -log.(1.0 .- 2.0 .* g_newton .* m ./ (r .* c_light^2)) ./ 2.0
    alpha = sol[3,:] .- beta[idx_border] .- sol[3,idx_border]

    return TovSolution(mass, radius, r ,m, p, e, alpha, beta)
end

function getMassRadiusRelation(eos::EquationOfState; n_points::Int64 = 256, pressure_min::Float64 = 1e34)
    
    mass = Vector{Float64}()
    radius = Vector{Float64}()

    # get
    if eos.unit_system != CGS()
        eos = convertUnits(eos, CGS())
    end
    central_pressure = eos["pressure"]
    
    # truncate pressure values such that they are above the minimum pressure
    min_pressure = convertValue(pressure_min, PressureDim(), CGS(), eos.unit_system)
    central_pressure = central_pressure[central_pressure .>= min_pressure]

    #refine pressure grid in case there are not enough points
    while length(central_pressure) < n_points
        central_pressure = refineGrid1d(central_pressure)
    end

    #solve the differential equation and extract mass and radius
    for p_c in central_pressure
        
        solution = solve(setupTOVeq(eos, p_c)[1], Tsit5(), reltol = 1e-8, dt = .01)

        idx_border = findfirst(diff(solution[1,:]) .== 0)
        m = solution[1,idx_border] 
        r = solution.t[idx_border] / 1e5
        append!(mass, [m / 1.988416e33 ])
        append!(radius, [r])
    end

    return mass, radius, central_pressure
end
