# Tollmann-Openheimer-Volkoff (TOV) equation 

The Tollmann-Openheimer-Volkoff (TOV) equation is given as follows 

$$
\begin{aligned}
    \frac{\partial m(r)}{\partial r} &= 4\pi r^2\frac{\rho(r)}{c^2} \\
    \frac{\partial P(r)}{\partial r} &= -(P(r)+\rho(r))\frac{G_N(4 \pi r^3 P(r) + c^2 m(r))}{(rc^2 - 2G_N m(r)) c^2r}\\
\end{aligned}
$$

To be fully solved, they need to be supplemented by the Equation of state $\rho(r) = \rho(P(r))$. Initial conditions can then be set by putting

$$
\begin{aligned}
    m(0) &= 0  \\
    P(0) &= P_c
\end{aligned}
$$

One free parameter remains, that must be choosen in accordance with the equation of state. This is the most efficient manuveur to calculate the mass-radius relations of neutron stars, by integrating the equations until $P(r=R)=0$ is reached, and mass and radius can be determined. 

If we are also interested into the corresponding metric 

$$
g = -e^{2\beta(r)}\mathrm{d}t^2 + e^{2\alpha(r)}\mathrm{d}r^2 + \mathrm{d}\Omega(\theta, \phi) 
$$

we can use the following linear differential equations for $\alpha$ and $\beta$.

$$
\begin{aligned}
\frac{\partial \alpha(r)}{\partial r} &= \frac{G_N(4 \pi r^3 P(r) + c^2 m(r))}{(rc^2 - 2G_N m(r)) c^2r} \\
    \frac{\partial \beta(r)}{\partial r} &=  \frac{G_N(4 \pi r^3 \rho(r) - c^2 m(r))}{(rc^2 - 2G_N m(r)) c^2r}
\end{aligned}
$$

The equation for $\alpha$ can be solved along the others with initial conditions $\alpha(0) = 0$. The physical $\alpha(r)$ must satisfy that $\alpha(R) = -\beta(R)$  (Since we must match with Scharzschild spacetime outside the neutron star). Due to the linearity of the equation, we can enforce this by simply adding the integration constant. 
For $\beta(r)$, it is more efficient to use the relation $e^{-2\beta(r)} = 1-\frac{2 G_N m(r)}{c^2 r}$ with the mass function. This gives

$$
\beta(r) = -\frac{1}{2}\ln\left(1-\frac{2 G_N m(r)}{c^2 r}\right) \,.
$$

A numerical problem can occur at the beginning. Here, the equations have the problem that they approach an undefined limit, since $\lim_{r\to 0} m(r) = 0$. However, we know the behaviour of the derivarives and we know that $\rho(r =0) = \rho(P_0) < \infty$ and $P(r=0) = P_0 < \infty$. Hence, we obtain that

$$
\lim_{r\to 0}  \frac{\partial m(r)}{\partial r} = \lim_{r\to 0} 4\pi r^2\frac{\rho(r)}{c^2}
$$

From this we can further deduce, using L'Hospital, that 

$$
\begin{aligned}
    \lim_{r\to 0}  \frac{m(r)}{r} &= \lim_{r\to 0} 4\pi r^2\frac{\rho(r)}{c^2} = 0\\
    \lim_{r\to 0}  \frac{m(r)}{r^2} &= \lim_{r\to 0} 2\pi r\frac{\rho(r)}{c^2} = 0
\end{aligned}
$$

Then, for the derivative of $\alpha$, and thus also the pressure, we get 

$$
\begin{aligned}
    \lim_{r\to 0} \frac{\partial \alpha(r)}{\partial r}
    &= \lim_{r\to 0} \frac{G_N(4 \pi r^3 P(r) + c^2 m(r))}{(rc^2 - 2G_N m(r)) c^2r} \\
    &= \lim_{r\to 0} \frac{G_N\left(4 \pi r P(r) + c^2 \frac{m(r)}{r^2}\right)}{\left(c^2 - 2G_N \frac{m(r)}{r}\right) c^2} = 0
\end{aligned}
$$

## Code Implementation

The implementation is split across `tov.jl` and `tov_eq.jl`. 

* `tov.jl`: Main module file. Defines the physical constants, the `TovSolution` container, and exports the public functions. 
* `tov_eq.jl`: Contains the actual solver logic.

### Solver logic

#### `setupTOVeq(eos, p_central, r_max)`

Prepares the initial-value problem. 

* It first converts the EOS to CGS units if needed. * Builds `energy_density(pressure)` via interpolation.
* Checks that the central pressure is inside the EOS
domain
* Defines the ODE system solved by `DifferentialEquations.jl`.

The ode is defined in terms of a vector variable `u`, whose entries have the following meaning.

```julia
u[1] = m(r)
u[2] = P(r) / c^2
u[3] = alpha(r)
```

At `r = 0`, the code uses the regular limiting value
`dalpha_dr = 0` to avoid the singular-looking central expression, see above.

#### `solveTOVeq(tov_problem, solver; kwargs...)` 

Solve the ODE.

* runs the ODE solver
* reconstructs
`P(r)`, 
* finds the stellar surface where the mass stops changing
* masks the pressure outside the star with `NaN`. 

At the end it then computes

```math
\beta(r) = -\frac{1}{2}\log\left(1-\frac{2G_Nm(r)}{rc^2}\right),
```

from the mass profile, and shifts the integrated `alpha` by a constant so the interior metric matches the exterior Schwarzschild solution at the surface.

#### `getMassRadiusRelation(eos; ...)` 

Iterates over `solveTOVeq()` for a grid of central pressures and returns the corresponding masses, radii, and central-pressure
values. Masses are returned in solar masses and radii in kilometers.
