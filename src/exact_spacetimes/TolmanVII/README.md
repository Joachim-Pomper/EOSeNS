# Tolman VII Spacetime

This directory implements the Tolman VII interior solution.  

## Theory

Throughout this note we use

$$
\xi = \frac{r}{R},
\qquad
C = \frac{G M}{R c^2},
\qquad
\kappa = \frac{8\pi G}{c^4},
$$

where `R` is the stellar radius, `M` is the total gravitational mass, and `C` is
the compactness. The constant `kappa` keeps the relativistic coupling compact in
the density and pressure formulas.

The static, spherically symmetric line element is given by

$$
ds^2 = -f(r)c^2dt^2 + \frac{dr^2}{h(r)} + r^2d\Omega^2.
$$

### Parameters

The generalized Tolman VII model is controlled by

$$
(R, C, \mu).
$$

Here `mu` is the self-boundness parameter. It fixes the surface density relative
to the central density:

$$
\epsilon(R) = \epsilon_c(1 - \mu).
$$

Special cases:

- `mu = 1` gives the usual Tolman VII profile with zero surface density.
- `0 < mu < 1` gives a self-bound star with nonzero surface density.
- `mu = 0` gives the constant-density star.

### Density

The model has a parabolic energy-density profile:

$$
\epsilon(r) = \epsilon_c\left(1 - \mu \xi^2\right).
$$

The central density is fixed by requiring $m(R) = M$:

$$
\epsilon_c =
\frac{30 C}{\kappa R^2(5 - 3\mu)}.
$$

Equivalently,

$$
\epsilon(r) =
\frac{30 C}{\kappa R^2(5 - 3\mu)}
\left(1 - \mu \xi^2\right).
$$

### Mass and metric coefficient $h(r)$ 

The enclosed mass follows from integrating the density:

$$
c^2m(r) =
4\pi\int_0^r \epsilon(s)s^2\,ds
$$

With the density above,

$$
c^2m(r) =
Mc^2\,
\frac{5\xi^3 - 3\mu\xi^5}{5 - 3\mu}
=
\frac{8\pi C R}{\kappa}\,
\frac{5\xi^3 - 3\mu\xi^5}{5 - 3\mu}.
$$

This satisfies $m(0) = 0$ and $m(R) = M$.
The radial metric function can be written in terms of the enclosed mass:

$$
h(r) =
1 - \frac{\kappa c^2m(r)}{4\pi r}.
$$

Using the Tolman VII mass profile,

$$
h(r) =
1 -
2C\xi^2
\frac{5 - 3\mu\xi^2}{5 - 3\mu}.
$$

### Pressure and metric coefficient $f(r)$

The remaining metric function $f(r)$ is derived by the pressure TOV equation. A derivation of the rather involved formulas can be found [here](tolman_vii_derivation.md).  

The pressure is given by 

$$
p(\xi)
=
\frac{2}{\kappa R^2}
\left[
q\sqrt{h(\xi)}\cot z(\xi)
-
C\frac{5 - 3\mu\xi^2}{5 - 3\mu}
\right].
$$

and the remaining metric function as 

$$
f(\xi)
=
(1 - 2C)
\frac{\sin^2 z(\xi)}{\sin^2 z(1)}.
$$

We introduced the auxilary variable 

$$
z(\xi)
=
z(1)
+
\frac{1}{2}
\log\left[
\frac{
\sqrt{h(\xi)} + q\left(\xi^2 - \frac{5}{6\mu}\right)
}{
\sqrt{1 - 2C} + q\left(1 - \frac{5}{6\mu}\right)
}
\right].
$$

with

$$
\cot z(1)
=
\frac{C}{q\sqrt{1 - 2C}}.
$$

### Equation of State

Following the Tolman VII equation-of-state parametrization discussed by
Raghoonundun and Hobill, the density profile can be inverted for $\mu > 0$:

$$
\xi^2(\epsilon) =
\frac{1 - \epsilon / \epsilon_c}{\mu}.
$$

This is valid for

$$
\epsilon_c(1 - \mu) \leq \epsilon \leq \epsilon_c.
$$

For compact notation in the EOS, define

$$
h_\epsilon =
h\left(R\xi(\epsilon)\right),
\qquad
z_\epsilon =
z\left(\xi(\epsilon)\right).
$$

The Tolman VII barotropic equation of state is then

$$
p(\epsilon) =
\frac{2}{\kappa R^2}
\left[
q\sqrt{h_\epsilon}\cot z_\epsilon
-
C
\frac{5 - 3\mu\xi^2(\epsilon)}{5 - 3\mu}
\right].
$$

Equivalently, the EOS can be treated parametrically as

$$
\epsilon = \epsilon(\xi),
\qquad
p = p(\xi),
\qquad
0 \leq \xi \leq 1.
$$

This form is often the most convenient implementation path, because it reuses
the same numerically stable radial expressions for $\epsilon(r)$ and $p(r)$.
For $\mu = 0$, the density is constant, so the relation between pressure and
density is degenerate rather than a single-valued barotropic function.

### TOV Potentials

The exact-spacetime code provides the same dimensionless potentials used by
the numerical TOV solver:

$$
ds^2 =
-\exp\left(2\alpha(r)\right)c^2dt^2
+ \exp\left(2\beta(r)\right)dr^2
+ r^2d\Omega^2.
$$

They are defined from $f$ and $h$ by

$$
\alpha(r) = \frac{1}{2}\log f(r),
\qquad
\beta(r) = -\frac{1}{2}\log h(r).
$$

These $\alpha$ and $\beta$ are not auxiliary Tolman VII parameters. They are only
the metric potentials used by the TOV equation.

## Implementation Details