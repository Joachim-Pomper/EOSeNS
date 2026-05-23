# Derivation of the Tolman VII metric

Starting from the definition of the line element

$$
ds^2 = -f(r)c^2dt^2 + \frac{dr^2}{h(r)} + r^2d\Omega^2.
$$

and the energy density, the functions $m(r)$ and $h(r)$ are easily derived. Here we provide a detailed discussion of the derivation of $f(r)$ and $p(r)$. We use the foloowing parameter

$$
\xi = \frac{r}{R},
\qquad
C = \frac{G M}{R c^2},
\qquad
\kappa = \frac{8\pi G}{c^4},
$$

as introduced in [`README.md`](README.md).


##  Derivation of pressue and metric coefficient $f(r)$ 

The remaining metric function $f(r)$ is fixed by the pressure equation amont the TOV equations. 
For the metric convention stated above, the radial Einstein equation gives

$$
\kappa p(r)
=
\frac{h(r)-1}{r^2}
+
\frac{h(r)}{r}\frac{d}{dr}\log f(r).
$$

In terms of the dimensionless radius $\xi = r/R$,

$$
\kappa R^2 p(\xi)
=
\frac{h(\xi)-1}{\xi^2}
+
\frac{h(\xi)}{\xi}\frac{d}{d\xi}\log f(\xi).
$$

Using the Tolman VII form of $h$,

$$
\frac{h(\xi)-1}{\xi^2}
=
-2C\frac{5 - 3\mu\xi^2}{5 - 3\mu}.
$$

For $\mu > 0$, define

$$
q =
\sqrt{\frac{6\mu C}{5 - 3\mu}}.
$$

The useful observation is that the combination

$$
Y(\xi)
=
\sqrt{h(\xi)}
+
q\left(\xi^2 - \frac{5}{6\mu}\right)
$$

has the logarithmic derivative

$$
\frac{1}{2}\frac{d}{d\xi}\log Y(\xi)
=
\frac{q\xi}{\sqrt{h(\xi)}}.
$$

This motivates defining a phase variable $z(\xi)$ by normalizing this logarithm at
the surface:

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

Then

$$
\frac{dz}{d\xi}
=
\frac{q\xi}{\sqrt{h(\xi)}}.
$$

Writing $f$ as a squared sine of this phase,

$$
f(\xi) = A\sin^2 z(\xi),
$$

gives

$$
\frac{d}{d\xi}\log f(\xi)
=
2\frac{dz}{d\xi}\cot z(\xi)
=
\frac{2q\xi}{\sqrt{h(\xi)}}\cot z(\xi).
$$

The pressure is therefore

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

The surface of the star is defined by $p(1) = 0$. Since $h(1) = 1 - 2C$, this
fixes the surface phase:

$$
\cot z(1)
=
\frac{C}{q\sqrt{1 - 2C}}.
$$

Equivalently, on the principal branch,

$$
z(1)
=
\arctan\left(
\frac{q\sqrt{1 - 2C}}{C}
\right).
$$

The remaining normalization constant $A$ is fixed by matching to the exterior
Schwarzschild metric:

$$
f(1) = 1 - 2C.
$$

Thus

$$
A =
\frac{1 - 2C}{\sin^2 z(1)},
$$

and the final expression for the time-time metric coefficient is

$$
f(\xi)
=
(1 - 2C)
\frac{\sin^2 z(\xi)}{\sin^2 z(1)}.
$$

This construction makes the two physical boundary conditions manifest:

$$
p(1) = 0,
\qquad
f(1) = h(1) = 1 - 2C.
$$

## Constant density start limit

The constant-density case $\mu = 0$ should be taken as a separate limit. In that
case the interior Schwarzschild expressions are

$$
f(\xi)
=
\frac{1}{4}
\left[
3\sqrt{1 - 2C}
-
\sqrt{1 - 2C\xi^2}
\right]^2,
$$

and

$$
p(\xi)
=
\epsilon_c
\frac{
\sqrt{1 - 2C\xi^2} - \sqrt{1 - 2C}
}{
3\sqrt{1 - 2C} - \sqrt{1 - 2C\xi^2}
}.
$$
