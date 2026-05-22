# CompOSE Data Files

This directory stores manually downloaded and unpacked EoS tables from
[CompOSE](https://compose.obspm.fr). The notes below summarize the parts of the
CompOSE table format that are currently relevant for EOSeNS import code.

Content of this readme is based on the information provided in the CompOSE manual, 
which can be downloaded [here](https://compose.obspm.fr/manual/).

## Parameter Grid Files

CompOSE stores the independent variables of an EoS table in separate parameter-grid
files. For a generic EoS named `eos`, the standard independent quantities/grid files are:

- `eos.t` : temperature grid, `T`
- `eos.nb`: baryon number density grid, `n_b`
- `eos.yq`: charge-fraction grid, `Y_q`

These files are required by the CompOSE format even when the table is effectively
one-dimensional or two-dimensional. For cold beta-equilibrated neutron-star matter,
the file that matters most for the current EOSeNS importer is `eos.nb`.

The grid files use the same basic storage convention:

1. first row: minimum grid index
2. second row: maximum grid index
3. following rows: grid values for all indices from minimum to maximum

For example, the SLy4 `eos.nb` begins with:

```text
0
196
1.000000e-07
1.096478e-07
...
```

Here the density grid indices run from `0` to `196`, and the baryon density values
start on the third numeric line. The density unit is `fm^-3`.

## Thermodynamic Table: `eos.thermo`

`eos.thermo` is the required CompOSE table containing thermodynamic quantities.
This file represents a three-dimensional table i.e. it has 3 indices. 
The first row is a header with three entries:

```text
m_n  m_p  I_l
```

- `m_n`: neutron mass used by this EoS model, in `MeV`
- `m_p`: proton mass used by this EoS model, in `MeV`
- `I_l`: lepton flag

If electrons and/or muons are included, `I_l = 1`. Otherwise, CompOSE treats the
table as not containing leptons, and the effective lepton chemical potential is taken to
be zero.

Every following row gives the grid indices and thermodynamic values for one table
point:

```text
i_T  i_nb  i_Yq  Q1  Q2  Q3  Q4  Q5  Q6  Q7  N_add  q1  q2  ...  qN_add
```

The first three entries identify the grid indices:

- `i_T` : index into `eos.t`
- `i_nb`: index into `eos.nb`
- `i_Yq`: index into `eos.yq`

The seven required thermodynamic quantities are:

| quantity | meaning | unit |
| --- | --- | --- |
| `Q1` | `p / n_b` | `MeV` |
| `Q2` | `s / n_b` | dimensionless |
| `Q3` | `mu_b / m_n - 1` | dimensionless |
| `Q4` | `mu_q / m_n` | dimensionless |
| `Q5` | `mu_l / m_n` | dimensionless |
| `Q6` | `f / (n_b m_n) - 1` | dimensionless |
| `Q7` | `e / (n_b m_n) - 1` | dimensionless |

where:

- `p` is pressure
- `s` is entropy density
- `mu_b` is baryon chemical potential
- `mu_q` is charge chemical potential
- `mu_l` is effective lepton chemical potential
- `f` is free energy density
- `e` is internal energy density
- `n_b` is baryon number density
- `m_n` is the neutron mass from the first row of `eos.thermo`

`N_add` is the number of optional extra thermodynamic quantities stored after `Q7`.
Their meanings are model-specific and should be described in the corresponding
CompOSE data sheet. If duplicate rows with the same `(i_T, i_nb, i_Yq)` occur, the
manual specifies that the values from the row read last should be used.

## Current EOSeNS Reconstruction

Currently, `EOSeNS` treats only 1d tables of cold neutron star matter.
This means we assume $T=0$ and $Y_q=\mathrm{const.}$ througout the table. 

For cold barotropic imports, EOSeNS currently reads `eos.nb` and `eos.thermo` and
reconstructs:

```text
n_b = eos.nb[i_nb]
p   = Q1 * n_b
e   = (Q7 + 1) * n_b * m_n
```

The resulting `BarotropicEOS` uses `NuclearUnits()`:

- `number_density`: `fm^-3`
- `pressure`: `MeV/fm^3`
- `energy_density`: `MeV/fm^3`

The current importer is intended for one-dimensional cold/barotropic slices, such as
cold beta-equilibrated neutron-star matter. Full finite-temperature tables with
multiple `T` or `Y_q` values need a richer tabular EoS representation.
