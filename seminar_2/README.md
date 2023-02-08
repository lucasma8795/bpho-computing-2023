# Exoplanets and Kepler's Third Law

All credits for calculations go to Dr Andrew French, Winchester College.
http://www.eclecticon.info/index_htm_files/BPhO%20CompPhys%2006%20Planets.pdf

---

## 0. Table of contents
- [Exoplanets and Kepler's Third Law](#exoplanets-and-keplers-third-law)
	- [0. Table of contents](#0-table-of-contents)
	- [1. Kepler III in the Solar System](#1-kepler-iii-in-the-solar-system)
	- [2. Exoplanets](#2-exoplanets)

---
## 1. Kepler III in the Solar System
> Kepler's Third Law: the squares of the orbital periods of the planets are directly proportional to the cubes of the semi-major axes of their orbits. ([NASA](https://solarsystem.nasa.gov/resources/310/orbits-and-keplers-laws/))

We define the following symbols:
| Symbol | Name | Unit | Description |
| --- | --- | --- | --- |
| $P$ | Orbital period | $\textrm{days}$ | Time taken for astronomical object to complete one orbit around another object |
| $M$, $m$ | Mass | $\textrm{kg}$ | Mass of a star; mass of a planet |
| $a$ | [Semi-major axis](https://en.wikipedia.org/wiki/Semi-major_and_semi-minor_axes) | $\textrm{m}$ | Length of the longest radius in an ellipse |
| $b$ | Semi-minor axis | $\textrm{m}$ | Length of the shortest radius in an ellipse |
| $\epsilon$ | [Eccentricity](https://en.wikipedia.org/wiki/Eccentricity_(mathematics)) | -- | Uniquely characterizes the shape of a conic section |

...and the following constants:
| Constant | Name | Value |
| --- | --- | --- |
| $M_\odot$ | Mass of Sun | $1.989 \times 10^{30}\space\textrm{kg}$ |
| $m_\oplus$ | Mass of Earth | $5.972 \times 10^{24}\space\textrm{kg}$ |
| $m_J$ | Mass of Jupiter | $1.898 \times 10^{27}\space\textrm{kg}$ |
| $G$ | [Gravitational constant](https://en.wikipedia.org/wiki/Gravitational_constant) | $6.674 \times 10^{-11}\space\textrm{m}^3\space\textrm{kg}^{-1}\space\textrm{s}^{-2}$ |

...and the following units:
| Unit | Name | Value |
| --- | --- | --- |
| $1\space\textrm{YR}$ | Days in a year | $365.242\space\textrm{days}$ |
| $1\space\textrm{AU}$ | [Earth-sun distance](https://en.wikipedia.org/wiki/Astronomical_unit) | $1.496 \times 10^{11}\space\textrm{m}$ |

Kepler's Third Law can be now expressed as an equation:
```math
\begin{align*}
	P^2 & = \frac{4\pi^2}{G(M+m)}a^3 \\[5pt]
	P^2 & \propto a^3
\end{align*}
```

In the context of Earth in the Solar System:
```math
\begin{align*}
	\textrm{YR}^2 & = \frac{4\pi^2}{G(M_\odot+m_\oplus)}\textrm{AU}^3 \\[5pt]
	& \approx \frac{4\pi^2}{GM_\odot}\textrm{AU}^3 \\[5pt]
	\frac{4\pi^2}{GM_\odot} & \approx 2.974\times10^{-19}\space\textrm{m}^{-3}\space\textrm{s}^2 \\[5pt]
	& = 133385 \space\textrm{AU}^{-3}\space\textrm{days}^2 \\[5pt]
	& = 1.000 \space\textrm{AU}^{-3}\space\textrm{YR}^2
\end{align*}
```

![Kepler III linear regression](./images/kepler_III_solar.png "Kepler III linear regression")

The gradient is approximately $1.000$ as expected.

---
## 2. Exoplanets
We obtain exoplanet data from [The Extrasolar Planets Encyclopedia](http://www.exoplanet.eu/).
After pre-processing, 1012 exoplanets remain.

```
>> exoplanet_data
     exoplanet_name     star_name       M  M_errmin  M_errmax         m  m_errmin  m_errmax              P       P_errmin       P_errmax          a  a_errmin  a_errmax
0          11 Oph b        11 Oph  0.0162     0.005     0.005  21.00000    3.0000    3.0000  730000.000000  365000.000000  365000.000000  243.00000   55.0000   55.0000
1      2M 2140+16 b    2M 2140+16  0.0800     0.060     0.060  20.00000   20.0000   80.0000    7340.000000     584.000000     584.000000    3.53000    0.1500    0.1500
2      2M 2206-20 b    2M 2206-20  0.1300     0.050     0.050  30.00000   20.0000   70.0000    8686.000000      69.400000      69.400000    4.48000    0.4000    0.4000
3       2M1059-21 b     2M1059-21  0.8100     0.008     0.002  66.95000    4.8400    4.4100     690.000000       3.000000       3.400000    0.80000    0.0200    0.0100
4          51 Eri b        51 Eri  1.7500     0.050     0.050   2.60000    0.3000    0.3000   10260.000000    1800.000000    6300.000000   11.10000    1.3000    4.2000
..              ...           ...     ...       ...       ...       ...       ...       ...            ...            ...            ...        ...       ...       ...
983   gamma 1 Leo b   gamma 1 Leo  1.2300     0.210     0.210  63.88000   56.1000   72.9100     428.500000       1.250000       1.250000    1.19000    0.0200    0.0200
984  gamma Cephei b  gamma Cephei  1.4000     0.120     0.120  17.53900    0.7450    0.7450     905.640000       2.830000       2.830000    2.14590    0.0048    0.0048
985     kappa And b     kappa And  2.8000     0.200     0.100  13.00000    2.0000   12.0000  215000.000000  100000.000000  100000.000000  100.00000   46.0000   27.0000
986        pi Men b        pi Men  1.0940     0.039     0.039  12.60000    2.0000    2.0000    2088.330000       0.340000       0.340000    3.30800    0.0390    0.0390
987        pi Men c        pi Men  1.0940     0.039     0.039   0.01142    0.0012    0.0012       6.267852       0.000016       0.000016    0.06702    0.0005    0.0005
```
Note: 
Unit of `M` is $M_\odot$, and unit of `m` is $m_J$.
Unit of `a` is $\textrm{AU}$, and unit of `P` is $\textrm{YR}$.

We rearrange Kepler III to get the following:
```math
\begin{align*}
	P^2 & = \frac{4\pi^2}{G(M+m)} \\[5pt]
	\textrm{YR}^2 & = \frac{4\pi^2}{G(M_\odot+m_\oplus)} \approx \frac{4\pi^2}{GM_\odot} (M_\odot \gg m_\oplus) \\[5pt]
	\therefore \left(\frac{P}{\textrm{YR}}\right)^2 & = \left(\frac{a}{\textrm{AU}}\right)^3 \left(\frac{M_\odot}{M+m}\right) \\[5pt]
	3\log\left(\frac{a}{\textrm{AU}}\right) & = 2 \log\left(\frac{P}{\textrm{YR}}\right) + \log\left(\frac{M}{M_\odot} + \frac{m}{m_J}\times\frac{m_J}{M_\odot}\right)
\end{align*}
```

There exists a relationship containing $a$, $P$, $M$ and $m$.

![Kepler III on exoplanets](./images/kepler_III_exoplanets.png "Kepler III on exoplanets")

The gradient is, again, approximately $1.000$ as expected.