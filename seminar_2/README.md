# Exoplanets and Kepler's Third Law

All credits for calculations go to Dr Andrew French, Winchester College.
http://www.eclecticon.info/index_htm_files/BPhO%20CompPhys%2006%20Planets.pdf

---

## 0. Table of contents
- [Exoplanets and Kepler's Third Law](#exoplanets-and-keplers-third-law)
	- [0. Table of contents](#0-table-of-contents)
	- [1. Kepler III in the Solar System](#1-kepler-iii-in-the-solar-system)

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
	& = 1.000 \space\textrm{AU}^{-3}\space\textrm{YR}^2 \\[5pt]
\end{align*}
```

![Kepler III linear regression](./images/kepler_III_solar.png "Kepler III linear regression")

The gradient is approximately $1.000$ as expected.
