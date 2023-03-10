import signal
from math import floor, sqrt
from sys import maxsize as INF
from typing import List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats
import toml
from sklearn.linear_model import LinearRegression


# physics constants
G 			= 6.674e-11 # m3 kg-1 s-2
AU 			= 1.496e11 # m
YR 			= 365.2422 # days
M_sun 		= 1.989e30 # kg
M_jupiter 	= 1.898e27 # kg
M_earth 	= 5.972e24 # kg


def plot_kepler_III(planet_data:pd.DataFrame, planet_names:List[str]):
	# The square of the orbital period of a planet is directly proportional to
	# 	the cube of the semi-major axis of its orbit
	X = (planet_data.a ** 3).values.reshape((-1, 1))
	y = planet_data.P ** 2
	reg = LinearRegression(fit_intercept=False).fit(X, y)
	m = reg.coef_[0]

	plt.figure("Kepler's Third Law")
	plt.title("Kepler's Third Law")
	plt.xlabel(r"$(a / AU)^3$")
	plt.ylabel(r"$(P / YR)^2$")
	plt.scatter(X, y, s=50, marker="+", color="k", linewidths=0.8)
	# planet name annotations
	for idx, row in planet_data.iterrows():
		plt.annotate(planet_names[idx], (row.a**3, row.P**2))
	# regression line
	plt.plot(np.linspace(0, X.max(), 2), np.linspace(0, X.max(), 2)*m, linewidth=1, color="k", linestyle="--")


def plot_kepler_III_modified(df:pd.DataFrame):
	X = 3 * np.log(df.a)
	X_min = 3 * np.log(df.a - df.a_errmax)
	X_max = 3 * np.log(df.a + df.a_errmax)

	tmp = np.log(df.M + df.m * M_jupiter / M_sun)
	y = 2 * np.log(df.P / YR) + tmp
	y_min = 2 * np.log((df.P - df.P_errmax) / YR) + tmp
	y_max = 2 * np.log((df.P + df.P_errmax) / YR) + tmp

	plt.figure("Kepler's Third Law (modified)")
	plt.title("Kepler's Third Law on exoplanets")
	plt.xlabel(r"$3 \log\left(\frac{a}{AU}\right)$")
	plt.ylabel(r"$2 \log\left(\frac{P}{YR}\right) + \log\left(\frac{M}{M_\odot} + \frac{m}{m_J}\times\frac{m_J}{M_\odot}\right)$")
	
	plt.errorbar(X, y, xerr=[X-X_min, X_max-X], yerr=[y-y_min, y_max-y], marker="+", c="k", elinewidth=1, lw=0, zorder=1, capsize=2)
	
	plt.plot(np.linspace(X.min(), X.max(), 2), np.linspace(X.min(), X.max(), 2), color="r", linestyle="--", lw=1)
	plt.legend(["Kepler III", "Data"])


def main():
	# enable ctrl-c program termination
	signal.signal(signal.SIGINT, signal.SIG_DFL)

	# read in planet data
	with open("./data/planets.csv", "r") as fo:
		planet_data = pd.read_csv(fo)

	# read config file
	# CONFIG = toml.load("./config.toml")

	planet_names = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]

	# relationship between semi-major axis and orbital period
	plot_kepler_III(planet_data, planet_names)

	# read in exoplanet data
	with open("./data/exoplanets.csv", "r") as fo:
		exoplanet_data = pd.read_csv(fo)

	# draw relationship between a, P, M and m with modified Kepler III
	plot_kepler_III_modified(exoplanet_data)

	plt.show()


if __name__ == "__main__":
	main()