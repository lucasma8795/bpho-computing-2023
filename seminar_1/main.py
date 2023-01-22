import signal
from math import floor, sqrt
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
import scipy.stats


def find_rho(df:pd.DataFrame) -> Tuple[float, float]:
	# rho ln(S_0/S) = I_0 + S_0 - I - S
	# given a df of collected data, perform linear regression to find slope rho
	# returns (m, c) where the regression line is mx + c

	I_0, S_0 = df.iloc[0].I, df.iloc[0].S

	X = np.log(S_0 / df.S).values.reshape((-1, 1))
	y = I_0 + S_0 - df.I - df.S

	reg = LinearRegression(fit_intercept=False).fit(X, y)

	return reg.coef_[0], reg.intercept_, X, y


def find_RMSE(model:pd.DataFrame, historical:pd.DataFrame, dt:float) -> float:
	total_sq_err = 0

	for _, row in historical.iterrows():
		# find corresponding entry in model
		# TODO: interpolate between two intervals
		model_ix = floor(row.t // dt)
		model_ix = max(0, model_ix)
		model_ix = min(len(model)-1, model_ix)
		
		total_sq_err += (model.iloc[model_ix].I - row.I) ** 2
		total_sq_err += (model.iloc[model_ix].S - row.S) ** 2
		total_sq_err += (model.iloc[model_ix].D - row.D) ** 2
	
	return sqrt(total_sq_err / len(historical))


def search_alpha(alpha_min:float, alpha_max:float, rho:float, historical:pd.DataFrame, dt:float, t_max:float, EPS:float=1e-7) -> float:
	# performs a ternary search on RMSE(alpha) aiming to minimize; returns an optimal alpha
	# https://en.wikipedia.org/wiki/Ternary_search
	I_0, S_0 = historical.iloc[0].I, historical.iloc[0].S

	while alpha_max - alpha_min > EPS:
		left_third = alpha_min + (alpha_max - alpha_min) / 3
		right_third = alpha_max - (alpha_max - alpha_min) / 3

		left_model = simple_eyam_model(left_third, left_third/rho, dt, t_max, I_0, S_0)
		right_model = simple_eyam_model(right_third, right_third/rho, dt, t_max, I_0, S_0)

		if find_RMSE(left_model, historical, dt) > find_RMSE(right_model, historical, dt):
			alpha_min = left_third
		else:
			alpha_max = right_third

	return (alpha_min + alpha_max) / 2
	

def simple_eyam_model(alpha:float, beta:float, dt:float, t_max:float, I_0:float, S_0:float) -> pd.DataFrame:
	SZ = floor(t_max / dt) + 1 # quantize time of length t_max into chunks of length dt

	I = np.zeros(shape=(SZ), dtype=np.float64)
	S = np.zeros(shape=(SZ), dtype=np.float64)
	D = np.zeros(shape=(SZ), dtype=np.float64)
	t = np.zeros(shape=(SZ), dtype=np.float64)
	I[0], S[0] = I_0, S_0

	for n in range(SZ-1):
		S[n+1] = S[n] - beta * S[n] * I[n] * dt						# dS/dt = -beta * S * I
		I[n+1] = I[n] + (beta * S[n] * I[n] - alpha * I[n]) * dt	# dI/dt = beta * S * I - alpha * I
		D[n+1] = D[n] + alpha * I[n] * dt							# dD/dt = alpha * I
		t[n+1] = dt * (n+1)
	
	df = pd.DataFrame(data={"t": t, "I": I, "S": S, "D": D})
	return df

rng = np.random.default_rng()
def poisson(lam):
	if lam >= 0: return rng.poisson(lam)
	else: return -rng.poisson(-lam)


def stochastic_eyam_model(alpha:float, beta:float, dt:float, t_max:float, I_0:float, S_0:float, runs:int) -> pd.DataFrame:
	# The spread of infection is a random process.
	# We can use the expected values of S, I & D changes within time interval Δt
	#	to be the mean (& variance) of a Poisson distribution.
	SZ = floor(t_max / dt) + 1 # quantize time of length t_max into chunks of length dt
	SZ_TOTAL = SZ * runs
	I_x, I_y = np.zeros(shape=(SZ_TOTAL)), np.zeros(shape=(SZ_TOTAL))
	S_x, S_y = np.zeros(shape=(SZ_TOTAL)), np.zeros(shape=(SZ_TOTAL))
	D_x, D_y = np.zeros(shape=(SZ_TOTAL)), np.zeros(shape=(SZ_TOTAL))

	for run in range(runs):
		I_y[run*SZ], S_y[run*SZ] = I_0, S_0
	
		for n in range(SZ-1):
			idx = run * SZ + n
			I_x[idx+1], S_x[idx+1], D_x[idx+1] = dt*(n+1), dt*(n+1), dt*(n+1)

			x, y = poisson(alpha * I_y[idx] * dt), poisson(beta * S_y[idx] * I_y[idx] * dt)
			S_y[idx+1] = max(0, S_y[idx] - y)			# dS/dt = -beta * S * I
			I_y[idx+1] = max(0, I_y[idx] + (y - x))		# dI/dt = beta * S * I - alpha * I
			D_y[idx+1] = max(0, D_y[idx] + x)			# dD/dt = alpha * I
	
	return I_x, I_y, S_x, S_y, D_x, D_y


def draw_eyam_model(model:pd.DataFrame, historical:pd.DataFrame, alpha:float, beta:float, dt:float) -> None:
	plt.figure("Simple Eyam Model")
	plt.title(rf"Eyam Model ($\alpha={alpha:.4f}, \beta={beta:.4f}, \Delta t={dt}$)")
	plt.xlabel("t / months")
	plt.ylabel("Eyam population")

	plt.plot(model.t, model.I, color="r")
	plt.plot(model.t, model.S, color="g")
	plt.plot(model.t, model.D, color="b")

	plt.scatter(historical.t, historical.I, s=50, marker="+", color="r", linewidths=0.8)
	plt.scatter(historical.t, historical.S, s=50, marker="+", color="g", linewidths=0.8)
	plt.scatter(historical.t, historical.D, s=50, marker="+", color="b", linewidths=0.8)

	plt.legend(["Infected (I)", "Surviving (S)", "Dead (D)"])


def draw_rho_regression(X:np.ndarray, y:np.ndarray, m:float):
	plt.figure("Alpha/beta regression line")
	plt.title("Alpha/beta regression line")
	plt.xlabel(r"X = $\ln(\frac{S_0}{S})$")
	plt.ylabel(r"y = $I_0 + S_0 - I - S$")

	max_x = X[-1][0]
	plt.scatter(X, y, s=50, marker="+", color="k", linewidths=0.8)
	plt.plot(np.linspace(0, max_x, 2), np.linspace(0, max_x, 2)*m, linewidth=1, color="k", linestyle="--")


def draw_stochastic_kde(ax, model_X:np.ndarray, model_y:np.ndarray, line_X:pd.Series, line_y:pd.Series, title:str, y_label:str):
	ax.set_title(title)
	ax.set_xlabel("t / months")
	ax.set_ylabel(y_label)

	xmin, xmax = 0, model_X.max()
	ymin, ymax = 0, model_y.max()

	xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
	positions = np.vstack([xx.ravel(), yy.ravel()])
	values = np.vstack([model_X, model_y])
	kernel = scipy.stats.gaussian_kde(values)
	f = np.reshape(kernel(positions).T, xx.shape)

	ax.imshow(np.rot90(f), cmap=plt.cm.viridis, extent=[xmin, xmax, ymin, ymax], aspect="auto")
	ax.plot(line_X, line_y, color="k", linestyle="--", linewidth=1)

	
def main():
	signal.signal(signal.SIGINT, signal.SIG_DFL)

	# symbols:
	# 	S -> susceptibles, I -> infectives, D -> dead, R -> recovered

	# read in historical data points of S, I & D
	with open("./data.csv") as fo:
		historical = pd.read_csv(fo)

	# define eyam model constants
	I_0, S_0 = historical.iloc[0].I, historical.iloc[0].S
	dt, t_max = 0.1, 5 # in months

	# define rho := alpha/beta
	rho, _, reg_X, reg_y = find_rho(historical)
	print(f"{rho=:.5f}")

	# find alpha by minimizing RMSE with historical data
	alpha = search_alpha(alpha_min=1, alpha_max=5, rho=rho, historical=historical, dt=dt, t_max=t_max)
	beta = alpha / rho
	print(f"{alpha=:.5f}, {beta=:.5f}")

	params = {
		"alpha": alpha, "beta": beta,
		"dt": dt, "t_max": t_max,
		"I_0": I_0, "S_0": S_0
	}

	# calculate & display simple eyam model
	model = simple_eyam_model(**params)
	draw_eyam_model(model=model, historical=historical, alpha=alpha, beta=beta, dt=dt)

	# calculate stochastic model
	runs = 10000
	I_x, I_y, S_x, S_y, D_x, D_y = stochastic_eyam_model(runs=runs, **params)

	# display stochastic model
	fig, axs = plt.subplots(2, 2, num="Stochastic Eyam Model")
	fig.set_figheight(7.5)
	fig.set_figwidth(10)
	fig.tight_layout(pad=4)
	fig.suptitle(rf"Stochastic Eyam Model ($\alpha={alpha:.4f}, \beta={beta:.4f}, \Delta t={dt}$, runs={runs})")

	if runs >= 1000: print(f"Drawing stochastic model with {runs} runs. This will take a while...")
	draw_stochastic_kde(
		ax=axs[0][0],
		model_X=I_x, model_y=I_y,
		line_X=model.t, line_y=model.I,
		title="Infectives",
		y_label="Infectives (I)"
	)
	draw_stochastic_kde(
		ax=axs[0][1],
		model_X=S_x, model_y=S_y,
		line_X=model.t, line_y=model.S,
		title="Susceptibles",
		y_label="Susceptibles (S)"
	)
	draw_stochastic_kde(
		ax=axs[1][0],
		model_X=D_x, model_y=D_y,
		line_X=model.t, line_y=model.D,
		title="Dead",
		y_label="Dead (D)"
	)
	axs[1][1].axis("off")

	# display alpha/beta regression line
	draw_rho_regression(X=reg_X, y=reg_y, m=rho)

	plt.show()


if __name__ == "__main__":
	main()