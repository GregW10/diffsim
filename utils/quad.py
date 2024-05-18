import scipy.integrate
import numpy as np
import sys

def func(x: float) -> float:
	return x**2 - 4

def integrate_riemann(f, a, b, num=1000000) -> float:
	ra = b - a
	dx = ra/num
	res = 0.0
	x = a
	i = 1
	while i <= num:
		res += f(x)*dx
		x = a + i*dx
		i += 1
	return res

def integrate_trapezoidal(f, a, b, num=1000000) -> float:
	ra = b - a
	dx = ra/num
	res = 0.0
	i = 0
	while i < num:
		x_lo = a + i*dx
		i += 1
		x_hi = a + i*dx
		f_lo = f(x_lo)
		f_hi = f(x_hi)
		diff = abs(f_hi - f_lo)
		mx = f_hi if f_hi >= f_lo else f_lo
		res += dx*(mx - 0.5*diff)
	return res


def main():
	f = func
	a = -2.0
	b = 2.0
	result = scipy.integrate.quad(f, a, b)
	print(f"Result from scipy = {result}")
	riemann = integrate_riemann(f, a, b, 1000000 if len(sys.argv) == 1 else int(sys.argv[1]))
	print(f"My Riemann result =     {riemann}")
	trap = integrate_trapezoidal(f, a, b, 1000000 if len(sys.argv) == 1 else int(sys.argv[1]))
	print(f"My trapezoidal result = {trap}")

if __name__ == "__main__":
	main()
