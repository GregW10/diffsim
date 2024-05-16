import scipy.integrate as integrate
import numpy as np
import time


def func(x):
    return np.sin(10000*x)


def main():
    a = 0
    b = 3
    start = time.time_ns()
    res = integrate.quad(func, a, b)
    end = time.time_ns()
    print(f"Result = {res}\nTime taken = {end/1000000000} s")


if __name__ == "__main__":
    main()
