import math


def dbl(func, xa, xb, ya, yb, numx=1_000_000, numy=1_000_000):
    x_range = xb - xa
    y_range = yb - ya
    dx = x_range/numx
    dy = y_range/numy
    dxdy = dx*dy
    x_start = xa + dx/2
    y_start = ya + dy/2
    res = 0
    x_counter = 0
    y_counter = 0
    while y_counter < numy:
        x_counter = 0
        y = y_start + y_counter*dy
        while x_counter < numx:
            x = x_start + x_counter*dx
            res += dxdy*func(x, y)
            x_counter += 1
        y_counter += 1
    return res


def dblquad(func, xa, xb, ya, yb):
    pass


def func(x, y):
    return x**2*y**2


def main():
    xa = 0
    xb = 4
    ya = 0
    yb = 4
    numx = 20_000
    numy = 20_000
    res = dbl(func, xa, xb, ya, yb, numx, numy)
    print(f"Result: {res}")


if __name__ == "__main__":
    main()
