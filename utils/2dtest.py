import math
import scipy.integrate


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


def quad(func, a, b, tol=1e-9):
    def quad_rec(funcc, aa, bb, dxoo2, dxoo4, faa, fbb, toll):
        estimatee = dxoo2*(faa + fbb)
        mm = (aa + bb)/2
        fmm = funcc(mm)
        iress = dxoo4*(faa + 2*fmm + fbb)
        if abs(estimatee - iress) > toll:
            return quad_rec(funcc, aa, mm, dxoo4, dxoo4/2, faa, fmm, toll/2) + \
                   quad_rec(funcc, mm, bb, dxoo4, dxoo4/2, fmm, fbb, toll/2)
        return iress
    dxo2 = (b - a)/2.0
    dxo4 = dxo2/2.0
    fa = func(a)
    fb = func(b)
    estimate = dxo2*(fa + fb)
    m = (a + b)/2
    fm = func(m)
    ires = dxo4*(fa + 2*fm + fb)
    if abs(estimate - ires) > tol:
        return quad_rec(func, a, m, dxo4, dxo4/2, fa, fm, tol/2) + \
               quad_rec(func, m, b, dxo4, dxo4/2, fm, fb, tol/2)
    return ires


def dblquad(func, ya, yb, gfun, hfun, tol=1e-5):
    def dblquad_rec(funcc, yaa, ybb, dyoo2, dyoo4, fyaa, fybb, gfunn, hfunn, toll):
        ymm = (yaa + ybb)/2.0
        fymm = quad(lambda x: funcc(x, ymm), gfunn(ym), hfunn(ym), toll)
        estimatee = dyoo2*(fyaa + fybb)
        iress = dyoo4*(fyaa + 2*fymm + fybb)
        if abs(estimatee - iress) > toll:
            return dblquad(funcc, yaa, ymm, gfunn, hfunn, toll/2) + \
                   dblquad(funcc, ymm, ybb, gfunn, hfunn, toll/2)
        return iress
    dy = yb - ya
    dyo2 = dy/2.0
    dyo4 = dy/4.0
    ym = (ya + yb)/2.0
    fya = quad(lambda x: func(x, ya), gfun(ya), hfun(ya), tol)
    fym = quad(lambda x: func(x, ym), gfun(ym), hfun(ym), tol)
    fyb = quad(lambda x: func(x, yb), gfun(yb), hfun(yb), tol)
    estimate = dyo2*(fya + fyb)
    ires = dyo4*(fya + 2*fym + fyb)
    if abs(estimate - ires) > tol:
        return dblquad_rec(func, ya, ym, dyo4, dyo4/2, fya, fym, gfun, hfun, tol/2) + \
               dblquad_rec(func, ym, yb, dyo4, dyo4/2, fym, fyb, gfun, hfun, tol/2)
    return ires


def func1d(x):
    return x**4 + 5*x**3 - 4*x**2 + 6*x - 7


def func2d(x, y):
    return x**2*y**2


def gfunc(y):
    return 1 - math.sqrt(1 - (y - 1)**2)


def hfunc(y):
    return 1 + math.sqrt(1 - (y - 1)**2)


def main():
    xa = 0
    xb = 2
    ya = 0
    yb = 2
    numx = 5_000
    numy = 5_000
    res2d = dbl(func2d, xa, xb, ya, yb, numx, numy)
    print(f"Result 2D: {res2d}, real result: {scipy.integrate.dblquad(func2d, ya, yb, lambda y: ya, lambda y: yb)[0]}")
    res1d = quad(func1d, xa, xb)
    print(f"Result 1D: {res1d}, real result: {scipy.integrate.quad(func2d, ya, yb)[0]}")
    res2dq = dblquad(func2d, ya, yb, gfunc, hfunc)
    # print(f"Result 2D AQ: {res2dq}, real result: {scipy.integrate.dblquad(func2d, ya, yb, gfunc, hfunc)[0]}")


if __name__ == "__main__":
    main()
