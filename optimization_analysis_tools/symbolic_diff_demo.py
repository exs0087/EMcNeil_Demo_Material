"""
Demo: Symbolic Gradient & Hessian for Test Functions

This script uses Sympy to derive and lambdify the gradient and Hessian
of Himmelblau's and Rosenbrock functions. It validates these symbolic
derivatives against finite-difference approximations at sample points.

Usage:
    python symbolic_diff_demo.py
"""

import numpy as np
import sympy as sp

def finite_diff_grad(f, x, h=1e-6):
    grad = np.zeros_like(x, dtype=float)
    for i in range(len(x)):
        xp = x.copy(); xm = x.copy()
        xp[i] += h; xm[i] -= h
        grad[i] = (f(xp) - f(xm)) / (2*h)
    return grad

# Define symbols
x1, x2 = sp.symbols('x1 x2')
# Define test functions
f_himmel = (x1**2 + x2 - 11)**2 + (x1 + x2**2 - 7)**2
f_rosen  = 100*(x2 - x1**2)**2 + (1 - x1)**2

# Symbolic derivatives
grad_himmel = [sp.diff(f_himmel, v) for v in (x1, x2)]
hess_himmel = sp.hessian(f_himmel, (x1, x2))
grad_rosen  = [sp.diff(f_rosen, v) for v in (x1, x2)]
hess_rosen  = sp.hessian(f_rosen, (x1, x2))

# Lambdify
f_himmel_num     = sp.lambdify((x1, x2), f_himmel, 'numpy')
grad_himmel_num  = sp.lambdify((x1, x2), grad_himmel, 'numpy')
hess_himmel_num  = sp.lambdify((x1, x2), hess_himmel, 'numpy')
f_rosen_num      = sp.lambdify((x1, x2), f_rosen, 'numpy')
grad_rosen_num   = sp.lambdify((x1, x2), grad_rosen, 'numpy')
hess_rosen_num   = sp.lambdify((x1, x2), hess_rosen, 'numpy')

if __name__ == "__main__":
    pts = [np.array([3.0, 2.0]), np.array([1.2, 1.2])]
    for f_num, grad_num, hess_num, name in [
        (f_himmel_num, grad_himmel_num, hess_himmel_num, "Himmelblau"),
        (f_rosen_num,  grad_rosen_num,  hess_rosen_num,  "Rosenbrock")
    ]:
        for x in pts:
            g_sym = np.array(grad_num(*x))
            H_sym = np.array(hess_num(*x))
            g_fd  = finite_diff_grad(lambda y: f_num(*y), x)
            print(f"--- {name} at {x} ---")
            print("Symbolic grad:    ", g_sym)
            print("Finite-diff grad: ", g_fd)
            print("Symbolic Hessian:", H_sym, "\n")
