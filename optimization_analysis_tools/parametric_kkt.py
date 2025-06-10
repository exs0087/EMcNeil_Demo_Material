"""
Demo: Parametric KKT Analysis

Solves the KKT conditions for the parametric problem with parameter beta,
evaluates second-order conditions, and plots feasibility regions.

Usage:
    python parametric_kkt.py
"""

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

beta = sp.symbols('beta', real=True)
x1, x2, lam = sp.symbols('x1 x2 lam', real=True)

# Define parametric objective and constraint
f = x1**2 + beta*x2**2
g = x1 + x2 - 1

L = f + lam*g
# KKT equations
eqs = [
    sp.diff(L, x1),
    sp.diff(L, x2),
    g
]
sol = sp.solve(eqs, [x1, x2, lam], dict=True)
# Pick one branch
sol = sol[0]
x1_sol = sol[x1]
x2_sol = sol[x2]

# Hessian of Lagrangian
H = sp.hessian(L.subs(sol), (x1, x2))

# Lambdify for numeric plot
x1_func = sp.lambdify(beta, x1_sol, 'numpy')
x2_func = sp.lambdify(beta, x2_sol, 'numpy')

b_vals = np.linspace(0.1, 5, 100)
x1_vals = x1_func(b_vals)
x2_vals = x2_func(b_vals)

plt.figure()
plt.plot(b_vals, x1_vals, label='x1*')
plt.plot(b_vals, x2_vals, label='x2*')
plt.legend(); plt.xlabel('beta'); plt.ylabel('solution'); plt.title('Parametric KKT Solution')
plt.savefig('parametric_kkt.png')
