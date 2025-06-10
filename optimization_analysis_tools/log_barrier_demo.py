"""
Demo: Log-Barrier Method

Constructs a log-barrier for the constraint x1 + x2 - 1 = 0 and minimizes
the barrier-augmented objective for increasing barrier coefficients.

Usage:
    python log_barrier_demo.py
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import sympy as sp

# Define symbolic barrier function
x1, x2, t = sp.symbols('x1 x2 t')
f = 2*x1**2 + x2**2 - 2*x1 - 6*x2
g = x1 + x2 - 1
phi = f - t*sp.log(-g)

# Lambdify barrier objective and gradient
phi_num = sp.lambdify((x1,x2,t), phi, 'numpy')
grad_phi = sp.lambdify((x1,x2,t), [sp.diff(phi, x1), sp.diff(phi, x2)], 'numpy')

def barrier_obj(x, t):
    return phi_num(x[0], x[1], t)

def barrier_grad(x, t):
    return np.array(grad_phi(x[0], x[1], t), dtype=float)

def run_demo():
    ts = [1,10,100,1000]
    sol = []
    for t in ts:
        res = minimize(lambda x: barrier_obj(x,t), np.array([0.5,0.5]),
                       jac=lambda x: barrier_grad(x,t),
                       method='L-BFGS-B', bounds=[(None,None),(None,None)])
        sol.append(res.x)
    sol = np.array(sol)

    plt.figure()
    plt.plot(ts, sol[:,0], 'o-', label='x1')
    plt.plot(ts, sol[:,1], 's-', label='x2')
    plt.xscale('log')
    plt.legend(); plt.xlabel('Barrier t'); plt.title('Log-Barrier Convergence')
    plt.savefig('log_barrier_demo.png')

if __name__ == "__main__":
    run_demo()
