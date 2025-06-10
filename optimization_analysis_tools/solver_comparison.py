"""
Demo: Solver Comparison on Himmelblau & Rosenbrock (Updated)

Runs SciPy solvers from random starts on two test functions
and plots their optimization paths on contour maps.

Usage:
    python solver_comparison.py
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# Test functions
def himmelblau(x):
    x1, x2 = x
    return (x1**2 + x2 - 11)**2 + (x1 + x2**2 - 7)**2

def rosenbrock(x):
    x1, x2 = x
    return 100*(x2 - x1**2)**2 + (1 - x1)**2

# Solvers to compare
methods = [
    ("BFGS",       {"method": "BFGS"}),
    ("L-BFGS-B",   {"method": "L-BFGS-B"}),
    ("Trust-Constr", {"method": "trust-constr"})
]

def run_demo(func, title):
    # Contour
    X = np.linspace(-3, 3, 400)
    Y = np.linspace(-3, 3, 400)
    XX, YY = np.meshgrid(X, Y)
    Z = np.vectorize(lambda a,b: func([a,b]))(XX, YY)

    plt.figure()
    plt.contour(XX, YY, Z, levels=50, cmap="coolwarm")
    starts = np.random.uniform(-2.5, 2.5, size=(5,2))

    for name, opts in methods:
        for x0 in starts:
            traj = []
            def callback(xk, *args):
                traj.append(xk.copy())
            res = minimize(func, x0, callback=callback, **opts)
            xs = np.array(traj)
            plt.plot(xs[:,0], xs[:,1], label=name, alpha=0.7)
    plt.title(f"{title} Optimization Paths")
    plt.legend()
    plt.savefig(f"{title.lower().replace(' ','_')}_paths.png")
    plt.close()

if __name__ == "__main__":
    run_demo(himmelblau, "Himmelblau")
    run_demo(rosenbrock, "Rosenbrock")
