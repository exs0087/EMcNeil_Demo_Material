"""
Demo: Penalty Method Convergence

Shows how minimizing f(x) + rho/2 * (g(x))^2 for increasing rho
approaches the constrained solution of f(x) subject to g(x)=0.

Usage:
    python penalty_convergence.py
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

def f(x):
    return 2*x[0]**2 + x[1]**2 - 2*x[0] - 6*x[1]

def g(x):
    return x[0] + x[1] - 1

def penalty_obj(x, rho):
    return f(x) + 0.5 * rho * g(x)**2

def run_demo():
    rhos = [1, 10, 100, 1000, 10000]
    sols = []
    for rho in rhos:
        res = minimize(lambda x: penalty_obj(x, rho), np.array([0,0]))
        sols.append(res.x)
    sols = np.array(sols)

    plt.figure()
    plt.plot(rhos, sols[:,0], 'o-', label='x1')
    plt.plot(rhos, sols[:,1], 's-', label='x2')
    plt.xscale('log')
    plt.xlabel('rho'); plt.ylabel('solution components')
    plt.legend()
    plt.title('Penalty Convergence to Constrained Solution')
    plt.savefig('penalty_convergence.png')

if __name__ == "__main__":
    run_demo()
