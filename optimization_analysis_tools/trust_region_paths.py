"""
Demo: Trust-Region Solver Paths

Runs SciPy trust-constr solver from random starts on the constrained quadratic:
    minimize f(x)
    subject to x1 + x2 = 1

Plots the solver trajectories on objective contours.

Usage:
    python trust_region_paths.py
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

def f(x):
    return 2*x[0]**2 + x[1]**2 - 2*x[0] - 6*x[1]

cons = {'type':'eq', 'fun': lambda x: x[0]+x[1]-1}

def run_demo():
    xs = np.linspace(-1,2,400)
    ys = np.linspace(-1,2,400)
    X, Y = np.meshgrid(xs, ys)
    Z = np.vectorize(lambda a,b: f([a,b]))(X, Y)

    plt.figure()
    plt.contour(X,Y,Z,30)

    for _ in range(5):
        x0 = np.random.uniform(-0.5,1.5,2)
        traj = []
        def cb(xk, *args):
            traj.append(xk.copy())
        res = minimize(f, x0, constraints=cons, method='trust-constr', callback=cb)
        traj = np.array(traj)
        plt.plot(traj[:,0], traj[:,1], '-o')

    plt.title('Trust-Region Paths')
    plt.savefig('trust_region_paths.png')

if __name__ == "__main__":
    run_demo()
