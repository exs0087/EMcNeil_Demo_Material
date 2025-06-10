"""
Demo: KKT Point Visualization

This script visualizes the solution to:
    minimize f(x) = 0.5 * x^T * Q * x + c^T * x
    subject to A x = b

It plots the quadratic contours, constraint lines, and marks the KKT point
with gradient arrows.

Usage:
    python kkt_visualization.py
"""

import numpy as np
import matplotlib.pyplot as plt

def objective(x):
    Q = np.array([[3, 0], [0, 1]])
    c = np.array([-1, -2])
    return 0.5 * x.T @ Q @ x + c.T @ x

def grad_obj(x):
    Q = np.array([[3, 0], [0, 1]])
    c = np.array([-1, -2])
    return Q @ x + c

# Constraint: x1 + x2 = 1
def constraint(x):
    return x[0] + x[1] - 1

def plot_demo():
    # Create grid
    xs = np.linspace(-1, 2, 400)
    ys = np.linspace(-1, 2, 400)
    X, Y = np.meshgrid(xs, ys)
    Z = np.vectorize(lambda a, b: objective(np.array([a,b])))(X, Y)

    # Plot contours
    plt.figure()
    cs = plt.contour(X, Y, Z, levels=30)
    plt.clabel(cs, inline=1, fontsize=8)

    # Plot constraint line
    x_line = np.linspace(-0.5, 1.5, 100)
    y_line = 1 - x_line
    plt.plot(x_line, y_line, 'r--', label='x1+x2=1')

    # Compute KKT solution analytically
    # Solve Q x + c + lambda * A^T = 0, A x = 1
    Q = np.array([[3,0],[0,1]])
    c = np.array([-1,-2])
    A = np.array([[1,1]])
    # KKT system
    KKT = np.block([[Q, A.T],[A, np.zeros((1,1))]])
    b = -np.concatenate([c, [1]])
    sol = np.linalg.solve(KKT, b)
    x_kkt = sol[:2]

    # Plot KKT point
    plt.plot(x_kkt[0], x_kkt[1], 'ko', label='KKT point')

    # Plot gradient arrow
    grad = grad_obj(x_kkt)
    plt.arrow(x_kkt[0], x_kkt[1], -grad[0]*0.1, -grad[1]*0.1,
              head_width=0.05, head_length=0.1, fc='k')

    plt.xlabel('x1'); plt.ylabel('x2')
    plt.legend()
    plt.title('KKT Visualization')
    plt.savefig('kkt_visualization.png')

if __name__ == "__main__":
    plot_demo()
