"""
Warm-Start Strategies Demo

Demonstrates the effect of different initial guesses on convergence speed and solution quality
of SciPy's trust-constr solver on the Rosenbrock function.

Usage:
    python warm_start_strategies.py
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

def rosen(x):
    """Rosenbrock function."""
    return sum(100.0*(x[1:] - x[:-1]**2.0)**2.0 + (1 - x[:-1])**2.0)

def make_callback(history):
    def callback(xk, *args):
        history.append(xk.copy())
    return callback

def run_demo():
    n = 5
    guesses = {
        'zeros': np.zeros(n),
        'sinusoidal': np.sin(np.linspace(0, 2*np.pi, n)),
        'random': np.random.uniform(-1.5, 1.5, n)
    }

    iter_counts = {}
    final_objs = {}

    for label, x0 in guesses.items():
        history = []
        res = minimize(
            rosen,
            x0,
            method='trust-constr',
            callback=make_callback(history),
            options={'maxiter': 1000}
        )
        iter_counts[label] = len(history)
        final_objs[label] = res.fun
        print(f"{label}: iterations = {len(history)}, final objective = {res.fun:.6f}")

    # Plot iteration counts
    plt.figure()
    labels = list(iter_counts.keys())
    counts = [iter_counts[l] for l in labels]
    plt.bar(labels, counts)
    plt.ylabel('Iterations to Convergence')
    plt.title('Warm Start Strategies Convergence')
    plt.savefig('hw5_warm_start_iterations.png')
    print("Saved plot: hw5_warm_start_iterations.png")

    # Plot objective vs iterations for one example
    plt.figure()
    for label, x0 in guesses.items():
        history = []
        make_callback(history)
        # Rerun to capture history
        _ = minimize(
            rosen,
            x0,
            method='trust-constr',
            callback=make_callback(history),
            options={'maxiter': 1000}
        )
        objs = [rosen(x) for x in history]
        plt.plot(objs, label=label)
    plt.xlabel('Iteration')
    plt.ylabel('Objective Value')
    plt.title('Convergence of Objective')
    plt.legend()
    plt.savefig('warm_start_objective.png')
    print("Saved plot: warm_start_objective.png")

if __name__ == "__main__":
    run_demo()
