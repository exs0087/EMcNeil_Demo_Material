"""
Demo: JAX vs. Complex-Step Gradient Comparison

This script compares JAX automatic differentiation gradient of a projectile
range function with a complex-step approximation over varying step sizes,
and plots the absolute error vs. step size.

Usage:
    python grad_cs_comparison.py
"""

import numpy as np
import matplotlib.pyplot as plt
import jax
import jax.numpy as jnp

# Define projectile range function: range(theta) = v0^2/g * sin(2*theta)
def range_jax(theta):
    g = 9.81      # m/s^2
    v0 = 100.0    # m/s
    return (v0**2 / g) * jnp.sin(2 * theta)

# JAX gradient function
grad_jax_fn = jax.grad(range_jax)

# NumPy version for complex-step
def range_np(theta):
    g = 9.81
    v0 = 100.0
    return (v0**2 / g) * np.sin(2 * theta)

if __name__ == "__main__":
    # Point at which to evaluate gradient
    theta0 = np.pi / 6  # 30 degrees

    # Compute JAX gradient
    grad_jax_val = grad_jax_fn(theta0)

    # Range of complex-step sizes
    hs = np.logspace(-16, -1, 30)
    errors = []

    for h in hs:
        cs_grad = np.imag(range_np(theta0 + 1j * h)) / h
        errors.append(abs(cs_grad - grad_jax_val))

    # Plot error vs step size
    plt.figure()
    plt.loglog(hs, errors, marker='o', base=10)
    plt.xlabel('Complex-step size $h$')
    plt.ylabel('Absolute error |grad_cs - grad_jax|')
    plt.title('JAX vs. Complex-Step Gradient Error')
    plt.grid(True, which='both', ls='--')
    plt.savefig('grad_cs_comparison.png')
    print("JAX gradient at theta0:", grad_jax_val)
    print("Plot saved to grad_cs_comparison.png")
