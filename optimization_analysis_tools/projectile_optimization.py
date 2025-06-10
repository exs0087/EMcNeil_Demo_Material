"""
Projectile Optimization Demo (Fixed Promotions)

Defines a projectile motion component in OpenMDAO and finds the optimal launch
angle by minimizing the negative of the range. Uses model.add_objective with scaler
to perform maximization.

Usage:
    python projectile_optimization.py
"""

import numpy as np
import matplotlib.pyplot as plt
from openmdao.api import Problem, Group, IndepVarComp, ExplicitComponent, ScipyOptimizeDriver

class ProjectileComp(ExplicitComponent):
    def setup(self):
        self.add_input('angle', val=0.5, desc='Launch angle [rad]')
        self.add_output('range', val=0.0, desc='Projectile range [m]')
        self.declare_partials(of='range', wrt='angle', method='fd')
    def compute(self, inputs, outputs):
        g = 9.81      # m/s^2
        v0 = 100.0    # m/s
        theta = inputs['angle']
        outputs['range'] = (v0**2 / g) * np.sin(2 * theta)

if __name__ == "__main__":
    prob = Problem()
    model = prob.model = Group()

    # Independent variable for launch angle (output promotion)
    model.add_subsystem('angle_input',
                        IndepVarComp('angle', val=0.5, desc='Launch angle [rad]'),
                        promotes_outputs=['angle'])

    # Projectile component
    model.add_subsystem('proj',
                        ProjectileComp(),
                        promotes=['angle', 'range'])

    # Declare design variable and objective
    model.add_design_var('angle', lower=0.0, upper=np.pi/2)
    # Maximize range by minimizing negative range
    model.add_objective('range', scaler=-1.0)

    # Setup the optimizer driver
    prob.driver = ScipyOptimizeDriver()
    prob.driver.options['optimizer'] = 'SLSQP'
    prob.driver.options['maxiter'] = 200

    prob.setup()
    prob.run_driver()

    angle_opt = prob.get_val('angle')[0]
    range_opt = prob.get_val('range')[0]
    print(f"Optimal launch angle: {angle_opt:.4f} rad ({angle_opt*180/np.pi:.2f} deg)")
    print(f"Maximum range: {range_opt:.2f} m")

    # Plot range vs angle
    angles = np.linspace(0, np.pi/2, 100)
    ranges = (100**2/9.81) * np.sin(2 * angles)
    plt.figure()
    plt.plot(angles*180/np.pi, ranges, label='Range Curve')
    plt.axvline(angle_opt*180/np.pi, color='r', linestyle='--', label='Optimal')
    plt.xlabel('Launch Angle (deg)')
    plt.ylabel('Range (m)')
    plt.title('Projectile Range vs Launch Angle')
    plt.legend()
    plt.savefig('projectile_optimization.png')
    print("Plot saved to projectile_optimization.png")
