# EMcNeil_Demo_Material

This repository showcases a collection of advanced aerospace guidance, navigation, and control (GNC) projects developed by Erin McNeil as part of graduate studies in Aerospace Engineering at Georgia Tech. It includes simulations, estimation filters, trajectory optimization, and satellite control models implemented in Python, MATLAB, Simulink, and C++.

---

## 🚀 Project Highlights

### 🔭 Apollo 11 Final Approach: EKF vs UKF vs Particle Filter Trade Study
- **Tools**: Simulink, MATLAB
- **Description**: A lunar lander simulation with proximity sensor data fused using Extended Kalman Filter (EKF), Unscented Kalman Filter (UKF), and Particle Filter. Compares estimation accuracy near lunar surface during Apollo-style descent.
- **Focus**: Sensor fusion, optimal filtering, lunar navigation.

### 🛰️ Satellite Attitude Estimation with MEKF
- **Tools**: MATLAB
- **Description**: Simulates a satellite’s 3-axis attitude estimation using gyroscopes and a Multiplicative Extended Kalman Filter (MEKF).
- **Focus**: Attitude kinematics, estimation accuracy, optimal filtering.

### 🔧 Satellite Attitude PD Controller with Reaction Wheels
- **Tools**: MATLAB
- **Description**: Implements a PD control system to command and stabilize satellite attitude.
- **Focus**: Reaction wheel control, stability tuning, attitude tracking.

### 🔧 Satellite Sliding Mode Controller
- **Tools**: MATLAB
- **Description**: Implements a sliding-mode control system to command and stabilize satellite attitude.
- **Focus**: Reaction wheel control, stability tuning, attitude tracking.

### 🧲 Magnetorquer-Based De-spin Controller (LEO)
- **Tools**: MATLAB
- **Description**: Simulates Earth’s magnetic field and uses a magnetorquer to despin a satellite in LEO.
- **Focus**: B-dot controller, magnetic moment interaction, de-spin strategy.

### ⚙️ Satellite Architecture – Python & C++ 
- **Tools**: Python, C++, CMake
- **Description**: Ongoing port of legacy MATLAB-based spacecraft control architecture into modern C++/Python for embedded or simulation applications.
- **Focus**: Modular control design, simulation infrastructure, future Python toolkit.

### ✈️ Supersonic Fixed-Wing Aircraft Trajectory Optimization
- **Tools**: Python (NumPy, SciPy, Matplotlib)
- **Description**: Implements full-space and reduced-space trajectory optimization using surrogate models for aerodynamics, thrust, and atmosphere. Features adjoint gradient computation and complex-step Jacobian verification.
- **Focus**: Optimal control, minimum-time trajectories, scalable GNC architecture.

### 📡 Orbital Mechanics Simulations & Analysis Tools

MATLAB, Simulink, Python, and C++ projects demonstrating key mission and vehicle analyses:

* **Apollo 11 Final Approach**: EKF, UKF & Particle Filter sensor fusion
* **Satellite Attitude & Control**:

  * MEKF Estimation (3‑axis gyros)
  * PD Control with Reaction Wheels
  * Sliding‑Mode Control
  * Magnetorquer De‑spin (LEO)
* **Supersonic Fixed‑Wing Trajectory Optimization**
* **Two‑Body Orbit Propagation**: 3D RK4 integrator & Earth sphere plot
* **Earth–Moon Lagrange Points**: ΔV estimates & Jacobi analysis (L₂, L₄)
* **Earth–Mars Transfer**: Lambert solver & phase/synodic analysis
* **Lunar Orbit Capture**: Hyperbolic approach & circularization ΔV
* **MRO Aerobraking**: Three drag passes in Mars’ exponential atmosphere
* **Mars Atmospheric Skimming Flyby**: Turning-angle vs. altitude

*All tools under* `orbital_mechanics_simulations_and_analysis_tools/`\_

### ⚙️ Optimization Analysis Tools

Python and OpenMDAO examples showcasing optimization fundamentals:

* **Symbolic Differentiation**: Sympy gradients & Hessians (Himmelblau, Rosenbrock)
* **Solver Comparison**: BFGS, L‑BFGS‑B & trust‑constr paths on test functions
* **KKT Visualization**: Contours, constraint, and gradient arrows
* **Penalty Method Convergence**: Penalized → constrained solution as ρ↑
* **Trust‑Region Paths**: Solver trajectories from random starts
* **Log‑Barrier Method**: Symbolic barrier + L‑BFGS‑B convergence
* **Parametric KKT**: Solution vs. parameter β
* **AD vs. Complex‑Step**: JAX gradient vs. complex-step error
* **KS Minimization**: OpenMDAO ThermalAnalysis + KSAggregator with SLSQP
* **Projectile Optimization**: Maximize range in OpenMDAO
* **Warm‑Start Strategies**: Impact of initial guess on convergence

*All examples under* `optimization_analysis_tools/`\_


---

## 🔧 Technologies Used

- MATLAB / Simulink
- Python (NumPy, SciPy, Matplotlib)
- C++ with CMake
- Trajectory and dynamics simulation
- Orbital Mechanics: MATLAB, Simulink, Python (NumPy, SciPy), C++/CMake
- Optimization: Sympy, JAX, SciPy, OpenMDAO, Matplotlib
- Control & Estimation: Kalman Filters, control laws, adjoint methods

---

## 📂 Repository Structure

```
EMcNeil_Demo_Material/
│
├── Apollo_11_Final_Approach_EKF_UKF_ParticleFilter_TradeStudy/
├── Satellite_Architecture/
├── Satellite_Attitude_Optimal_Estimation_MEKFilter__3-axis_Gyroscopes/
├── Satellite_Attitude_PD_Controller__with_Reaction_Wheels/
├── Satellite_LEO_Magnetorquer_Despin_Controller/
└── Supersonic_Fixed-wing_Aircraft_Trajectory_Optimization/
├── orbital_mechanics_simulations_and_analysis_tools/
│   ├── Atmospheric_Drag_Pass_Mars/
│   ├── Earth-Moon_Lagrange_Points_DeltaV_Jacobi/
│   ├── Earth-to-Mars_Interplanetary_Transfer/
│   ├── Gauss_Problem_via_Universal_Variables/
│   ├── Lunar_Capture_and_Circularization_Maneuver/
│   ├── Two-Body_ODE_Integrator_3D_Orbit_Propagation/
│   └── Venus_Fly-by_Mission_Design/
│
├── optimization_analysis_tools/
│   ├── symbolic_diff_demo.py
│   ├── solver_comparison.py
│   ├── kkt_visualization.py
│   ├── penalty_convergence.py
│   ├── trust_region_paths.py
│   ├── log_barrier_demo.py
│   ├── parametric_kkt.py
│   ├── grad_cs_comparison.py
│   ├── projectile_optimization.py
│   └── warm_start_strategies.py
│
└── README_EMcNeil_Demo_Material.md
```

---

## 📄 License

This repository is shared for academic and research purposes. Licensing for specific components (e.g., Python modules or C++ libraries) may be added individually as development continues.

---

## 👤 Author

Erin McNeil  
M.S. Aerospace Engineering, Georgia Institute of Technology  
GitHub: [github.com/exs0087](https://github.com/exs0087)  
LinkedIn:(https://www.linkedin.com/in/ees0087/)

---

## 🛰️ Intended Audience

This repository is intended for:
- PhD advisors in aerospace engineering
- GNC researchers and engineers
- Aerospace software developers
- Graduate students studying estimation and control
