##### Author: Erin McNeil        #####
##### Date: 04-26-2025           #####
#------------------------------------#
##### AE 6310 Final Project      #####
##### Full-Space Method          #####
#------------------------------------#

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from scipy.spatial import Delaunay
from numpy.linalg import lstsq
from scipy.optimize import minimize, NonlinearConstraint

#########################################
# SURROGATE MODEL BUILDING FUNCTIONS
#########################################

class CustomLinear2DInterpolator:
    def __init__(self, points, values):
        self.points = np.asarray(points)
        self.values = np.asarray(values)
        self.tri = Delaunay(self.points)
        self.coefficients = self._precompute_plane_coefficients()

    def _precompute_plane_coefficients(self):
        coeffs = []
        for simplex in self.tri.simplices:
            pts = self.points[simplex]  # shape (3, 2)
            vals = self.values[simplex]  # shape (3,)
            A = np.column_stack((pts, np.ones(3)))  # shape (3, 3)
            sol, *_ = lstsq(A, vals, rcond=None)  # fit plane: ax + by + c
            coeffs.append(sol)
        return np.array(coeffs)

    def __call__(self, x):
        x = np.atleast_2d(x)
        results = np.zeros(x.shape[0], dtype=np.complex128)
        for i, xi in enumerate(x):
            xi_real = np.real(xi)
            simplex = self.tri.find_simplex(xi_real)
            if simplex == -1:
                results[i] = np.nan
                continue
            a, b, c = self.coefficients[simplex]
            f_real = a * xi_real[0] + b * xi_real[1] + c
            if np.iscomplexobj(xi):
                grad = np.array([a, b])
                f_imag = np.dot(grad, np.imag(xi))
                results[i] = f_real + 1j * f_imag
            else:
                results[i] = f_real
        return results[0] if results.shape[0] == 1 else results


class ThrustSurrogateCS:
    def __init__(self, mach_data, alt_data, thrust_table):
        self.points, self.values = self._process_data(mach_data, alt_data, thrust_table)
        self.interp = CustomLinear2DInterpolator(self.points, self.values)
        self.plot_surface()  # Automatically generate plot after building the model

    def _process_data(self, mach_data, alt_data, thrust_table):
        points = []
        values = []
        for i, mach in enumerate(mach_data):
            for j, alt in enumerate(alt_data):
                val = thrust_table[i, j]
                if not np.isnan(val):
                    points.append([mach, alt])
                    values.append(val * 1000.0)  # Convert to lbf
        return np.array(points), np.array(values)

    def predict(self, x):
        x = np.atleast_2d(x)
        results = np.zeros(x.shape[0], dtype=np.complex128)

        for i, xi in enumerate(x):
            val = self.interp(xi)

            if np.isnan(val):
                # Fallback: nearest neighbor (real part only)
                xi_real = np.real(xi)
                distances = np.linalg.norm(self.points - xi_real, axis=1)
                val = self.values[np.argmin(distances)]

            results[i] = val

        return results[0] if results.size == 1 else results

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    def plot_surface(self, resolution=50):
        """
        Automatically plot a 3D surface of the thrust model upon creation.
        """
        mach_range = (self.points[:, 0].min(), self.points[:, 0].max())
        alt_range = (self.points[:, 1].min(), self.points[:, 1].max())

        mach_vals = np.linspace(*mach_range, resolution)
        alt_vals = np.linspace(*alt_range, resolution)
        Mach_grid, Alt_grid = np.meshgrid(mach_vals, alt_vals)
        Thrust_grid = np.zeros_like(Mach_grid)

        for i in range(resolution):
            for j in range(resolution):
                val = self.predict([[Mach_grid[i, j], Alt_grid[i, j]]])
                Thrust_grid[i, j] = np.real(val.item() if hasattr(val, "item") else val[0])

        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111, projection='3d')
        surf = ax.plot_surface(Mach_grid, Alt_grid, Thrust_grid,
                               cmap='viridis', edgecolor='k', linewidth=0.3, alpha=0.9)
        ax.set_xlabel("Mach")
        ax.set_ylabel("Altitude (kft)")
        ax.set_zlabel("Thrust (lbf)")
        ax.set_title("Thrust Surface from Surrogate Model")
        plt.tight_layout()
        plt.show()

def test_interp_gradient_raw(thrust_model):
    print("\n=== Raw Interpolator Complex-Step Test ===")

    x_test = np.array([1.2, 25.0])
    p = np.array([1.0, 1.0]) / np.sqrt(2)
    h = 1e-30

    x_pert = x_test + 1j * h * p
    f_pert = thrust_model.predict(x_pert)

    print("f_pert:", f_pert)
    print("Directional derivative (CS):", f_pert.imag / h)

def plot_interp_sensitivity(thrust_model, alt_kft=25.0):
    mach_vals = np.linspace(0.1, 2.5, 200)
    thrust_vals = []

    for mach in mach_vals:
        f = thrust_model.predict([[mach, alt_kft]])
        thrust_vals.append(np.real(f))

    plt.plot(mach_vals, thrust_vals, label=f'Altitude = {alt_kft} kft')
    plt.xlabel("Mach")
    plt.ylabel("Thrust (lbf)")
    plt.title("Thrust vs Mach (from Interpolator)")
    plt.grid(True)
    plt.legend()
    plt.show()


class CS1DCubicInterpolator:
    def __init__(self, x, y):
        self.x = np.array(x, dtype=float)
        self.y = np.array(y, dtype=float)
        self.n = len(x)
        self._compute_coefficients()

    def _compute_coefficients(self):
        x, y, n = self.x, self.y, self.n
        h = np.diff(x)
        A = np.zeros((n, n))
        rhs = np.zeros(n)

        A[0, 0] = A[-1, -1] = 1

        for i in range(1, n - 1):
            A[i, i - 1] = h[i - 1]
            A[i, i] = 2 * (h[i - 1] + h[i])
            A[i, i + 1] = h[i]
            rhs[i] = 3 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1])

        c = np.linalg.solve(A, rhs)

        b = np.zeros(n - 1)
        d = np.zeros(n - 1)
        for i in range(n - 1):
            b[i] = (y[i + 1] - y[i]) / h[i] - h[i] * (2 * c[i] + c[i + 1]) / 3
            d[i] = (c[i + 1] - c[i]) / (3 * h[i])

        self.a = y[:-1]
        self.b = b
        self.c = c[:-1]
        self.d = d
        self.h = h

    def __call__(self, x_new):
        x_new = np.atleast_1d(x_new)
        result = np.zeros_like(x_new, dtype=np.complex128)
        for j, xj in enumerate(x_new):
            i = np.searchsorted(self.x, np.real(xj)) - 1
            i = np.clip(i, 0, self.n - 2)
            dx = xj - self.x[i]
            result[j] = (self.a[i] +
                         self.b[i] * dx +
                         self.c[i] * dx ** 2 +
                         self.d[i] * dx ** 3)
        return result if result.size > 1 else result[0]

    def plot(self, resolution=200):
        x_dense = np.linspace(self.x.min(), self.x.max(), resolution)
        y_dense = self.__call__(x_dense)
        plt.plot(self.x, self.y, "ro", label="Data")
        plt.plot(x_dense, np.real(y_dense), "b-", label="Spline")
        plt.title("Custom Cubic Spline Interpolation")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.legend()
        plt.grid(True)
        plt.show()


class IdentityScaler:
    def transform(self, X):
        return X

class FunctionSurrogate:
    def __init__(self, func):
        self.func = func

    def predict(self, X):
        if isinstance(X, pd.DataFrame):
            X = X.values
        return self.func(X.ravel())


def build_piecewise_surrogate_model(data: pd.DataFrame, target: str, breakpoint_index):
    # Extract predictors (X) and target (y)
    X = data.loc[:, data.columns != target].values
    y = data[target].values

    # For simplicity, assume one predictor; if more, choose the appropriate column.
    # Sort data by X:
    sort_idx = np.argsort(X.ravel())
    X = X[sort_idx].ravel()  # make X a 1D array
    y = y[sort_idx]

    # Define the breakpoint value as the x-value corresponding to the breakpoint_index.
    x_break = X[breakpoint_index]

    # For the "linear" segment, use the first breakpoint_index+1 points.
    x_linear = X[:breakpoint_index + 1]
    y_linear = y[:breakpoint_index + 1]

    # For the "spline" segment, use the remaining points (including the breakpoint point)
    x_spline = X[breakpoint_index:]
    y_spline = y[breakpoint_index:]

    # Build the cubic spline on the second segment.
    spline_model = CS1DCubicInterpolator(x_spline, y_spline)

    def surrogate(x_new):
        """
        Evaluate the piecewise surrogate while preserving complex perturbations.
        For an input value x_new (which may be complex), this function:
          - If real(x_new) is in the linear segment (<= x_break), performs linear interpolation
            for the real part and estimates the imaginary part using a local derivative.
          - If real(x_new) is above x_break, uses the cubic spline which handles complex inputs.
        """
        x_new = np.atleast_1d(x_new)
        y_pred = np.empty_like(x_new, dtype=np.complex128)
        for idx, xi in enumerate(x_new):
            xi_real = np.real(xi)
            if xi_real <= x_break:
                # Compute the real part via linear interpolation.
                y_real = np.interp(xi_real, x_linear, y_linear)
                # Estimate a local slope using a finite difference over the linear data.
                i = np.searchsorted(x_linear, xi_real)
                if i == 0:
                    i = 0
                elif i >= len(x_linear):
                    i = len(x_linear) - 2
                else:
                    i = i - 1
                slope = (y_linear[i + 1] - y_linear[i]) / (x_linear[i + 1] - x_linear[i])
                y_imag = slope * np.imag(xi)
                y_pred[idx] = y_real + 1j * y_imag
            else:
                # For the spline segment, evaluate using the cubic spline directly.
                y_pred[idx] = spline_model(xi)
        return y_pred[0] if y_pred.size == 1 else y_pred

    # Evaluate the surrogate on the training inputs for diagnostics.
    yhat = surrogate(X)
    mse = np.mean((np.real(y) - np.real(yhat))**2)  # compare only real parts since training data is real
    def r2_score_manual(y, yhat):
        y = np.real(y)
        yhat = np.real(yhat)
        ss_res = np.sum((y - yhat) ** 2)
        ss_tot = np.sum((y - np.mean(y)) ** 2)
        return 1 - ss_res / ss_tot
    r2 = r2_score_manual(y, np.real(yhat))
    results = {'mse': mse, 'r2': r2}
    print(f"Piecewise surrogate for {target}: MSE = {mse:.4f}, R2 = {r2:.4f}")

    # Optionally, plot the surrogate versus the data.
    plt.figure(figsize=(8, 5))
    plt.scatter(X, y, c='red', marker='o', label='Data Points')
    x_dense = np.linspace(X.min(), X.max(), 200)
    plt.plot(x_dense, np.real(surrogate(x_dense)), 'b-', lw=2, label='Piecewise Model')
    plt.axvline(x=x_break, color='gray', linestyle='--', label=f'Breakpoint (x = {x_break:.2f})')
    plt.xlabel('Predictor (e.g. Mach)')
    plt.ylabel(target)
    plt.title(f'Piecewise Surrogate for {target}')
    plt.legend()
    plt.grid(True)
    plt.show()

    wrapped_model = FunctionSurrogate(surrogate)
    return wrapped_model, IdentityScaler(), results

class SplineSurrogate:
    def __init__(self, spline):
        self.spline = spline

    def predict(self, X):
        X = np.atleast_1d(X)
        X_flat = X.ravel()

        X_real = np.real(X_flat).astype(np.float64)
        X_imag = np.imag(X_flat)

        f_real = self.spline(X_real)

        if np.any(X_imag != 0):
            h = 1e-8
            f_pert = self.spline(X_real + h)
            df_dx = (f_pert - f_real) / h
            return f_real + 1j * df_dx * X_imag
        return f_real


def build_spline_surrogate_model(data: pd.DataFrame, target: str):
    """
    Build a surrogate model using a cubic spline interpolant for a given target.

    Args:
        data: pandas.DataFrame containing the independent variable(s) and target.
        target: str, the name of the target variable column.

    Returns:
        model: A SplineSurrogate object with a .predict() method.
        scaler: An IdentityScaler (kept for compatibility).
        results: A dictionary with training MSE and R² scores.
    """
    # Extract the independent variable(s) and target.
    # This function assumes there is one independent variable.
    X = data.loc[:, data.columns != target].values
    y = data[target].values

    # Sort the data by the independent variable (assuming one-dimension).
    sort_idx = np.argsort(X.ravel())
    X_sorted = X[sort_idx].ravel()
    y_sorted = y[sort_idx]

    # Build the cubic spline interpolant.
    spline_model = CS1DCubicInterpolator(X_sorted, y_sorted)

    model = SplineSurrogate(spline_model)
    scaler = IdentityScaler()

    # Evaluate the spline on the training points.
    yhat = model.predict(X_sorted)
    def mean_squared_error_manual(y, yhat):
        return np.mean((np.real(y) - np.real(yhat)) ** 2)
    mse = mean_squared_error_manual(y_sorted, np.real(yhat))
    def r2_score_manual(y, yhat):
        y = np.real(y)
        yhat = np.real(yhat)
        ss_res = np.sum((y - yhat) ** 2)
        ss_tot = np.sum((y - np.mean(y)) ** 2)
        return 1 - ss_res / ss_tot
    r2 = r2_score_manual(y_sorted, np.real(yhat))

    results = {'mse': mse, 'r2': r2}

    print(f"Cubic Spline surrogate for {target}:\n  MSE = {mse:.4f}, R2 = {r2:.4f}")

    # Optionally, plot the surrogate fit versus the data.
    plt.figure(figsize=(8, 5))
    plt.scatter(X_sorted, y_sorted, c='red', label='Data Points')
    # Create a dense grid for a smooth curve.
    X_dense = np.linspace(X_sorted.min(), X_sorted.max(), 200)
    plt.plot(X_dense, model.predict(X_dense), 'b-', lw=2, label='Cubic Spline Fit')
    plt.xlabel(list(data.loc[:, data.columns != target].columns)[0])
    plt.ylabel(target)
    plt.title(f"Cubic Spline Surrogate for {target}")
    plt.legend()
    plt.grid(True)
    plt.show()

    return model, scaler, results


#########################################
# DATA PREPARATION FOR SURROGATES
#########################################

# 1. Thrust Data
Mach_data_thrust = np.array([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8])
Alt_kft = np.array([0, 5, 10, 15, 20, 25, 30, 40, 50, 70])  # in thousands of ft
T_table = np.array([
    [24.2, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
    [28.0, 24.6, 21.1, 18.1, 15.2, 12.8, 10.7, np.nan, np.nan, np.nan],
    [28.3, 25.2, 21.9, 18.7, 15.9, 13.4, 11.2, np.nan, np.nan, np.nan],
    [30.8, 27.2, 23.8, 20.5, 17.3, 14.7, 12.3, 8.1, 4.9, np.nan],
    [34.5, 30.3, 26.6, 23.2, 19.8, 16.8, 14.1, 9.4, 5.6, 1.1],
    [37.9, 34.3, 30.4, 26.8, 23.3, 19.8, 16.8, 11.2, 6.8, 1.4],
    [36.1, 38.0, 34.9, 31.3, 27.3, 23.6, 20.1, 13.4, 8.3, 1.7],
    [np.nan, 36.6, 38.5, 36.1, 31.6, 28.1, 24.2, 16.2, 10.0, 2.2],
    [np.nan, np.nan, np.nan, 38.7, 35.7, 32.0, 28.1, 19.3, 11.9, 2.9],
    [np.nan, np.nan, np.nan, np.nan, np.nan, 34.6, 31.1, 21.7, 13.3, 3.1]
])
thrust_records = []
for i in range(len(Mach_data_thrust)):
    for j in range(len(Alt_kft)):
        if not np.isnan(T_table[i, j]):
            thrust_records.append({'Mach': Mach_data_thrust[i], 'Alt_kft': Alt_kft[j], 'Thrust': T_table[i, j] * 1000.0})

df_thrust = pd.DataFrame(thrust_records)
model_thrust = ThrustSurrogateCS(Mach_data_thrust, Alt_kft, T_table)
# print("Thrust Data (first few rows):")
# print(df_thrust.head())

# 2. CL_alpha Data
Mach_data_coeffs = np.array([0, 0.4, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8])
C_Lalpha_data = np.array([3.44, 3.44, 3.44, 3.58, 4.44, 3.44, 3.01, 2.86, 2.44])
df_CL = pd.DataFrame({'Mach': Mach_data_coeffs,
                      'CL_alpha': C_Lalpha_data})
# print("\nCL_alpha Data:")
# print(df_CL.head())

# 3. C_D0 Data
C_D0_data = np.array([0.013, 0.013, 0.013, 0.014, 0.031, 0.041, 0.039, 0.036, 0.035])
df_CD = pd.DataFrame({'Mach': Mach_data_coeffs,
                      'C_D0': C_D0_data})
# print("\nC_D0 Data:")
# print(df_CD.head())

# 4. Density Data
alt_data_ft = np.array([0, 5000, 10000, 15000, 20000, 25000, 30000,
                        35000, 40000, 45000, 50000, 55000, 60000, 65000, 70000], dtype=float)
rho_data = np.array([0.00237717, 0.00204814, 0.00175528, 0.00149698, 0.00126732,
                     0.00105621, 0.00090928, 0.00081988, 0.00066011, 0.00051340,
                     0.00041350, 0.000337, 0.000274, 0.000223, 0.000182], dtype=float)
df_rho = pd.DataFrame({'Altitude_ft': alt_data_ft,
                       'Density': rho_data})
# print("\nDensity Data:")
# print(df_rho.head())

# 5. Speed of Sound Data
a_data = np.array([1116.45, 1097.14, 1080.04, 1063.62, 1046.87,
                  1029.82, 1013.28, 997.46, 975.59, 954.22,
                  933.03, 912.84, 893.27, 874.44, 857.41], dtype=float)
df_a = pd.DataFrame({'Altitude_ft': alt_data_ft,
                     'Speed_of_Sound': a_data})
# print("\nSpeed of Sound Data:")
# print(df_a.head())

# 6. Eta Data
eta_data = np.array([0.54, 0.54, 0.54, 0.75, 0.79, 0.78, 0.89, 0.93, 0.93])
df_eta = pd.DataFrame({'Mach': Mach_data_coeffs,
                       'Eta': eta_data})
# print("\nEta Data:")
# print(df_eta.head())

#########################################
# Build Surrogate Models and Assemble Dictionary
#########################################
surrogate_models = {}
model_thrust = ThrustSurrogateCS(Mach_data_thrust, Alt_kft, T_table)
surrogate_models['Thrust'] = (model_thrust, IdentityScaler())

model_CL_spline, scaler_CL_spline, _ = build_piecewise_surrogate_model(df_CL, target='CL_alpha', breakpoint_index=3)
surrogate_models['CL_alpha'] = (model_CL_spline, scaler_CL_spline)

model_CD_spline, scaler_CD_spline, _ = build_piecewise_surrogate_model(df_CD, target='C_D0', breakpoint_index=3)
surrogate_models['C_D0'] = (model_CD_spline, scaler_CD_spline)

model_rho_spline, scaler_rho_spline, results_rho = build_spline_surrogate_model(df_rho, target='Density')
surrogate_models['Density'] = (model_rho_spline, scaler_rho_spline)

model_a_spline, scaler_a_spline, results_a = build_spline_surrogate_model(df_a, target='Speed_of_Sound')
surrogate_models['Speed_of_Sound'] = (model_a_spline, scaler_a_spline)

model_eta_poly, scaler_eta_poly, _  = build_piecewise_surrogate_model(df_eta, target='Eta', breakpoint_index=2)
surrogate_models['Eta'] = (model_eta_poly, scaler_eta_poly)

# =====================================
# SCALING FACTORS FOR STATES
# =====================================
scaling_factors = {
  'v'    : 500.0,
  'gamma': 0.5,
  'h'    : 5e4,
  'r'    : 1e5,
  'm'    : 1e3
}

residual_scaling = np.array([
  20.0,    # velocity
  2.0,     # gamma
  200.0,   # altitude
  1e7,     # range
  1e4      # mass
])

def scale_states(q):
    q = np.asarray(q)
    return np.array([
        q[0] / scaling_factors['v'],
        q[1] / scaling_factors['gamma'],
        q[2] / scaling_factors['h'],
        q[3] / scaling_factors['r'],
        q[4] / scaling_factors['m']
    ])

def unscale_states(q_scaled):
    q_scaled = np.asarray(q_scaled)
    return np.array([
        q_scaled[0] * scaling_factors['v'],
        q_scaled[1] * scaling_factors['gamma'],
        q_scaled[2] * scaling_factors['h'],
        q_scaled[3] * scaling_factors['r'],
        q_scaled[4] * scaling_factors['m']
    ])

def get_fullspace_scaling_vector(nt, blk, scaling_factors):
    """
    Return a vector of length (nt+1)*blk + 1 that contains the scaling
    factors applied to each variable in the full decision vector.
    """
    scale_vec = np.ones((nt + 1) * blk + 1)
    for k in range(nt + 1):
        idx = k * blk
        scale_vec[idx + 1] = scaling_factors['v']
        scale_vec[idx + 2] = scaling_factors['gamma']
        scale_vec[idx + 3] = scaling_factors['h']
        scale_vec[idx + 4] = scaling_factors['r']
        scale_vec[idx + 5] = scaling_factors['m']
    # The final time (last element) is left as-is
    return scale_vec


#########################################
# DYNAMICS AND TRAJECTORY INTEGRATION USING THE TRAPEZOIDAL METHOD
#########################################

class SurrogateVehicleDynamics:
    """
    Vehicle dynamics using surrogate models.
    """
    def __init__(self, surrogate_models):
        self.surrogate_models = surrogate_models
        self.g = 32.174       # ft/s^2
        self.S = 530.0        # ft^2
        self.Isp = 1600.0     # sec
        self.mass0 = 1304.0   # slugs
        self.gamma0 = np.radians(1.7)     # rad
        self.v_const = 500.0  # ft/s #temp -- remove later

    def clamp_complex(self, z, lo, hi):
        """
        Clamp real part of complex-valued variable z between lo and hi, preserving imaginary part.
        """
        real = np.clip(np.real(z), lo, hi)
        imag = np.imag(z)
        return real + 1j * imag

    def getAtmDensity(self, h):
        h = self.clamp_complex(h, 0.1, 70000)
        model, scaler = self.surrogate_models['Density']
        X_new = pd.DataFrame({'Altitude_ft': [h]})
        val = model.predict(scaler.transform(X_new))
        return val.item() if hasattr(val, "item") else val[0]

    def getAtmSoundSpeed(self, h):
        h = self.clamp_complex(h, 0.1, 70000.0)
        model, scaler = self.surrogate_models['Speed_of_Sound']
        X_new = pd.DataFrame({'Altitude_ft': [h]})
        val = model.predict(scaler.transform(X_new))
        return val.item() if hasattr(val, "item") else val[0]

    def getThrust(self, h, v):
        h = self.clamp_complex(h, 0.1, 70000.0)
        v = self.clamp_complex(v, 0.1, 2000.0)
        alt_kft = h / 1000.0
        a = self.getAtmSoundSpeed(h)
        M = self.clamp_complex(v / a, 0.01, 1.8)
        model, scaler = self.surrogate_models['Thrust']
        X_new = pd.DataFrame({'Mach': [M], 'Alt_kft': [alt_kft]})
        val = model.predict(scaler.transform(X_new))
        return val.item() if hasattr(val, "item") else val[0]

    def getCL_alpha(self, v, h):
        h = self.clamp_complex(h, 0.1, 70000.0)
        v = self.clamp_complex(v, 0.1, 2000.0)
        a = self.getAtmSoundSpeed(h)
        M = self.clamp_complex(v / a, 0.01, 1.8)
        model, scaler = self.surrogate_models['CL_alpha']
        X_new = pd.DataFrame({'Mach': [M]})
        val = model.predict(scaler.transform(X_new))
        return val.item() if hasattr(val, "item") else val[0]

    def getCL(self, alpha, v, h):
        h = self.clamp_complex(h, 0.1, 70000.0)
        v = self.clamp_complex(v, 0.1, 2000.0)
        cl_alpha = self.getCL_alpha(v, h)
        return cl_alpha * alpha

    def getCD(self, alpha, v, h):
        h = self.clamp_complex(h, 0.1, 70000.0)
        v = self.clamp_complex(v, 0.1, 2000.0)
        a = self.getAtmSoundSpeed(h)
        M = self.clamp_complex(v / a, 0.01, 1.8)
        X_new = pd.DataFrame({'Mach': [M]})

        model_cd, scaler_cd = self.surrogate_models['C_D0']
        val = model_cd.predict(scaler_cd.transform(X_new))
        val = np.atleast_1d(val)
        base_CD0 = val[0]
        model_eta, scaler_eta = self.surrogate_models['Eta']
        eta_val = model_eta.predict(scaler_eta.transform(X_new))
        eta_val = np.atleast_1d(eta_val)
        eta_val = eta_val[0]
        cl_alpha = self.getCL_alpha(v, h)
        return base_CD0 + eta_val * (cl_alpha * alpha) ** 2

    def getEta(self, v, h):
        h = self.clamp_complex(h, 0.1, 70000.0)
        v = self.clamp_complex(v, 0.1, 2000.0)
        a = self.getAtmSoundSpeed(h)
        M = self.clamp_complex(v / a, 0.01, 1.8)
        model, scaler = self.surrogate_models['Eta']
        X_new = pd.DataFrame({'Mach': [M]})
        val = model.predict(scaler.transform(X_new))
        return val.item() if hasattr(val, "item") else val[0]

    def getInitConditions(self):
        v0 = 380.0
        gamma0 = np.radians(1.7)  # changed from degrees to radians
        h0 = 0.0
        r0 = 0.0
        m0 = self.mass0
        return [v0, gamma0, h0, r0, m0]

    def getTargetConditions(self):
        vf = 986.5
        gammaf = 0.0  # still in radians now
        hf = 65617.0
        return [vf, np.radians(gammaf), hf, None, None]

    def getTargetScale(self):
        # final‐state tolerances:
        return [
            1.0,  # v  (ft/s)
            1, # gamma (radians)
            1.0,  # h  (ft)
            1,  # r  (ft)
            1  # m  (slug)
        ]

    def getNumStates(self):
        return 5

    def getNumControlStates(self):
        return 1

    def compute_f(self, q, alpha):
        """
        Compute the dynamics f(q, alpha).
        q = [v, gamma, h, r, m]
        """
        v, gamma, h, r, m = q

        rho = self.getAtmDensity(h)  # slugs/ft^3
        a = self.getAtmSoundSpeed(h)  # ft/sec
        M = v / a
        CL = self.getCL(alpha, v, h)
        CD = self.getCD(alpha, v, h)
        qinf = 0.5 * rho * v ** 2  # slugs/(ft * sec^2)
        L = qinf * self.S * CL  # lbf
        D = qinf * self.S * CD # lbf
        T = self.getThrust(h, v)

        f0 = (T * np.cos(alpha) - D) / m - self.g * np.sin(gamma)
        f1 = ((T * np.sin(alpha) + L) / (m * v) - (self.g / v) * np.cos(gamma))
        f2 = v * np.sin(gamma)
        f3 = v * np.cos(gamma)
        f4 = -T / (self.g * self.Isp)

        return np.array([f0, f1, f2, f3, f4])

    def computeSystemResidual(self, qdot, q, α, res):

        f0, f1, f2, f3, f4 = self.compute_f(q, α)
        res[0] = qdot[0] - f0
        res[1] = qdot[1] - f1
        res[2] = qdot[2] - f2
        res[3] = qdot[3] - f3
        res[4] = qdot[4] - f4

    def computeSystemJacobian(self, α_coef, β_coef, qdot, q, α, J):

        dh = 1e-8
        q_c = q.astype(np.complex128)
        for j in range(5):
            dq = np.zeros(5, dtype=np.complex128)
            dq[j] = 1j * dh
            f_pert = np.array(self.compute_f(q_c + dq, α))
            J[:, j] = - (f_pert.imag / dh)

    def computeSystemControlJacobian(self, qdot, q, α, Jx):
        dh = 1e-8
        f_pert = np.array(self.compute_f(q, α + 1j*dh))
        Jx[:] = - (f_pert.imag / dh)


# =====================================
# 3) FULL-SPACE TRAJECTORY CLASS
# =====================================
class MinTimeFullSpace:
    def __init__(self, system, num_time_steps, t):
        self.sys = system
        self.nt = num_time_steps
        self.t = t
        self.n  = system.getNumStates()

    def evalObj(self, xfull):
        return xfull[-1]

    def evalObjGrad(self, xfull):
        g = np.zeros_like(xfull)
        g[-1] = 1.0
        return g

    def computeResidual(self, xfull):
        cs = np.iscomplexobj(xfull)
        tf = np.real_if_close(xfull[-1])
        self.t = np.linspace(0, tf, self.nt + 1)

        blk = self.n + 1
        X   = xfull[:-1].reshape(self.nt + 1, blk)
        alphas = X[:, 0]
        Qs      = X[:, 1:]

        R = []

        q0_phys = unscale_states(Qs[0])
        ic_res  = q0_phys - np.array(self.sys.getInitConditions())
        R.extend(ic_res / residual_scaling)

        # trapezoid dynamics (5 rows per interval)
        for k in range(1, self.nt + 1):
            dt      = self.t[k] - self.t[k-1]
            qk_phys    = unscale_states(Qs[k])
            qkm1_phys  = unscale_states(Qs[k-1])

            qdot_phys = (qk_phys - qkm1_phys) / dt
            qmid_phys = 0.5 * (qk_phys + qkm1_phys)
            amid      = 0.5 * (alphas[k] + alphas[k-1])

            f_mid = self.sys.compute_f(qmid_phys, amid)
            f_res = qdot_phys - f_mid

            R.extend(f_res / residual_scaling)

        # terminal constraints on v,γ,h (3 rows)
        qf_phys = unscale_states(Qs[-1])
        vf, gf, hf, *_ = self.sys.getTargetConditions()
        sv, sg, sh, *_ = self.sys.getTargetScale()
        R.append((qf_phys[0] - vf) / sv)
        R.append((qf_phys[1] - gf) / sg)
        R.append((qf_phys[2] - hf) / sh)

        R = np.array(R, dtype=np.complex128)
        return R.real if not cs else R

    def computeJacobian(self, xfull):
        cs = np.iscomplexobj(xfull)
        x0 = xfull.astype(np.complex128)
        R0 = self.computeResidual(x0)  # complex‐aware
        m  = R0.size
        N  = x0.size

        J = np.zeros((m, N), dtype=np.complex128)
        dh = 1e-8

        for j in range(N):
            xp = x0.copy()
            xp[j] += 1j * dh
            Rp = self.computeResidual(xp)
            J[:, j] = Rp.imag / dh

        J = J.real if not cs else J
        print(f"in computeJacobian J  = {J}")
        return J


    def visualize(self, alphas, Q, interactive=False):
        t = np.linspace(0, 1, Q.shape[0])
        labels = ['v (ft/s)', 'γ (deg)', 'h (ft)', 'r (ft)', 'm (slug)']
        fig, axs = plt.subplots(Q.shape[1] + 1, 1, figsize=(10, 15), sharex=True)

        for i in range(Q.shape[1]):
            vals = Q[:, i].copy()
            if i == 0: vals *= scaling_factors['v']
            if i == 1: vals = np.degrees(vals * scaling_factors['gamma'])
            if i == 2: vals *= scaling_factors['h']
            if i == 3: vals *= scaling_factors['r']
            if i == 4: vals *= scaling_factors['m']

            axs[i].plot(t, vals)
            axs[i].set_ylabel(labels[i])
            axs[i].grid(True)

        axs[-1].plot(t, np.degrees(alphas))
        axs[-1].set_ylabel('α (deg)')
        axs[-1].set_xlabel('Normalized Time')
        axs[-1].grid(True)

        plt.tight_layout()
        if interactive:
            plt.show()

# --------------------------- Test Setup Functions --------------------------------------------#

# Toggle warm-start behavior
USE_WARMSTART = False  # set False to always use initial-guess generators


def make_initial_guess(dynamics, nt, blk, mode='shaped'):
    """
    Build an initial guess vector [alpha, q_scaled] for each time step plus final time.
    """
    q0      = np.array(dynamics.getInitConditions())
    qf_phys = np.array([986.5, 0.0, 65617.0, 275000.0, 1174.0])
    tf0     = 370 # tof initial guess
    tau_grid = np.linspace(0, 1, nt+1)
    x0      = np.zeros((nt+1)*blk + 1)

    for k, tau in enumerate(tau_grid):
        # select alpha pattern
        if mode == 'shaped':
            alpha_k = np.radians(2.0 + 0.5 * np.sin(2.5 * np.pi * tau))  # smooth sin wave at 2
            # alpha_k = np.radians(0 + 0.5 * np.sin(2.5 * np.pi * tau))  # smooth sin wave at 0
            # alpha_k = -1.0 + 1.5*np.sin(np.pi*tau) + 0.5*np.sin(3*np.pi*tau)
        elif mode == 'flat0':        alpha_k = 0.0
        elif mode == 'flatneg2':     alpha_k = np.radians(-2.0)
        elif mode == 'random':       alpha_k = np.clip(np.radians(np.random.normal(0,3)), np.radians(-8), np.radians(8))
        else:                        raise ValueError(f"Unknown mode '{mode}'")

        # shape between q0 and qf_phys
        shape    = (np.sin(np.pi*tau/2))**2
        q_phys   = (1-shape)*q0 + shape*qf_phys
        q_scaled = scale_states(q_phys)

        idx       = k*blk
        x0[idx]   = alpha_k
        x0[idx+1: idx+1+len(q0)] = q_scaled

    x0[-1] = tf0
    return x0


def setup_solver(nt, dynamics, t):
    """
    Initialize solver, variable bounds, and constraint.
    """
    solver = MinTimeFullSpace(dynamics, nt, t)
    n      = dynamics.getNumStates()
    blk    = n + 1

    # alpha bounds in radians, states unbounded
    bnds = [(-20*np.pi/180, 20*np.pi/180)]*(nt+1) + \
           [(None,None)]*((nt+1)*n) + \
           [(100,450)]

    con = NonlinearConstraint(
        solver.computeResidual,
        lb=np.zeros(n + nt*n + 3),
        ub=np.zeros(n + nt*n + 3),
        jac=solver.computeJacobian
    )
    return solver, blk, bnds, con

def run_mode(mode, solver, blk, bnds, con, dynamics, nt):
    """
    Run one optimization mode and record histories.
    """
    history = {'alphas':[], 'Q':[], 'obj':[], 'res':[]}
    objective_history = []

    # Define output directory for plots
    plot_dir = f"C:\\Users\\ees00\\Final Project AE 6310\\plots_from_full_space_{mode}"
    os.makedirs(plot_dir, exist_ok=True)

    def callback(xk, state):
        print(f"in callback")
        X = xk[:-1].reshape(nt+1, blk)
        alphas_k = X[:,0].copy()
        Q_k      = X[:,1:].copy()
        history['alphas'].append(alphas_k)
        history['Q'].append(Q_k)
        history['obj'].append(solver.evalObj(xk))
        history['res'].append(np.linalg.norm(solver.computeResidual(xk)))

        # save for warm-start
        np.savez('warmstart_4-25.npz', alphas=alphas_k, Q=Q_k, tf=xk[-1])

        obj = solver.evalObj(xk)
        objective_history.append(obj)
        print(f"Iteration {len(objective_history)}: objective = {obj:.6f}")

        # build time grid from current final time
        tf_k = xk[-1]
        t = np.linspace(0, tf_k, nt + 1)

        u = alphas_k
        Q = Q_k

        fig, ax = plt.subplots(6, 1, figsize=(8, 12))

        # before your plotting block
        scales = np.array([
            scaling_factors['v'],
            scaling_factors['gamma'],
            scaling_factors['h'],
            scaling_factors['r'],
            scaling_factors['m']
        ])
        # un–scale the whole Q array at once
        Q_phys = Q * scales

        # alpha in radians - convert to degrees
        alpha_deg = np.degrees(u)

        # plot
        ax[0].plot(t, alpha_deg, label='α (deg)')
        ax[1].plot(t, Q_phys[:, 0], label='v (ft/s)')
        ax[2].plot(t, np.degrees(Q_phys[:, 1]), label='γ (deg)')
        ax[3].plot(t, Q_phys[:, 2], label='h (ft)')
        ax[4].plot(t, Q_phys[:, 3], label='r (ft)')
        ax[5].plot(t, Q_phys[:, 4], label='m (slug)')

        for a in ax:
            a.legend()
            a.set_xlim(t[0], t[-1])
            ydata = a.lines[0].get_ydata()
            if np.allclose(ydata, ydata[0]):
                a.set_ylim(ydata[0] - 1, ydata[0] + 1)
            a.grid(True)

        plt.tight_layout()
        plt.suptitle(f"Full-Space Trajectory at Iteration {len(objective_history)}", y=1.02)

        # Save the figure
        filename = os.path.join(plot_dir, f"trajectory_iter_{len(objective_history):03d}.png")
        plt.savefig(filename)
        plt.pause(0.01)
        plt.close(fig)

        # early stopping
        if len(history['obj'])>10 and abs(history['obj'][-1]-history['obj'][-2])<1e-6:
            print("Plateau detected.")

    print(f"\n=== Mode: {mode} ===")

    # choose x0: warmstart or initial-guess generator
    if USE_WARMSTART and os.path.exists('warmstart.npz'):
        data = np.load('warmstart_4-25.npz')
        alphas_ws = data['alphas']  # shape (nt+1,)
        Q_ws      = data['Q']       # shape (nt+1,n)
        tf_ws     = data['tf']
        X_ws = np.hstack([alphas_ws.reshape(-1,1), Q_ws])
        x0  = np.concatenate([X_ws.ravel(), [tf_ws]])
        print("Loaded warmstart_4-25.npz for initial x0.")
    else:
        x0 = make_initial_guess(dynamics, nt, blk, mode)

    # visualize initial guess or warm-start
    X0 = x0[:-1].reshape(nt+1, blk)
    solver.visualize(X0[:,0], X0[:,1:])
    plt.suptitle(f"Init Guess ({mode})");
    plt.show()

    # number of decision variables
    num_vars = x0.size

    print(f"stepping in to minimize:")

    res = minimize(
        solver.evalObj, x0,
        method='trust-constr',
        jac=solver.evalObjGrad,
        bounds=bnds,
        constraints=[con],
        callback=callback,
        options={
            'initial_tr_radius': 30.0,
            'maxiter': 20,
            'gtol': 1, #tolerence for the constraint violations
            'xtol': 1e-3,
            'barrier_tol': 1e-4,
            'verbose': 3
        }
    )

    # --- number of control parameters from the reduced‐space problem ---
    num_ctrl_pts = 37
    # -----------------------------------------------

    # pull α out of the full‐space solution
    Xf_full = res.x[:-1].reshape(nt + 1, blk)
    alphas_full = Xf_full[:, 0]  # α in radians
    tf_full = res.x[-1]  # final time

    # normalized time grids
    t_full_norm = np.linspace(0, 1, nt + 1)
    t_reduced_norm = np.linspace(0, 1, num_ctrl_pts)

    # interpolate α→ reduced grid
    alpha_reduced = np.interp(t_reduced_norm, t_full_norm, alphas_full)

    # build reduced‐space design vector [α₁…αₙ, tf]
    x0_reduced = np.concatenate([alpha_reduced, [tf_full]])

    # save for the reduced‐solver’s load_design_vector()
    np.save("x_opt_reduced_space_4-25.npy", x0_reduced)
    print(" Warm-start saved to x_opt_reduced_space_4-25.npy")


    df_log = pd.DataFrame({
        "Objective": history["obj"],
        "ResidualNorm": history["res"]
    })
    df_log.to_csv(f"fullspace_convergence_{mode}.csv", index=False)
    print(f" Saved optimizer convergence log to fullspace_convergence_{mode}.csv")

    return history, res


def plot_results(objective_history, res, solver, blk, dynamics, nt, mode):
    """
    Plot and print final and objective_history results.
    """
    tf_opt = res.x[-1]
    print(f"Final time ({mode}): {tf_opt:.2f}s")

    Xf = res.x[:-1].reshape(nt+1, blk)
    alphas_f, Qf = Xf[:,0], Xf[:,1:]

    # final state printout
    q_end = unscale_states(Qf[-1])
    v, γ, h, r, m = q_end
    print(f"v={v:.2f}, γ={np.degrees(np.real(γ)):.2f}°, h={h:.2f}, r={r:.2f}, m={m:.2f}")

    # plot trajectories & objective_history
    solver.visualize(alphas_f, Qf)

    fig = plt.gcf()
    fig.savefig(f"trajectory_final_{mode}.png")
    print(f"🖼 Final trajectory plot saved to trajectory_final_{mode}.png")

    plt.title(f"Optimal Trajectory ({mode})")
    plt.show()
    plt.figure();
    plt.plot(objective_history['obj'], 'o-');
    plt.title('Objective');
    plt.grid();
    plt.show()

    plt.figure();
    plt.semilogy(objective_history['res'], 'o-');
    plt.title('Optimizer Objective Function Residual Norm');
    plt.grid(); plt.show()

    # dynamic pressure & L/D
    t = np.linspace(0, tf_opt, nt+1)
    qinf, LD = [], []
    for k in range(nt+1):
        q_phys = unscale_states(Qf[k]); αk = alphas_f[k]
        v_, γ_, h_, _, _ = q_phys
        ρ = dynamics.getAtmDensity(h_)
        qinf.append(0.5*ρ*v_**2)
        CL = dynamics.getCL(αk, v_, h_)
        CD = dynamics.getCD(αk, v_, h_)
        LD.append(CL/CD)
    plt.figure(); plt.plot(t, qinf); plt.title('Dynamic Pressure'); plt.grid();
    plt.savefig(f"dynamic_pressure_{mode}.png")
    plt.figure(); plt.plot(t, LD); plt.title('L/D Ratio'); plt.grid();
    plt.savefig(f"LD_ratio_{mode}.png")

    def export_trajectory_to_csv(t, Q, u, filename="fullspace_trajectory.csv"):
        Q_phys = Q * np.array([
            scaling_factors['v'],
            scaling_factors['gamma'],
            scaling_factors['h'],
            scaling_factors['r'],
            scaling_factors['m']
        ])
        df = pd.DataFrame(Q_phys, columns=['v (ft/s)', 'gamma (rad)', 'h (ft)', 'r (ft)', 'm (slug)'])
        df['time (s)'] = t
        df['alpha (deg)'] = np.degrees(u)
        df.to_csv(filename, index=False)
        print(f" Trajectory saved to {filename}")

    export_trajectory_to_csv(t, Qf, alphas_f)

    def print_final_state_table(Q_final, labels):
        print("\n Final State Table:")
        for i, label in enumerate(labels):
            print(f"{label:<20}: {Q_final[i]:.4f}")

    labels = ['Velocity (ft/s)', 'Gamma (rad)', 'Altitude (ft)', 'Range (ft)', 'Mass (slug)']
    print_final_state_table(unscale_states(Qf[-1]), labels)


def test_jacobian_fd(solver, x0, p, hs=None, filename="jacobian_fd_errors.png"):
    # default sweep of step sizes
    if hs is None:
        hs = np.logspace(-8, -2, 7)

    # compute base residual once
    R0 = solver.computeResidual(x0)

    # compute analytic directional derivative once
    J = solver.computeJacobian(x0)

    print("\nColumn-wise relative errors (FD vs J[:,j]):")
    h_fd = hs[-1]
    for j in range(J.shape[1]):
        ej = np.zeros_like(x0);
        ej[j] = 1.0
        Rph = solver.computeResidual(x0 + h_fd * ej)
        dRfd = (Rph - R0) / h_fd
        rel = np.linalg.norm(dRfd - J[:, j]) / (np.linalg.norm(J[:, j]) + 1e-20)
        print(f"  col {j:3d} → rel_err = {rel:.2e}")

    Jp = J @ p
    norm_Jp = np.linalg.norm(Jp) + 1e-20

    rel_errors = []
    for h in hs:
        # one extra residual per h
        Rph = solver.computeResidual(x0 + h * p)
        dRdp_fd = (Rph - R0) / h

        err = np.linalg.norm(dRdp_fd - Jp)
        rel = err / norm_Jp
        rel_errors.append(rel)

    rel_errors = np.array(rel_errors)

    #  Plot and save
    plt.figure(figsize=(5, 4))
    plt.loglog(hs, rel_errors, 'o-', base=10)
    plt.xlabel("Finite‐diff step size $h$")
    plt.ylabel(r"Relative error $\|\Delta_{FD}-Jp\|/\|Jp\|$")
    plt.title("Jacobian FD vs Analytic at x0 = 0")
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(filename)
    plt.show()
    plt.close()

    return hs, rel_errors


# =====================================================
#     MAIN
# =====================================================
if __name__ == '__main__':
    # build once
    dynamics = SurrogateVehicleDynamics(surrogate_models)
    nt = 370
    n = dynamics.getNumStates()
    blk = n + 1
    N = (nt + 1) * blk + 1
    t = np.linspace(0, 2.0, N + 1)
    solver, blk, bnds, con = setup_solver(nt, dynamics, t)

    for mode in ['shaped','flat0','flatneg2','random']:
        hist, res = run_mode(mode, solver, blk, bnds, con, dynamics, nt)
        plot_results(hist, res, solver, blk, dynamics, nt, mode)




    # ====== PROBLEM SETUP for FD vs. CS Jacobian Check ======
    h = 1e-30
    x_test = 0
    val = spline(x_test + 1j * h)
    print("Value at 0:", val)
    print("CS Derivative:", val.imag / h)

    x0 = make_initial_guess(dynamics, nt, blk, mode='shaped')
    p = np.random.randn(x0.size)
    p /= np.linalg.norm(p)

    hs, rel_errors = test_jacobian_fd(solver, x0, p)

    for h, e in zip(hs, rel_errors):
        print(f"h={h:.1e} → rel_error={e:.2e}")

