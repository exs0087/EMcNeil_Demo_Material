##### Author: Erin McNeil                   #####
##### Date: 04-26-2025                      #####
#-----------------------------------------------#
##### AE 6310 Final Project                 #####
##### Reduced Space with Adjoint Derivative #####
#-----------------------------------------------#

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

    # Wrap the surrogate function in a class with a .predict() method.
    wrapped_model = FunctionSurrogate(surrogate)
    return wrapped_model, IdentityScaler(), results

class SplineSurrogate:
    def __init__(self, spline):
        self.spline = spline

    def predict(self, X):
        # Ensure input is 1D numpy array (flattened if needed)
        X = np.atleast_1d(X)
        X_flat = X.ravel()

        # Extract real and imaginary parts
        X_real = np.real(X_flat).astype(np.float64)
        X_imag = np.imag(X_flat)

        # Evaluate spline on real part
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

    # Wrap the spline model in a class with a predict method.
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
        return [10.0, 1.0, 1e2, 1.0, 1.0]

    def getNumStates(self):
        return 5

    def getNumControlStates(self):
        return 1

    def safe_gamma(self, gamma, limit=np.radians(89.9)):
        gamma_real = np.real(gamma)
        gamma_clamped = np.clip(gamma_real, -limit, limit)
        return gamma_clamped + 1j * np.imag(gamma)

    def computeSystemResidual(self, qdot, q, x, res):
        v, gamma, h, r, m = q
        alpha = x

        # Safely clamp gamma (only real part)
        gamma = self.safe_gamma(gamma)

        # Clamp state and control values early to prevent surrogate model errors
        h = self.clamp_complex(h, 0.1, 70000.0)
        v = self.clamp_complex(v, 0.1, 2000.0)
        alpha = self.clamp_complex(alpha, np.radians(-20), np.radians(20))

        rho = self.getAtmDensity(h)  # slugs/ft^3
        a = self.getAtmSoundSpeed(h)  # ft/sec
        CL = self.getCL(alpha, v, h)
        CD = self.getCD(alpha, v, h)
        qinf = 0.5 * rho * v ** 2
        L = qinf * self.S * CL
        D = qinf * self.S * CD
        T = self.getThrust(h, v)

        f0 = (T * np.cos(alpha) - D) / m - self.g * np.sin(gamma)
        f1 = ((T * np.sin(alpha) + L) / (m * v) - (self.g / v) * np.cos(gamma))
        f2 = v * np.sin(gamma)
        f3 = v * np.cos(gamma)
        f4 = -T / (self.g * self.Isp)

        res[:] = qdot - np.array([f0, f1, f2, f3, f4], dtype=np.complex128)

    def computeSystemJacobian(self, alpha, beta, qdot, q, x, J):
        n = self.getNumStates()
        dh = 1e-30

        q_complex = q.astype(np.complex128)
        qdot_complex = qdot.astype(np.complex128)
        res = np.zeros(n, dtype=np.complex128)

        for i in range(n):
            e = np.zeros(n)
            e[i] = 1.0

            q_pert = q_complex + 1j * dh * alpha * e
            qdot_pert = qdot_complex + 1j * dh * beta * e

            self.computeSystemResidual(qdot_pert, q_pert, x, res)
            J[:, i] = np.imag(res) / dh

        return

    def computeSystemControlJacobian(self, qdot, q, x, Jx):
        n = self.getNumStates()
        dh = 1e-30

        qdot_complex = qdot.astype(np.complex128)
        q_complex = q.astype(np.complex128)
        x_pert = x + 1j * dh

        res = np.zeros(n, dtype=np.complex128)
        self.computeSystemResidual(qdot_complex, q_complex, x_pert, res)
        Jx[:] = np.imag(res) / dh

        return

class MinTimeTrajectory:
    def __init__(self, system, num_time_steps, num_ctrl_pts):
        self.system = system
        self.num_states = self.system.getNumStates()
        self.num_time_steps = num_time_steps
        self.num_ctrl_pts = num_ctrl_pts
        self.num_design_vars = self.num_ctrl_pts + 1
        self.ctrl_nodes = 0.5 - 0.5 * np.cos(np.linspace(0, np.pi, self.num_ctrl_pts))
        self.qinit = np.array(self.system.getInitConditions(), dtype=float)
        self.qtarget = self.system.getTargetConditions()
        self.num_constraints = sum(1 for qval in self.qtarget if qval is not None)
        self.max_newton_iters = 100
        self.newton_tol = 1e-6
        self.knots = 0.5 - 0.5 * np.cos(np.linspace(0, np.pi, self.num_ctrl_pts))
        self.fig = None
        self.ax = None


    def evalInterp(self, u):
        """Given a parametric point, u, evaluate the Lagrange basis"""
        N = np.ones(self.num_ctrl_pts)
        for i in range(self.num_ctrl_pts):
            for j in range(self.num_ctrl_pts):
                if i != j:
                    N[i] *= (u - self.knots[j]) / (self.knots[i] - self.knots[j])

        return N

    def getControl(self, u, x):
        """Return the control value given the parametric point u and x input"""
        return np.dot(self.evalInterp(u), x[:-1])

    def addControlDeriv(self, u, dfdctrl, dfdx):
        """Add the derivative of a function w.r.t. the control to the derivative w.r.t. x"""
        N = self.evalInterp(u)
        dfdx[:-1] += dfdctrl * N
        return

    def evalObj(self, x):
        fobj = x[-1]

        return fobj

    def evalObjGradient(self, x):
        """Evaluate the objective gradient"""
        g = np.zeros(x.shape, dtype=x.dtype)
        g[-1] = 1.0

        return g

    def evalCon(self, x):
        """Evaluate the constraints"""

        # Compute the full trajectory based on the input
        self.q = self.computeTrajectory(x)

        # Compute the constraints
        con = np.zeros(self.num_constraints, dtype=x.dtype)

        count = 0
        scale = self.system.getTargetScale()
        for state, qval in enumerate(self.qtarget):
            if qval is not None:
                con[count] = (self.q[-1, state] - qval) / scale[state]
                count += 1

        print('Objective = ', x[-1], 'Constraints = ', con)
        return con

    def evalConGradient(self, x, interactive=True):
        """Evaluate the constraint gradient"""

        # Compute the full trajectory
        q = self.computeTrajectory(x)

        # Create a vector to store the derivative of the constraints
        dcdx = np.zeros((self.num_constraints, self.num_design_vars), dtype=x.dtype)

        count = 0
        scale = self.system.getTargetScale()
        for state, qval in enumerate(self.qtarget):
            if qval is not None:
                self.computeAdjointDeriv(q, x, state, dcdx[count, :])
                dcdx[count, :] /= scale[state]
                count += 1

        self.visualize(x, q, interactive=interactive)

        return dcdx

    def computeTrajectory(self, x):
        """
        Given the design variable vector, compute the state history.

        x = [alpha_0, alpha_1, ..., alpha_{num_ctrl_pts-1}, tf]

        Args:
            x: The design variable vector

        Returns:
            q: The state variables
        """

        # Set the time steps
        self.t = np.linspace(0, x[-1], self.num_time_steps + 1)

        # Allocate space for the state variables
        q = np.zeros((self.num_time_steps + 1, self.num_states), dtype=x.dtype)

        # Set the initial conditions.
        q[0, :] = self.qinit

        # Compute the residual and Jacobian
        res = np.zeros(self.num_states, dtype=q.dtype)
        J = np.zeros((self.num_states, self.num_states), dtype=q.dtype)

        # Integrate forward in time
        for k in range(1, self.num_time_steps + 1):
            # Copy the starting point for the first iteration
            q[k, :] = q[k-1, :]

            residuals_list = []
            # Solve the nonlinear equations for q[k]
            for j in range(self.max_newton_iters):
                # Compute the state values at the mid-point
                alpha = 0.5
                qk = alpha*(q[k, :] + q[k-1, :])

                # Compute the time derivative approximation
                beta = 1.0/(self.t[k] - self.t[k-1])
                qkdot = beta*(q[k, :] - q[k-1, :])

                # Compute the control value at the mid-point
                xk = self.getControl((k - 0.5)/self.num_time_steps, x)

                # Compute the residuals
                self.system.computeSystemResidual(qkdot, qk, xk, res)

                # Compute the system Jacobian matrix
                self.system.computeSystemJacobian(alpha, beta, qkdot, qk, xk, J)

                # Solve for the update using Newton's method
                try:
                    update = np.linalg.solve(J, res)
                except:
                    print('Jacobian factorization failed')
                    print(J)
                    exit(0)
                q[k, :] -= update

                if not np.isfinite(q[k, :]).all():
                    print(f" Non-finite state at step {k}: {q[k, :]}")
                    raise RuntimeError("Trajectory integration diverged.")

                # Check for convergence
                rnorm = np.sqrt(np.dot(res, res))
                residuals_list.append(rnorm)
                if rnorm < self.newton_tol:
                    break
            else:
                rnorm = np.sqrt(np.dot(res, res))

        return q

    def computeAdjointDeriv(self, q, x, state, dfdx):
        """
        Compute the derivative of the final state with the specified state
        index with respect to the control.
        """

        # Set the time steps
        self.t = np.linspace(0, x[-1], self.num_time_steps + 1)

        # Zero-out the contributions to the state variables
        dfdx[:] = 0.0

        # Set the right-hand-side for the adjoint equations
        res = np.zeros(self.num_states, dtype=dfdx.dtype)
        res[state] = 1.0 # df/du

        # The Jacobian matrix
        J = np.zeros((self.num_states, self.num_states), dtype=dfdx.dtype)

        # The Jacobian matrix w.r.t. the control
        Jx = np.zeros(self.num_states, dtype=dfdx.dtype)

        # Integrate the adjoint in reverse
        for k in range(self.num_time_steps, 0, -1):
            # Compute the state values at the mid-point
            alpha = 0.5
            qk = alpha*(q[k, :] + q[k-1, :])

            # Compute the time derivative approximation
            beta = 1.0/(self.t[k] - self.t[k-1])
            qkdot = beta*(q[k, :] - q[k-1, :])

            # Compute the control value
            xk = self.getControl((k - 0.5)/self.num_time_steps, x)

            # Compute the Jacobian matrix
            self.system.computeSystemJacobian(alpha, beta, qkdot, qk, xk, J)

            # Compute the adjoint variables
            adjoint = -np.linalg.solve(J.T, res)

            # Compute the control input Jacobian
            self.system.computeSystemControlJacobian(qkdot, qk, xk, Jx)

            # Add the contribution from the total derivative
            dfdctrl = np.dot(Jx, adjoint)

            # Add the control derivative
            self.addControlDeriv((k - 0.5)/self.num_time_steps, dfdctrl, dfdx)

            # Add the contributions from the time derivative
            dfdx[-1] -= beta * (1.0 / self.num_time_steps) * np.dot(qkdot, adjoint)

            # Compute the right-hand-side for the next adjoint
            self.system.computeSystemJacobian(alpha, -beta, qkdot, qk, xk, J)

            # Update the right-hand-side for the adjoint
            res = np.dot(J.T, adjoint)

        return

    def getControlTrajectory(self, x):
        """Compute the control trajectory for visualization"""
        u = np.zeros(self.num_time_steps + 1)

        for k in range(self.num_time_steps + 1):
            u[k] = self.getControl(k / self.num_time_steps, x)

        return u

    def visualize(self, x, q, interactive=True):
        """
        Visualize the result.

        During an optimization, use interactive=True so that the
        plot will update as the optimization problem is solved.
        """

        if self.fig is None:
            self.fig, self.ax = plt.subplots(6, 1, figsize=(8, 12))
            if interactive:
                plt.ion()
                plt.show()

        for a in self.ax:
            a.clear()

        # Set the time steps
        t = np.linspace(0, x[-1], self.num_time_steps + 1)

        u = self.getControlTrajectory(x)
        self.ax[0].plot(t, np.radians(u), label='alpha')
        self.ax[1].plot(t, q[:, 0], label='v')
        self.ax[2].plot(t, np.radians(q[:, 1]), label='gamma')
        self.ax[3].plot(t, q[:, 2], label='h')
        self.ax[4].plot(t, q[:, 3], label='r')
        self.ax[5].plot(t, q[:, 4], label='m')

        for a in self.ax:
            a.legend()

        if interactive:
            plt.pause(0.01)

        self.fig.canvas.draw()

        return


def save_design_vector(x, filename="x_opt.npy"):
    """Save the optimized design variable vector to a file."""
    np.save(filename, x)
    print(f"Saved design vector to {os.path.abspath(filename)}")

def load_design_vector(filename="x_opt.npy"):
    """Load a saved design variable vector."""
    if os.path.exists(filename):
        x_loaded = np.load(filename)
        print(f"Loaded design vector from {os.path.abspath(filename)}")
        return x_loaded
    else:
        raise FileNotFoundError(f"No saved design vector found at {filename}")

def sweep_and_plot_cs_vs_adjoint(problem, x, p,
                                 hs=None,
                                 filename=None):
    """
    For each step size in hs, compute:
      adjoint_result = p^T · J   (once, via your evalConGradient)
      cs_result      = (evalCon(x + i*h*p).imag / h)
    then plot rel_error = ||adjoint_result - cs_result|| / ||cs_result||.
    """
    # default sweep
    if hs is None:
        hs = np.logspace(-12, -2, 5)

    # compute adjoint direction once
    dcdx           = problem.evalConGradient(x, interactive=False)
    adjoint_result = dcdx.dot(p)

    rel_errors = []
    for h in hs:
        cs_dir = problem.evalCon(x + 1j*h*p).imag / h
        err    = np.linalg.norm(adjoint_result - cs_dir)
        rel    = err / (np.linalg.norm(cs_dir) + 1e-16)
        rel_errors.append(rel)

    rel_errors = np.array(rel_errors)

    # now plot
    plt.figure(figsize=(5,4))
    plt.loglog(hs, rel_errors, 'o-', base=10)
    plt.xlabel("Complex-step size $h$")
    plt.ylabel("Relative error\n$\|Jp - \delta_{CS}\|/\\|\delta_{CS}\\|$")
    plt.title("Adjoint vs. CS Directional Derivative")
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.tight_layout()
    if filename:
        plt.savefig(filename, dpi=200)
    plt.show()

    return hs, rel_errors


def plot_dynamic_pressure_and_LD(problem, x_opt, Q):
    rho_vals = np.array([problem.system.getAtmDensity(h) for h in Q[:, 2]])
    v_vals = Q[:, 0]
    dyn_press = 0.5 * rho_vals * v_vals ** 2
    L_vals, D_vals = [], []
    for q, alpha in zip(Q, problem.getControlTrajectory(x_opt)):
        CL = problem.system.getCL(alpha, q[0], q[2])
        CD = problem.system.getCD(alpha, q[0], q[2])
        qinf = 0.5 * problem.system.getAtmDensity(q[2]) * q[0] ** 2
        L_vals.append(qinf * problem.system.S * CL)
        D_vals.append(qinf * problem.system.S * CD)
    LD = np.array(L_vals) / np.array(D_vals)
    t = np.linspace(0, x_opt[-1], Q.shape[0])
    plt.plot(t, dyn_press, label='Dynamic Pressure')
    plt.plot(t, LD, label='Lift-to-Drag')
    plt.xlabel("Time (s)")
    plt.legend()
    plt.title("Dynamic Pressure and L/D vs Time")
    plt.grid(True)
    plt.show()


def print_final_state_table(q_final, target):
    labels = ["v", "gamma", "h"]
    print("| State | Target | Actual | Residual |")
    print("|-------|--------|--------|----------|")
    for i, label in enumerate(labels):
        if target[i] is not None:
            residual = q_final[i] - target[i]
            print(f"| {label} | {target[i]:.2f} | {q_final[i]:.2f} | {residual:.2f} |")


def plot_trajectory_summary(t, q, u):
    labels = ['Velocity (ft/s)', 'Gamma (deg)', 'Altitude (ft)', 'Range (ft)', 'Mass (slugs)', 'Alpha (deg)']
    data = [q[:, 0], np.degrees(q[:, 1]), q[:, 2], q[:, 3], q[:, 4], np.degrees(u)]

    fig, axs = plt.subplots(3, 2, figsize=(12, 8))
    axs = axs.ravel()

    for i, ax in enumerate(axs):
        ax.plot(t, data[i])
        ax.set_title(labels[i])
        ax.grid(True)

    fig.tight_layout()
    plt.suptitle("Optimized Trajectory Results", fontsize=14, y=1.02)
    plt.show()


def print_final_state_table(q, target, labels=None):
    print("\nFinal State Comparison:")
    print("{:<15} {:>12} {:>12} {:>12}".format("State", "Final Value", "Target", "Error"))
    for i, val in enumerate(target):
        if val is not None:
            name = labels[i] if labels else f"State {i}"
            final_val = q[-1, i]
            error = final_val - val
            print(f"{name:<15} {final_val:12.4f} {val:12.4f} {error:12.4f}")

def export_trajectory_to_csv(t, q, u, filename="trajectory_results.csv"):
    df = pd.DataFrame(q, columns=['v', 'gamma', 'h', 'r', 'm'])
    df['t'] = t
    df['alpha_deg'] = np.degrees(u)
    df.to_csv(filename, index=False)
    print(f"Trajectory saved to {os.path.abspath(filename)}")

def plot_residual_convergence(residuals, k):
    """
    Plot the convergence of residuals during Newton iterations.

    Args:
        residuals: List or array of residual norms at each Newton iteration.
        k: Time step index for labeling purposes.
    """
    plt.figure(figsize=(6, 4))
    plt.semilogy(residuals, marker='o')
    plt.xlabel("Newton Iteration")
    plt.ylabel("Residual Norm")
    plt.title(f"Residual Convergence at Time Step {k}")
    plt.grid(True, which="both", ls="--")
    plt.tight_layout()
    plt.show()


# =====================================================
#     MAIN
# =====================================================
if __name__ == "__main__":
    # ====== Surrogate Model Check ======
    test_interp_gradient_raw(model_thrust)
    plot_interp_sensitivity(model_thrust)

    x = np.array([0.0, 1.0, 2.0, 3.0])
    y = np.array([0.0, 1.0, 0.0, 1.0])
    spline = CS1DCubicInterpolator(x, y)

    h = 1e-30
    x_test = 1.5
    val = spline(x_test + 1j * h)
    print("Value at 1.5:", val)
    print("CS Derivative:", val.imag / h)

    # ====== PROBLEM SETUP for CS vs. ADJOINT TEST ======
    tf_guess = 50.0  # seconds

    dynamics = SurrogateVehicleDynamics(surrogate_models)
    num_ctrl_pts = 10
    num_time_steps = int(tf_guess)
    problem = MinTimeTrajectory(dynamics, num_time_steps, num_ctrl_pts)

    num_design_vars = problem.num_design_vars
    x = np.zeros(num_design_vars)

    def make_shaped_initial_guess(num_ctrl_pts, tf_guess=tf_guess):
        t_grid = np.linspace(0, 1, num_ctrl_pts)
        alpha_guess = np.radians(2.75* t_grid)
        # alpha_guess = np.radians(2+2.75 * t_grid) #alternate AOA
        return np.concatenate([alpha_guess, [tf_guess]])

    x = make_shaped_initial_guess(num_ctrl_pts)

    plt.plot(np.linspace(0, x[-1], num_ctrl_pts), np.degrees(x[:-1]), label="Initial α guess")
    plt.xlabel("Time (s)")
    plt.ylabel("Angle of Attack (deg)")
    plt.title("Initial Control Guess (shaped)")
    plt.grid(True)
    plt.legend()
    plt.show()

    # ====== CS vs ADJOINT VALIDATION ======
    h = 1e-50
    p = np.random.randn(num_design_vars)
    p /= np.linalg.norm(p)  # Normalize the direction

    print("\n--- CS vs. Adjoint Gradient Check ---")

    print("\n--- Adjoint Gradient Check ---")
    dcdx = problem.evalConGradient(x, interactive=False)  # Adjoint-based full Jacobian
    adjoint_result = np.dot(dcdx, p)  # Adjoint directional derivative
    print("\n--- CS Check ---")
    cs_result = problem.evalCon(x + h * 1j * p).imag / h  # Complex-step directional derivative

    print("\n--- Begin Comparison ---")
    print(f"\n{'Constraint':>10s} {'Adjoint':>25s} {'CS Result':>25s} {'Rel. Error':>20s}")
    for i, (adj, cs) in enumerate(zip(adjoint_result, cs_result)):
        rel_err = np.abs((adj - cs) / (cs + 1e-16))  # avoid div-by-zero
        print(f"{i:10d} {adj:25.15e} {cs:25.15e} {rel_err:20.10e}")

    print("\n--- Plot Adjoint vs. CS Relative Error ---")
    rel_errors = []
    hp = []
    hs = np.logspace(-50, -30, 7)
    for h in hs:
        print(f"\n h = {h}")
        print(f"\n{'Constraint':>10s} {'Adjoint':>25s} {'CS Result':>25s} {'Rel. Error':>20s}")
        p = np.random.randn(num_design_vars)
        p /= np.linalg.norm(p)  # Normalize the direction
        dcdx = problem.evalConGradient(x, interactive=False)  # Adjoint-based full Jacobian
        adjoint_result = np.dot(dcdx, p)  # Adjoint directional derivative
        cs_result = problem.evalCon(x + h * 1j * p).imag / h  # Complex-step directional derivative
        for i, (adj, cs) in enumerate(zip(adjoint_result, cs_result)):
            rel_err = np.abs((adj - cs) / (cs + 1e-16))  # avoid div-by-zero
            print(f"{i:10d} {adj:25.15e} {cs:25.15e} {rel_err:20.10e}")
        rel_errors.append(rel_err)
        hp.append(h)


    print("h values:       ", hp)
    print("relative errors:", rel_errors)

    plt.figure(figsize=(5,4))
    plt.loglog(hp, rel_errors, 'o-', base=10)
    plt.xlabel("Complex‐step size $h$")
    plt.ylabel("Relative error $\|Jp - \delta_{CS}\| / \|\delta_{CS}\|$")
    plt.title("Adjoint vs. Complex‐Step Gradient Accuracy")
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.tight_layout()
    filename = "reduced_space_Adjoint_v_CS_plot.png"
    plt.savefig(filename)
    plt.show()


# ====== PROBLEM SETUP for OPTIMIZER ======
    tf_guess = 370.0  # seconds

    dynamics = SurrogateVehicleDynamics(surrogate_models)
    num_ctrl_pts = 37
    num_time_steps = int(tf_guess)
    problem = MinTimeTrajectory(dynamics, num_time_steps, num_ctrl_pts)

    num_design_vars = problem.num_design_vars
    x0 = np.zeros(num_design_vars)


    def make_shaped_initial_guess(num_ctrl_pts, tf_guess):
        t_grid = np.linspace(0, 1, num_ctrl_pts)
        alpha_guess = np.radians(2.0 + 0.5 * np.sin(2.5 * np.pi * t_grid))
        return np.concatenate([alpha_guess, [tf_guess]])


    try:
        x0 = load_design_vector("x_opt_reduced_space.npy")
    except FileNotFoundError:
        print("Could not load saved design vector. Using default initial guess.")
        x0 = make_shaped_initial_guess(num_ctrl_pts, tf_guess)

    plt.plot(np.linspace(0, x0[-1], num_ctrl_pts), np.degrees(x0[:-1]), label="Initial α guess")
    plt.xlabel("Time (s)")
    plt.ylabel("Angle of Attack (deg)")
    plt.title("Initial Control Guess (shaped)")
    plt.grid(True)
    plt.legend()
    plt.show()

    Q = problem.computeTrajectory(x0)
    print("Final state from initial guess:")
    print(f"  v = {Q[-1, 0]:.2f} ft/s")
    print(f"  gamma = {np.degrees(Q[-1, 1]):.2f} deg")
    print(f"  h = {Q[-1, 2]:.2f} ft")
    print(f"  r = {Q[-1, 3]:.2f} ft")
    print(f"  m = {Q[-1, 4]:.2f} slug")

    print("x[:5]:", x0[:5])
    print("q[:5, :]:", Q[:5, :])

    t = np.linspace(0, x0[-1], Q.shape[0])
    u = problem.getControlTrajectory(x0)

    fig, ax = plt.subplots(6, 1, figsize=(8, 12))

    ax[0].plot(t, np.degrees(u), label='alpha (deg)')
    ax[1].plot(t, Q[:, 0], label='v (ft/s)')
    ax[2].plot(t, np.degrees(Q[:, 1]), label='gamma (deg)')  # Convert from radians
    ax[3].plot(t, Q[:, 2], label='h (ft)')
    ax[4].plot(t, Q[:, 3], label='r (ft)')
    ax[5].plot(t, Q[:, 4], label='m (slug)')

    for a in ax:
        a.legend()
        a.set_xlim(t[0], t[-1])
        ydata = a.lines[0].get_ydata()
        if np.allclose(ydata, ydata[0]):
            a.set_ylim(ydata[0] - 1, ydata[0] + 1)
        a.grid(True)

    plt.tight_layout()
    plt.suptitle("Trajectory from Initial Guess", y=1.02)
    plt.show()



# ----- TEST OPTIMIZER - REDUCED SPACE ----- #

    # Set up the nonlinear constraint for the final state variables.
    lb = np.zeros(problem.num_constraints)
    ub = np.zeros(problem.num_constraints)
    con = NonlinearConstraint(problem.evalCon, lb, ub, jac=problem.evalConGradient)
    bounds = [(np.radians(-20.0), np.radians(20.0))] * (num_design_vars - 1) + [(100.0, None)]

    # Define output directory for plots
    plot_dir = r"C:\Users\ees00\Final Project AE 6310\plots_from_reduced_space"
    os.makedirs(plot_dir, exist_ok=True)

    objective_history = []
    def early_stop_callback(xk, state):
        obj = problem.evalObj(xk)
        objective_history.append(obj)

        print(f"Iteration {len(objective_history)}: objective = {obj:.6f}")

        constraint_violation = np.linalg.norm(problem.evalCon(xk))
        objective_history.append((len(objective_history), obj, constraint_violation))

        try:
            Q = problem.computeTrajectory(xk)
            t = np.linspace(0, xk[-1], Q.shape[0])
            u = problem.getControlTrajectory(xk)

            fig, ax = plt.subplots(6, 1, figsize=(8, 12))

            ax[0].plot(t, np.degrees(u), label='alpha (deg)')
            ax[1].plot(t, Q[:, 0], label='v (ft/s)')
            ax[2].plot(t, np.degrees(Q[:, 1]), label='gamma (deg)')
            ax[3].plot(t, Q[:, 2], label='h (ft)')
            ax[4].plot(t, Q[:, 3], label='r (ft)')
            ax[5].plot(t, Q[:, 4], label='m (slug)')

            for a in ax:
                a.legend()
                a.set_xlim(t[0], t[-1])
                ydata = a.lines[0].get_ydata()
                if np.allclose(ydata, ydata[0]):
                    a.set_ylim(ydata[0] - 1, ydata[0] + 1)
                a.grid(True)

            plt.tight_layout()
            plt.suptitle(f"Trajectory at Iteration {len(objective_history)}", y=1.02)

            # Save the figure
            filename = os.path.join(plot_dir, f"trajectory_iter_{len(objective_history):03d}.png")
            plt.savefig(filename)
            plt.pause(0.01)
            plt.close(fig)

        except Exception as e:
            print(f"Plotting failed: {e}")


    # Running the optimizer

    res = minimize(
        problem.evalObj, x0,
        method='trust-constr',
        jac=problem.evalObjGradient,
        hess=lambda x: np.zeros((num_design_vars, num_design_vars)),
        bounds=bounds,
        constraints=[con],
        callback=early_stop_callback,
        options={
            'initial_tr_radius': 30.0,
            'maxiter': 50,
            'gtol': 1,
            'xtol': 1e-3,
            'barrier_tol': 1e-4,
            'verbose': 3
        }
    )


    # Extract optimized values
    x_opt = res.x
    alphas = x_opt[:-1]
    tf_opt = x_opt[-1]
    Q = problem.computeTrajectory(x_opt)

    # Final state extraction
    V_idx, gamma_idx, h_idx, r_idx, m_idx = 0, 1, 2, 3, 4
    q_final = Q[-1]
    v_final = q_final[V_idx]
    gamma_final_rad = q_final[gamma_idx]
    h_final = q_final[h_idx]
    r_final = q_final[r_idx]
    m_final = q_final[m_idx]
    gamma_final_deg = np.degrees(gamma_final_rad)

    print("\n Final Optimized Trajectory State:")
    print(f"  Velocity (v)                = {v_final:.2f} ft/s")
    print(f"  Flight Path Angle (γ)       = {gamma_final_deg:.4f} deg")
    print(f"  Altitude (h)                = {h_final:.2f} ft")
    print(f"  Downrange (r)               = {r_final:.2f} ft")
    print(f"  Mass (m)                    = {m_final:.2f} slug")
    print(f"\n⏱ Final time: {tf_opt:.2f} sec")

    # Plot final trajectory
    t = np.linspace(0, x0[-1], Q.shape[0])
    u = problem.getControlTrajectory(x0)

    fig, ax = plt.subplots(6, 1, figsize=(8, 12))

    ax[0].plot(t, np.degrees(u), label='alpha (deg)')
    ax[1].plot(t, Q[:, 0], label='v (ft/s)')
    ax[2].plot(t, np.degrees(Q[:, 1]), label='gamma (deg)')  # Convert from radians
    ax[3].plot(t, Q[:, 2], label='h (ft)')
    ax[4].plot(t, Q[:, 3], label='r (ft)')
    ax[5].plot(t, Q[:, 4], label='m (slug)')

    for a in ax:
        a.legend()
        a.set_xlim(t[0], t[-1])
        ydata = a.lines[0].get_ydata()
        if np.allclose(ydata, ydata[0]):
            a.set_ylim(ydata[0] - 1, ydata[0] + 1)
        a.grid(True)

    plt.tight_layout()
    plt.suptitle("Final Trajectory")
    plt.show()

    print("Optimization result:")
    print(res)

    ## --------------Warm Re-start---------------- ##
    x0 = x_opt.copy()

    res = minimize(
        problem.evalObj, x0,
        method='trust-constr',
        jac=problem.evalObjGradient,
        hess=lambda x: np.zeros((num_design_vars, num_design_vars)),
        bounds=bounds,
        constraints=[con],
        callback=early_stop_callback,
        options={
            'initial_tr_radius': 30.0,
            'maxiter': 50,
            'gtol': 1,
            'xtol': 1e-3,
            'barrier_tol': 1e-4,
            'verbose': 3
        }
    )

    try:
        #  Save optimized design vector
        save_design_vector(res.x, filename="x_opt_reduced_space.npy")
    except Exception as e:
        # Handle any error that occurs
        print(f" Failed to save design vector: {e}")

    #  Extract optimized values
    x_opt = res.x
    alphas = x_opt[:-1]
    tf_opt = x_opt[-1]
    Q = problem.computeTrajectory(x_opt)

    #  Time and control trajectory
    t = np.linspace(0, tf_opt, problem.num_time_steps + 1)
    u = problem.getControlTrajectory(x_opt)

    try:
        #  Plot results
        plot_trajectory_summary(t, Q, u)
    except Exception as e:
        # Handle any error that occurs
        print(f"Failed to plot results: {e}")

    try:
        #  Print final state summary
        labels = ['Velocity', 'Gamma (rad)', 'Altitude', 'Range', 'Mass']
        print_final_state_table(Q, problem.qtarget, labels)
    except Exception as e:
        # Handle any error that occurs
        print(f"Failed to print final state summary: {e}")


    # Export to CSV
    export_trajectory_to_csv(t, Q, u)

    plot_dynamic_pressure_and_LD(problem, x_opt, Q)

    labels = ['Velocity (ft/s)', 'Gamma (rad)', 'Altitude (ft)', 'Range (ft)', 'Mass (slug)']
    print_final_state_table(Q, problem.qtarget, labels)

    # Save convergence history to CSV
    df_log = pd.DataFrame(objective_history, columns=["Iteration", "Objective", "ConstraintNorm"])
    df_log.to_csv("optimizer_convergence.csv", index=False)
    print("Saved optimizer convergence log to optimizer_convergence.csv")


