# import numpy as np

# # Example data for radius and chemical potential
# r = np.array([0,0.000557306,0.00105409,0.00149683,0.00189142,0.002243,0.00255607,0.00283456,0.00308195,0.00330138,0.00349571,0.00366745,0.003819,0.00395269,0.00407052,0.00417406])  # Radii
# mu = np.array([531.496,531.512,531.493,531.089,529.396,524.8,514.937,497.108,468.909,429.314,378.683,316.56,244.749,166.967,85.6483,0])  # Chemical potential

# # Compute gradient using NumPy's gradient function and find velocity k/eta*gradient
# mu_gradient = np.gradient(mu, r)
# velocity = 10E-11*mu_gradient 

# print("Chemical potential gradient:", mu_gradient)

# import matplotlib.pyplot as plt

# #plt.plot(r, mu, label='Chemical Potential')
# # plt.plot(r, mu_gradient, label='Gradient of Chemical Potential', linestyle='--')
# # plt.xlabel('Radius (r)')
# # plt.ylabel('Value')
# # plt.legend()
# # plt.show()

# # Curve Fitting and Statistical Analysis
# from scipy.optimize import curve_fit  # For fitting models to data
# from sklearn.metrics import r2_score, mean_squared_error  # For calculating goodness-of-fit metrics

# #Using calculated chemical potential gradients for fitting a functin for radial velocity values.
# # Define the sigmoid (logistic) function
# def logistic_function(r, L, k, x0):
#     """
#     Logistic function.
#     L: maximum value (carrying capacity)
#     k: steepness of the curve
#     x0: midpoint (where the curve transitions)
#     """
#     return L / (1 + np.exp(-k * (r - x0)))

# # Provide initial guesses for parameters
# initial_guesses_logistic = [-8e5, 1000, 0.002]  # L, k, x0
# radius = r
# chemical_potential_gradient = mu_gradient
# # Fit the data using the logistic function
# popt_logistic, pcov_logistic = curve_fit(
#     logistic_function, radius, chemical_potential_gradient, p0=initial_guesses_logistic
# )

# # Extract the fitted parameters
# L, k, x0 = popt_logistic

# # Calculate fitted values using the logistic function
# fitted_values_logistic = logistic_function(radius, L, k, x0)

# # Calculate goodness of fit measures
# r2_logistic = r2_score(chemical_potential_gradient, fitted_values_logistic)
# mse_logistic = mean_squared_error(chemical_potential_gradient, fitted_values_logistic)

# # Plot the data and the logistic fit
# plt.figure(figsize=(8, 6))
# plt.plot(radius, chemical_potential_gradient, 'bo-', label="Real Data")
# plt.plot(radius, fitted_values_logistic, 'r--', label="Fitted Curve (Logistic)")
# plt.xlabel("Radius (r)")
# plt.ylabel("Chemical potential gradient (Pa/m)")
# plt.legend()
# plt.title(f"Logistic Fit: R² = {r2_logistic:.4f}, MSE = {mse_logistic:.4e}")
# plt.grid(True)
# plt.show()

# # Output goodness of fit measures and logistic parameters
# #r2_logistic, mse_logistic, popt_logistic
# print("Fitted Parameters (Logistic Function):")
# print(f"L (Maximum Gradient Magnitude): {L:.4e}")
# print(f"k (Steepness): {k:.4e}")
# print(f"x0 (Midpoint): {x0:.4e}")


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#import ace_tools as tools  # Ensure ace_tools is installed
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score, mean_squared_error

# Define logistic function for fitting
def logistic_function(r, L, k, x0):
    return L / (1 + np.exp(-k * (r - x0)))

# Load the CSV file
file_path = "/Users/sharvari/Desktop/Kyle_Vining_ABM/Convection-diffusion_setup_julia_sharvari - temporal_chem_pot_30_comp_5hrs.csv"
#file_path = "/Users/sharvari/Downloads/Convection-diffusion_setup_julia_sharvari - temporal_chem_pot_preliminary_version.csv"
data = pd.read_csv(file_path)
mobility_data = pd.read_csv("/Users/sharvari/Desktop/Kyle_Vining_ABM/mobility_data.csv")
# Extract unique time points
unique_times = data["time(min)"].unique()

# Initial guesses for logistic function fitting
initial_guesses_logistic = [-8e5, 1000, 0.002]

# Prepare storage for fitting parameters along with R² and MSE
fit_results = []

# Set up the figure for plotting
plt.figure(figsize=(10, 6))

# Iterate through each unique time, compute gradient, fit logistic function, and plot
for time in unique_times:
    subset = data[data["time(min)"] == time]  # Extract subset for the given time
    radii = subset["distance(m)"].values
    mu_values = subset["μ(Pa)"].values
    #mobility data
    

    #velocity calculation

    radii = np.array(radii, dtype=float)
    mu_values = np.array(mu_values, dtype=float)
    # Compute chemical potential gradient
    mu_gradient = np.gradient(mu_values, radii)
    velocity = mu_gradient*mobility_data[mobility_data["time(min)"]==time]["mobility"].values
    print(velocity)
    # Fit logistic function with updated initial guesses
    try:
        popt, _ = curve_fit(logistic_function, radii, velocity, p0=initial_guesses_logistic, maxfev=5000)
        L, k, x0 = popt

        # Generate fitted values
        fitted_values = logistic_function(radii, L, k, x0)

        # Compute R² and MSE
        r2 = r2_score(velocity, fitted_values)
        mse = mean_squared_error(velocity, fitted_values)

        # Store results
        fit_results.append([time, L, k, x0, r2, mse])

        # Plot the data and logistic fit
        plt.plot(radii, velocity, 'o', markersize=3, label=f"Data {time} min")
        plt.plot(radii, fitted_values, '--', label=f"Fit {time} min")

    except RuntimeError:
        fit_results.append([time, np.nan, np.nan, np.nan, np.nan, np.nan])

# Formatting plot
plt.xlabel("Radius (m)")
plt.ylabel("fluid velocity (m/s)")
plt.legend()
#plt.title("Logistic Fits for Chemical Potential Gradient at Different Times")
plt.grid(True)
plt.show()

# Convert results to a DataFrame and save as CSV
fit_df = pd.DataFrame(fit_results, columns=["time(min)", "L", "k", "xlog", "R²", "MSE"])
#fit_df = fit_df.merge(mobility_data, on="time(min)", how="left")
fit_df["vel_at_4mm"]=fit_df["L"]/(1+np.exp(-fit_df["k"]*(0.004-fit_df["xlog"]))) #chem pot grad is maximum at the periphery (so observed vel is max at 0.004 m)
max_vel_at_4mm = np.abs(fit_df["vel_at_4mm"].max())
fit_df["xplateau"]=fit_df["xlog"] + 2.2/fit_df["k"] #plateau is at 1/3 of the max value
fit_df["L_s"]=fit_df["L"]/max_vel_at_4mm
fit_df["k_s"]= fit_df["k"]* 0.008 # Characteristic length scale X = 0.008 m chosen to map original domain [0, 0.004] m to non-dimensional range [0, 0.5]
fit_df["xlog_s"]=fit_df["xlog"]/ 0.008 # Characteristic length scale X = 0.008 m chosen to map original domain [0, 0.004] m to non-dimensional range [0, 0.5]

# Round all values to 3 decimal places
fit_df = fit_df.round({"k": 3, "x0": 3, "R²": 3, "L_s": 3, "k_s": 3, "xlog_s": 3})

# Save the updated CSV file
csv_output_path = "/Users/sharvari/Desktop/Kyle_Vining_ABM/jan27_logistic_fit_result_mar18.csv"
fit_df.to_csv(csv_output_path, index=False)

# Define non-dimensionalized radius range [0, 0.5]
r_non_dim = np.linspace(0, 0.5, 100)

# Initialize a matrix for non-dimensional logistic fits
non_dim_fit_matrix = np.zeros((len(r_non_dim), len(fit_df)))

# Loop over each time point and compute the non-dimensional logistic function
for i in range(len(fit_df)):
    L_s = fit_df["L_s"].iloc[i]
    L = fit_df["L"].iloc[i] 
    k_s = fit_df["k_s"].iloc[i]
    xlog_s = fit_df["xlog_s"].iloc[i]

    # Compute non-dimensionalized logistic function
    non_dim_fit_matrix[:, i] = L_s* max_vel_at_4mm/ (1 + np.exp(-k_s * (r_non_dim - xlog_s)))

# Convert matrix to DataFrame for visualization
non_dim_fit_df = pd.DataFrame(non_dim_fit_matrix, index=r_non_dim, columns=fit_df["time(min)"])

# Plot non-dimensionalized logistic fits
plt.figure(figsize=(8, 5))
for i, time in enumerate(fit_df["time(min)"].values):
    arrayedv=np.array(non_dim_fit_df.iloc[:, i])
    plt.plot(r_non_dim, arrayedv, label=f"Time {int(time)} min")

plt.xlabel("Non-Dimensional Radius (0 to 0.5)")
plt.ylabel("Velocity (m/s)")
plt.title("  Radial velocity profiles on the non-dimensionalized lattice")
plt.legend()
plt.grid(True)
plt.show()
