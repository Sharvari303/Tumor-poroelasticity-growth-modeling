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
file_path = "Convection-diffusion_setup_julia_sharvari - temporal_chem_pot_30_comp_5hrs.csv"
#file_path = "Convection-diffusion_setup_julia_sharvari - temporal_chem_pot_preliminary_version.csv"
data = pd.read_csv(file_path)
mobility_data = pd.read_csv("mobility_data.csv")
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
