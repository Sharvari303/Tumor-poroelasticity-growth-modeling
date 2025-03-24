# ===== Loading packages required =======
using ModelingToolkit, MethodOfLines, LinearAlgebra, OrdinaryDiffEq, DomainSets
using ModelingToolkit: Differential
using Plots, Dates, ColorSchemes
using CSV, DataFrames, Statistics, LsqFit
import LsqFit: curve_fit
plotlyjs()

# # ===== Parameters initialization ======
@parameters t, x, y
@parameters L, k, xlog, D
@variables u(..)
Dx = Differential(x)
Dy = Differential(y)
Dt = Differential(t)
Dxx = Differential(x)^2
Dyy = Differential(y)^2

# ==== Load the temporal L,k,x0 csv input file =======
df = CSV.read("/Users/sharvari/Desktop/Kyle_Vining_ABM/jan27_logistic_fit_result_mar18.csv", DataFrame)

# ==== Define constant conditions for PDE solvers ======
t_min = 0.0
t_max = 0.8 #each PDE solve runs for 0.8 advection time 
x_min = -0.5
x_max = 0.5
y_min = -0.5
y_max = 0.5
dx = 0.025
dy = 0.025
#order = 3
x0 = 0
y0 = 0

# ==== SET PECLET NUMBER FOR SIMULATION ====
global Pe=2.208

# ==== Define velocity fields and extraction of velocity paramters from csv file  ====
function get_velocity_params(t)
    global df
    # Find the row where adv_time_scaled matches input t
    idx = findlast(x -> x ≤ t, df.adv_time_scaled) #findfirst(isequal(t), df.adv_time_scaled)
    
    # If no match is found, return nothing
    if isnothing(idx)
        return nothing
    end
    
    # Extract L, k, and xlog values for the matching time
    L = df.L_s[idx] #L here also incorporates the mobility term (permeability/eta)
    k = df.k_s[idx]
    xlog = df.xlog_s[idx]
    D= 1/Pe #1/150 #df.Pe[idx]

    return ModelingToolkit.MTKParameters([L, k, xlog, D], (), (), ())
end

function vx( x, y, p)
    
    L, k, xlog = p 
    r = sqrt((x - 0)^2 + (y - 0)^2)
    radial_vel = -L / (1 + exp(-k * (r - xlog)))
    return ifelse(r == 0, 0.0, radial_vel * (x - 0) / r) #(x-0)/r gives the direction of the velocity vector
end

function vy(x, y, p)
    
    L, k, xlog = p 
    r = sqrt((x - 0)^2 + (y - 0)^2)
    radial_vel =  -L / (1 + exp(-k * (r - xlog)))
    return ifelse(r==0, 0.0, radial_vel * (y - 0) / r) #(y-0)/r gives the direction of the velocity vector
end

eq = Dt(u(t, x, y)) ~ D*(Dxx(u(t, x, y)) + Dyy(u(t, x, y)))  - vx(x,y,[L, k, xlog])*Dx(u(t, x, y)) -  vy(x,y,[L, k,xlog])*Dy( u(t, x, y)) 
#eq = Dt(u(t, x, y)) ~ D*(Dxx(u(t, x, y)) + Dyy(u(t, x, y))) 

x0 = 0.0
y0 = 0.0

bcs = [
    u(t_min, x, y) ~ 1, #0,  # Initial Condition
    # # Robin BCs for left (x_min) and right (x_max) boundaries
    # α * u(t, x_min, y) - β * Dx(u(t, x_min, y)) ~ 0,
    # α * u(t, x_max, y) + β * Dx(u(t, x_max, y)) ~ 0,

    # # Robin BCs for bottom (y_min) and top (y_max) boundaries
    # α * u(t, x, y_min) - β * Dy(u(t, x, y_min)) ~ 0,
    # α * u(t, x, y_max) + β * Dy(u(t, x, y_max)) ~ 0
    u(t, x_min, y) ~ 0,  # Dirichlet BCs for left (x_min) boundary
    u(t, x_max, y) ~ 0,  # Dirichlet BCs for right (x_max) boundary
    u(t, x, y_min) ~ 0,  # Dirichlet BCs for bottom (y_min) boundary
    u(t, x, y_max) ~ 0   # Dirichlet BCs for top (y_max) boundary 
]

domains = [
            t ∈ Interval(t_min, t_max),
            x ∈ Interval(x_min, x_max),
            y ∈ Interval(y_min, y_max)
]

# ==== SET PARAMETERS AT t=0 FROM CSV FILE HERE ====
@named pdesys = PDESystem([eq], bcs, domains, [t, x, y], [u(t, x, y)], [D,L,k,xlog]; defaults = Dict([D => 1/Pe, L => -15.175, k => 9.338, xlog => 0.784])) 
@time discretization = MOLFiniteDifference([x => dx, y => dy],t; advection_scheme= UpwindScheme()) #,grid_adaptivity=true, flux_limiter="MC") 
@time prob = discretize(pdesys, discretization)
#println(prob.u0)

# === First Solve (Outside Loop) ===
global global_time = 0.00
global_time_end = 19.9 #16.2 #5.4 #advection time scaled max time + 1 adv time scale

println("Starting first PDE solve from t = 0.0 to t = 0.4...")
@time sol = solve(prob,TRBDF2(), dt = 0.01, saveat=0.1, verbose = true) #SSPRK54(), TRBDF2(), AutoTsit5(Rosenbrock23())

# === storage initialization of all solutions wrt time ====
solution_array = Array{Float64, 3}(undef, 0, 41, 41) #41x41 for 0.025, 21x21 for 0.05, 11x11 for 0.1
t_store = []
Nt_local = length(sol[t])

#extracting first solution for first global time step
solu = sol[u(t,x, y)] 
#concatenate cumulative solution with global time step 
solution_array = cat(solution_array, solu[1:end-1,:,:], dims=1) 
for t_local in 1:Nt_local-1  # Skip last step while keeping track of time
    append!(t_store, sol.t[t_local])
end
println(t_store)

#first final_state extraction after first global time step - to be used when remake is called for the first time in the loop
solucop = copy(sol[u(t,x, y)][end,:,:])
solucopwobc = solucop[2:end-1,2:end-1]
final_state = reshape(solucopwobc,:)

# Update the global time
global_time += t_max

# === Loop Over Temporal Steps ===
while round(global_time, digits = 1) <= global_time_end
    global global_time
    global solution_array
    global t_store
    global final_state
    global Nt_local
    global nsolu
    global solucop
    global solucopwobc
 
    # Get the updated parameters
    new_params = get_velocity_params(global_time)

    # Update the problem with the new parameters and IC
    new_prob = remake(prob, u0 = final_state, p=new_params)

    println("Starting PDE solve from t = $global_time to t = $(global_time + t_max)...")

    # Solve the PDE with the updated parameters
    new_sol = solve(new_prob, TRBDF2(), dt=0.01, saveat=0.1, verbose=true)

    # Extracting solution for each global time step
    # Extract the solution at the last time step - to feed as IC to the remake problem
    nsolu = new_sol[u(t,x, y)]
    solucop = copy(new_sol[u(t,x, y)][end,:,:])
    solucopwobc = solucop[2:end-1,2:end-1]
    final_state = reshape(solucopwobc,:)
    
    #concatenate cumulative nsolu with global time step and cumulative global time
    solution_array = cat(solution_array, nsolu[1:end-1,:,:], dims=1)
    for t_local in 1:Nt_local-1  # Skip last step
        append!(t_store, round(new_sol.t[t_local] + global_time, digits = 1))
    end
    println(t_store)
    # Update the global time
    global_time += t_max 
    #save cumulative time dependant solutions in a single array        
end

# # === Post-Processing ===

# Ensure Pe is filename-safe
Pe_str = replace(string(Pe), "." => "_")  # Convert "67.6" to "67_6"
##Convert "0.4" to "0_4"

#Make a gif of the solution
anim = @animate for i=1:length(t_store)
     maxval= maximum(solution_array)
     surface(solution_array[i,:,:], camera=(55.0, 30.0), size=(500,500), zlabel=("z"), 
             zlims=(0.0, maxval), xlabel = ("x"), ylabel=("y"), clims=(0.0,maxval), 
            title="t = $(round(t_store[i]*12.08, digits=3)) min, Pe=$Pe" )
end

# Get today's date
today_date = Dates.format(Dates.now(), "yyyy-mm-dd")  # Format as YYYY-MM-DD
#@time gif(anim, "Scenario1_$(today_date)_Pe$(Pe_str).gif", fps=5)  #####---- saving gif of resultant solution here ----####

#save surface plots at time steps at selected time points.
save_times =  [0, 0.1, 0.5, 1.6, 4.0] #, 19.2, 24.8] #[0, 0.6, 1.0, 1.6,4.0, 19.2, 24.8] 
indices = [argmin(abs.(t_store .- t)) for t in save_times]

# Save profile plots
for (idx, t) in zip(indices, save_times)

    maxval= maximum(solution_array)
    t_str= replace(string(t), "." => "_")
    plt = surface(solution_array[idx, :, :], camera=(55.0, 30.0), size=(500, 500),
                  zlabel="C", zlims=(0.0, maxval), xlabel="x", ylabel="y",
                  clims=(0.0, maxval), title="t = $(round(t*12.08, digits=3)) min, Pe=$Pe", titlefontsize=12, guidefontsize=10, tickfontsize=10)
    savefig(plt, "Scenario1_surfaceplots__t$(t_str)_Pe$(Pe_str).png")  #####---- saving surface plot of resultant solution here ----####

end

# make the radial concentration profiles
# Create radial grid after sol computation
r_grid = [sqrt((x - x0)^2 + (y - y0)^2) for x in range(x_min, x_max, step=dx), y in range(y_min, y_max, step=dy)]


function compute_radial_profile(sol, r_grid, discrete_t)
    #remove boudnaries from r_grid
    r_grid_inner = r_grid[1:end, 1:end]

    radial_profile_data = []
    for t_idx in 1:length(discrete_t)
        #Extract current snapshot of solution

        u_snapshot = sol[t_idx, 1:end, 1:end] 
        # println(typeof(u_snapshot)) 
        # println(u_snapshot)
        
        # Group by radial distance
        radial_concentration = Dict{Float64, Vector{Float64}}()
        for i in eachindex(r_grid_inner[:, 1])  
            r = round(r_grid_inner[i, 21], digits=3) # Round radial distance to 3 decimal places
            value = round(u_snapshot[i, 21], digits =3)  # Round concentration to 3 decimal places
            if haskey(radial_concentration, r)
                push!(radial_concentration[r], value)
            else
                radial_concentration[r] = [value]
            end
        end

        # Compute mean concentration per radius (handle empty bins)
        r_values = sort(collect(keys(radial_concentration)))  # Sorted radii
        #println(radial_concentration)
    
        # Sort keys (radial distances) and process avg values
        for r in r_values  # Sort r in ascending order
            avg_concentrations = mean(radial_concentration[r])
            push!(radial_profile_data, (discrete_t[t_idx], r, avg_concentrations))
        end
    end
    return radial_profile_data
end

radial_data = compute_radial_profile(solution_array, r_grid, t_store)

# Save to CSV
radial_df = DataFrame(time=[row[1] for row in radial_data],
                      radius=[row[2] for row in radial_data],
                      concentration=[row[3] for row in radial_data])

filename = "Scenario1_radialconcprofiles_$(today_date)_Pe$(Pe_str).csv"    ### -----saving csv file with radial conc of resultant solution here   ----###             
# Check if the file exists and delete it
if isfile(filename)
    rm(filename)
    #println("Existing file '$filename' deleted.")
end
                      
CSV.write(filename, radial_df)

##radial_df = CSV.read("Scenario2_radialconcprofiles_2025-03-21_Pe2_208.csv", DataFrame)
plot()
@time begin
    # Define specific time points to plot
    selected_times_2 = [0,0.1, 0.3, 0.4, 1.0, 1.6]  ##1.3, 2.5, 3.5, 4.5, 5.3, 10.3, 15.3, 16.3] #[0.2, 0.4, 0.6, 0.8, 1.0]
    
    # Define colormap
    num_times = length(selected_times_2)
    colors = get(ColorSchemes.viridis, range(0, 1, length=num_times))  # Get gradient colors

    # Filter the DataFrame for the selected time points
    filtered_df = filter(row -> row.time in selected_times_2, radial_df)

    # Group by time and plot each selected time points
    grouped = groupby(filtered_df, :time)
    for (i,g) in enumerate(grouped)
        plot!(g.radius[1:end]*8, g.concentration[1:end], label="t=$(round(g.time[1]*12.08, digits=1)) min", linewidth =2.5)
        #plot!(g.radius[1:end-1], g.concentration[1:end-1], label="", lw=2, color=colors[i]) 
    end

    xlabel!("Radius (mm)")
    ylabel!("Concentration (non-dimensionalized)")
    title!("Temporal Radial Concentration Profiles for Pe$(Pe)")
   # legend()  # Add the legend back to distinguish profiles

   # Add a colorbar legend to show time progression
   #cb = Colorbar(fontsize=10, color=colors, label="Time (advection time scale)")
    plot!(guidefontsize=13,    # Font size for axis labels (xlabel, ylabel)
    titlefontsize=13,    # Font size for title
    tickfontsize=13,     # Font size for axis ticks
    legendfontsize=13,   # Font size for legend
    legend_position=(1.01,0.75))  # Move legend to prevent overlap)
   savefig("Scenario2_radial_concentration_profiles_$(today_date)_Pe$(Pe_str).png")  ### -----saving pic with radial conc of resultant solution here   ----###  
end

