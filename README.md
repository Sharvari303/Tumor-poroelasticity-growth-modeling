Reproducing Main Figure 4 and Supplementary Figure 8: Advection-Diffusion PDE Solutions

The script and input csvs are provided at this github link: https://github.com/Sharvari303/Tumor-poroelasticity-growth-modeling/tree/main

This section describes the usage of the Julia script (advection_diffusion_solver.jl) to replicate the spatiotemporal concentration profiles and heatmaps.

Simulation Logic

The script solves the 2D advection-diffusion equation:
∂C/∂t​=D∇2C−v⋅∇C

where the velocity field v is derived from the poroelastic fluid efflux calculated in the Finite Element (FEM) simulations. The simulation uses a moving-parameter approach, where the velocity field parameters (L,k,xlog​) are updated iteratively from logistic_fit_result.csv to match the time-dependent relaxation of the hydrogel.

To reproduce the specific panels in Figure 8, modify the following variables in the script:

Péclet Number (Pe):

    Set global Pe = 5.0 to reproduce the primary results shown in Figure 8.

    Vary Pe between 0.22 and 220.8 to see the transition from diffusion-dominated to advection-dominated transport.

Initial & Boundary Conditions:

    IC: u(t_min, x, y) ~ 1 (The domain starts fully saturated with growth factor).

    BC: u(t, x_min, y) ~ 0 (Dirichlet conditions at the boundaries represent the external sink/washout).

    2\. **Generated Outputs**

Running the script generates the following visualization files:

For Figure 8B (Heatmaps):

    The script produces PNG files named Scenario1_surfaceplots__t[TIME]_Pe5_0.png.

    These surface plots represent the 2D spatial distribution of growth factor concentration (C/C0​) across the 8×8 mm lattice at discrete time points.

For Figure 8C (Radial Profiles):

    The script calculates the mean concentration at each radial distance from the center.

    Output File: Scenario1_radialconcprofiles_[DATE]_Pe5_0.csv.

    Plot: Scenario2_radial_concentration_profiles_[DATE]_Pe5_0.png shows C/C0​ vs. Radius (mm) for multiple time points (e.g., t=0 to t=360 min).

Software Requirements:

Language: Julia (1.9+)

Key Libraries: MethodOfLines.jl (for PDE discretization), ModelingToolkit.jl (for symbolic modeling), and OrdinaryDiffEq.jl.

Input Dependency: Requires logistic_fit_result.csv in the same directory to provide the time-varying velocity parameters.

Supplementary Figure 4: Chemical Potential and Logistic Fit Parameters

These files contain the underlying mechanical driving forces calculated via Finite Element Method (FEM) and the subsequent mathematical fitting used for the advection-diffusion modeling.

File 1: supplementary_fig4-chemicalpotentialdata.csv

Description: Raw chemical potential (μ) distribution within the alginate gel across various radial distances and time points during the 30% compression relaxation.

Columns:

    time(min): The elapsed time of the simulation (Units: min).

    distance(m): Radial distance from the center of the 8mm diameter gel disc (Units: m).

    μ(Pa): Chemical potential value (Units: Pa).

File 2: supplementary_fig4-logisticfit.csv

Description: Parameters for the logistic functions fitted to the spatial gradients of the chemical potential. These parameters are used as inputs for the velocity field in the Julia transport script.

Columns:

    time(min): Time point for the fit (Units: min).

    L, k, xlog: Raw fitting parameters for the gradient profile.

    R², MSE: Goodness-of-fit metrics.

    L_s, k_s, xlog_s: Scaled parameters utilized directly in the Julia PDE solver to account for the mobility coefficient (M=k/η).

    adv_time_scaled: Non-dimensionalized time used for the advection-diffusion simulation steps.

Supplementary Figure 5: Network Densification and Permeability Decay

This file characterizes the physical state of the biopolymer network as it undergoes poroelastic relaxation.
File: supplementary_fig5.csv

Description: Temporal evolution of the solid phase density and the resulting hydraulic permeability. This data illustrates the "transport bottleneck" where increasing solid volume fraction (ϕs​) leads to a significant decrease in permeability (K).

Columns:

    Time(s): Simulation time (Units: seconds).

    Phi_s: Solid volume fraction of the alginate network (Dimensionless).

    K(m^2): Hydraulic permeability (Units: m$^2$).

Technical Context: The relationship between ϕs​ and K follows a power-law or exponential decay, which is a fundamental assumption in the poroelastic framework to explain restricted growth factor transport.

Figure 5, 6 and Supplementary Figure 9: Results from Agent-Based Modeling (PhysiCell)

This section describes the implementation of the agent-based models (ABM) used to simulate tumor growth response. The ABM utilizes the spatiotemporal results from the advection-diffusion PDE simulations as its biochemical input. (an example directory is PhysiCell_directory_poroelasticity.zip)

Simulation Overview

The ABM is implemented using PhysiCell (v1.10.4) to simulate how tumor cell proliferation is regulated by the transport limitations identified in the poroelastic continuum modeling.

How to Run the Model:

Environment: Requires a C++ compiler and OpenMP.

Ensure custom.cpp is in the /custom_modules/ directory. Place PhysiCell_settings.xml and cell_rules.csv in the /config/ directory.

Compilation: Execute make in the root directory.

Execution: Run the model with ./project ./config/PhysiCell_settings.xml (mac/linux)

 **2.Reproducing Specific Results (Tuning Parameters)**

To reproduce the tumor growth curves for a specific Péclet number (Pe), the PhysiCell environment must be synchronized with the concentration values obtained from the Julia PDE solver.

A. Mapping PDE Results to ABM Environment

Initial Conditions (IC): Map the GF initial concentration in the XML (under <microenvironment_setup>) to the initial condition used in the PDE run.

Boundary Conditions (BC): Map the Dirichlet BCs for the GF substrate in the XML to the growth factor concentration sampled from the PDE simulations at the specific Pe and time point of interest.

B. Biological Logic

Cell proliferation is assumed to be governed by a sigmoidal relationship between growth factor (GF) concentration and the cell growth rate, as defined in the cell_rules.csv parameters (Max rate: 2×10−4 min$^{-1}$, K0.5​: 1×10−8 M, Hill coefficient: 2.0).

Outputs:

cell_populations.csv: A temporal log generated by custom.cpp tracking Total_Cells, Alive, and Dead counts over the 10,080-minute (7-day) period.

Spatial Data: SVG and XML snapshots are generated every 240 minutes for visualizing the spatial distribution of cells relative to GF gradients.

