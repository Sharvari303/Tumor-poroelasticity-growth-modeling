Tumor Poroelasticity Growth Modeling
This repository contains the scripts and data required to reproduce the computational results from our study on growth factor transport in mechanically compressed hydrogels, including PDE-based advection-diffusion simulations and agent-based tumor growth modeling.

Repository Contents
PDE Simulation (Julia)

advection_diffusion_solver.jl — Solves the 2D advection-diffusion equation using a moving-parameter velocity field derived from poroelastic FEM simulations
logistic_fit_result.csv — Time-varying velocity field parameters (input to the Julia solver)

Data Files

supplementary_fig4-chemicalpotentialdata.csv — Raw chemical potential (μ) distribution across radial distances and time points during 30% compression relaxation
supplementary_fig4-logisticfit.csv — Logistic fit parameters for spatial gradients of chemical potential; used as velocity field inputs in the PDE solver
supplementary_fig5.csv — Temporal evolution of solid volume fraction (ϕₛ) and hydraulic permeability (K)

Agent-Based Model (PhysiCell)

PhysiCell_directory_poroelasticity.zip — Reference working directory containing the PhysiCell framework and all configuration files used to reproduce the tumor growth simulations


Reproducing the Results
Figure 8 — Advection-Diffusion PDE (Julia)
Requirements:

Julia 1.9+
Packages: MethodOfLines.jl, ModelingToolkit.jl, OrdinaryDiffEq.jl
logistic_fit_result.csv must be in the same directory as the script

Key parameters:

Set Pe = 5.0 to reproduce the primary results in Figure 8
Vary Pe between 0.22 and 220.8 to explore the diffusion-to-advection-dominated transition

Outputs:

Scenario1_surfaceplots__t[TIME]_Pe5_0.png — 2D heatmaps of C/C₀ across the 8×8 mm domain (Figure 8B)
Scenario1_radialconcprofiles_[DATE]_Pe5_0.csv — Radial concentration profiles
Scenario2_radial_concentration_profiles_[DATE]_Pe5_0.png — C/C₀ vs. radius plot for multiple time points (Figure 8C)


Figures 5, 6 & Supplementary Figure 9 — Agent-Based Model (PhysiCell)
Requirements:

C++ compiler with OpenMP
PhysiCell v1.10.4 (included in PhysiCell_directory_poroelasticity.zip)

Setup:

Place custom.cpp in /custom_modules/
Place PhysiCell_settings.xml and cell_rules.csv in /config/
Compile: make
Run: ./project ./config/PhysiCell_settings.xml (macOS/Linux)

Synchronizing with PDE results:

Map the GF initial concentration in PhysiCell_settings.xml (under <microenvironment_setup>) to the IC used in the Julia PDE run
Map the Dirichlet boundary conditions to the GF concentration sampled from the PDE output at the desired Pe and time point

Biological parameters (defined in cell_rules.csv):

Max proliferation rate: 2×10⁻⁴ min⁻¹
K₀.₅: 1×10⁻⁸ M
Hill coefficient: 2.0

Outputs:

cell_populations.csv — Total, alive, and dead cell counts over 10,080 minutes (7 days)
SVG and XML spatial snapshots every 240 minutes


Data Availability
All original data files and the Julia script are deposited on Dryad: [your Dryad DOI here].
The PhysiCell working directory (PhysiCell_directory_poroelasticity.zip) is provided in this repository only and is not included in the Dryad deposit.

License
Original Data & Scripts — CC0 1.0 (Public Domain)
All data files and the Julia analysis script created for this study are released under the CC0 1.0 Universal (Public Domain Dedication) license. You may use, copy, modify, and distribute them freely without restriction.
This includes:

advection_diffusion_solver.jl
logistic_fit_result.csv
supplementary_fig4-chemicalpotentialdata.csv
supplementary_fig4-logisticfit.csv
supplementary_fig5.csv

PhysiCell & Associated Files — BSD 3-Clause License
The PhysiCell framework (v1.10.4) and all files within PhysiCell_directory_poroelasticity.zip — including custom.cpp, PhysiCell_settings.xml, and cell_rules.csv — are distributed under the BSD 3-Clause License. See the PhysiCell repository for full license terms.
