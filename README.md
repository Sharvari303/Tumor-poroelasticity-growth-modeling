# Tumor-poroelasticity-growth-modeling

1. Chempotgrad.py fits a logistic function to chemical potential gradient/velocity data generated from the FEM simulations.

2. advection-diffusion-PDE-solver.jl solves the PDE numerically and generates spatiotemporal growth factor concentration profiles under different Peclet number and Scenario runs.

3. Sigmoid plot - Cellular growth rate vs GF concentration.
  
4. Sample PhysiCell directory used for tumor growth ABM simulations.


Reproducing Main Figure 4 and Supplementary Figure 8: Advection-Diffusion PDE Solutions

This section describes the usage of the Julia script (advection_diffusion_solver.jl) to replicate the spatiotemporal concentration profiles and heatmaps.

1. Simulation Logic

The script solves the 2D advection-diffusion equation:
∂C/∂t​=D∇2C−v⋅∇C

where the velocity field v is derived from the poroelastic fluid efflux calculated in the Finite Element (FEM) simulations. The simulation uses a moving-parameter approach, where the velocity field parameters (L,k,xlog​) are updated iteratively from logistic_fit_result.csv to match the time-dependent relaxation of the hydrogel.

2. Configuration for Figure 8

To reproduce the specific panels in Figure 8, modify the following variables in the script:

    Péclet Number (Pe):

        Set global Pe = 5.0 to reproduce the primary results shown in Figure 8.

        Vary Pe between 0.22 and 220.8 to see the transition from diffusion-dominated to advection-dominated transport.

    Initial & Boundary Conditions:

        IC: u(t_min, x, y) ~ 1 (The domain starts fully saturated with growth factor).

        BC: u(t, x_min, y) ~ 0 (Dirichlet conditions at the boundaries represent the external sink/washout).

3. Generated Outputs

Running the script generates the following visualization files:

    For Figure 8B (Heatmaps):

        The script produces PNG files named Scenario1_surfaceplots__t[TIME]_Pe5_0.png.

        These surface plots represent the 2D spatial distribution of growth factor concentration (C/C0​) across the 8×8 mm lattice at discrete time points.

    For Figure 8C (Radial Profiles):

        The script calculates the mean concentration at each radial distance from the center.

        Output File: Scenario1_radialconcprofiles_[DATE]_Pe5_0.csv.

        Plot: Scenario2_radial_concentration_profiles_[DATE]_Pe5_0.png shows C/C0​ vs. Radius (mm) for multiple time points (e.g., t=0 to t=360 min).

4. Software Requirements

    Language: Julia (1.9+)

    Key Libraries: MethodOfLines.jl (for PDE discretization), ModelingToolkit.jl (for symbolic modeling), and OrdinaryDiffEq.jl.

    Input Dependency: Requires logistic_fit_result.csv in the same directory to provide the time-varying velocity parameters.
