# Tumor Poroelasticity Growth Modeling

This repository contains the scripts and data required to reproduce the computational results from our study on growth factor transport in mechanically compressed hydrogels, including PDE-based advection-diffusion simulations and agent-based tumor growth modeling.

---

## Table of Contents
- [Repository Contents](#repository-contents)
- [Reproducing the Results](#reproducing-the-results)
- [Data Availability](#data-availability)
- [License](#license)

---

## Repository Contents

### PDE Simulation (Julia)
| File | Description |
|------|-------------|
| `advection_diffusion_solver.jl` | Solves the 2D advection-diffusion equation using a moving-parameter velocity field derived from poroelastic FEM simulations |
| `logistic_fit_result.csv` | Time-varying velocity field parameters (input to the Julia solver) |

### Data Files
| File | Description |
|------|-------------|
| `supplementary_fig4-chemicalpotentialdata.csv` | Raw chemical potential (μ) distribution across radial distances and time points during 30% compression relaxation |
| `supplementary_fig4-logisticfit.csv` | Logistic fit parameters for spatial gradients of chemical potential; used as velocity field inputs in the PDE solver |
| `supplementary_fig5.csv` | Temporal evolution of solid volume fraction (ϕₛ) and hydraulic permeability (K) |

### Agent-Based Model (PhysiCell)
| File | Description |
|------|-------------|
| `PhysiCell_directory_poroelasticity.zip` | Reference working directory containing the PhysiCell framework and all configuration files used to reproduce the tumor growth simulations |

---
## Reproducing the Results

### Figure 8 — Advection-Diffusion PDE (Julia)

**Requirements**
- Julia 1.9+
- Packages: `MethodOfLines.jl`, `ModelingToolkit.jl`, `OrdinaryDiffEq.jl`
- `logistic_fit_result.csv` must be in the same directory as the script

**Key Parameters**

| Parameter | Value | Description |
|-----------|-------|-------------|
| `Pe` | `5.0` | Reproduces primary results in Figure 8 |
| `Pe` range | `0.22 – 220.8` | Explores diffusion-to-advection-dominated transition |

**Initial & Boundary Conditions**
- **IC:** `u(t_min, x, y) = 1` — domain starts fully saturated with growth factor
- **BC:** `u(t, x_min, y) = 0` — Dirichlet condition representing external sink/washout

**Outputs**

| File | Figure | Description |
|------|--------|-------------|
| `Scenario1_surfaceplots__t[TIME]_Pe5_0.png` | Figure 8B | 2D heatmaps of C/C₀ across the 8×8 mm domain |
| `Scenario1_radialconcprofiles_[DATE]_Pe5_0.csv` | — | Radial concentration profile data |
| `Scenario2_radial_concentration_profiles_[DATE]_Pe5_0.png` | Figure 8C | C/C₀ vs. radius for multiple time points |

---

### Figures 5, 6 & Supplementary Figure 9 — Agent-Based Model (PhysiCell)

**Requirements**
- C++ compiler with OpenMP
- PhysiCell v1.10.4 (included in `PhysiCell_directory_poroelasticity.zip`)

**Setup & Execution**

```bash
# 1. Extract the zip and navigate to the root directory
unzip PhysiCell_directory_poroelasticity.zip
cd PhysiCell_directory_poroelasticity

# 2. Confirm file placement
#    custom.cpp            → /custom_modules/
#    PhysiCell_settings.xml → /config/
#    cell_rules.csv         → /config/

# 3. Compile
make

# 4. Run (macOS/Linux)
./project ./config/PhysiCell_settings.xml
```

**Synchronizing with PDE Results**

To reproduce tumor growth curves for a specific Péclet number (Pe):
1. Open `PhysiCell_settings.xml` and update `<microenvironment_setup>` to match the GF initial concentration used in the corresponding Julia PDE run
2. Set the Dirichlet boundary conditions for the GF substrate to the concentration values sampled from the PDE output at the desired Pe and time point

**Key Variables & Parameters**

The following biological parameters govern the ABM simulation. These can be 
configured directly in PhysiCell by the user. For a complete working directory, 
refer to `PhysiCell_directory_poroelasticity.zip`.

| Parameter | Value | Description |
|-----------|-------|-------------|
| Simulation Domain | 400 µm × 400 µm | Represents Regions A and B within the larger PDE grid |
| Grid Resolution (dx, dy) | 20 µm | Spatial discretization |
| Initial Cell Population | 200 cells | Seeded at simulation start |
| Simulation Time | 10,080 min | 7-day period |
| GF Diffusion Coefficient | 100 µm²/s | Growth factor diffusion in the domain |
| Cellular Apoptosis Rate | 5.3×10⁻⁵ min⁻¹ | Baseline cell death rate |

**Cell Proliferation Function** (defined in `cell_rules.csv`)

Cell growth is governed by a sigmoidal (Hill) relationship to local growth 
factor (GF) concentration:

| Parameter | Value |
|-----------|-------|
| Max Proliferation Rate | 2×10⁻⁴ min⁻¹ |
| Half-max Concentration (K₀.₅) | 10⁻⁹ M |
| Hill Coefficient | 2.0 |

**Outputs**

| File | Description |
|------|-------------|
| `cell_populations.csv` | Total, alive, and dead cell counts over 10,080 minutes (7 days) |
| SVG/XML snapshots | Spatial distribution of cells relative to GF gradients, saved every 240 minutes |

---

## Data Availability

All original data files and the Julia script are deposited on Dryad: **[Dryad DOI here]**

> **Note:** The PhysiCell working directory (`PhysiCell_directory_poroelasticity.zip`) is provided in this repository for reproducibility but is **not** included in the Dryad deposit, as PhysiCell is third-party software.

---

## License

### Original Data & Scripts — CC0 1.0 (Public Domain)

[![License: CC0](https://img.shields.io/badge/License-CC0_1.0-lightgrey.svg)](https://creativecommons.org/publicdomain/zero/1.0/)

The following files are released under the [CC0 1.0 Universal (Public Domain Dedication)](https://creativecommons.org/publicdomain/zero/1.0/) license. You may use, copy, modify, and distribute them freely without restriction.

- `advection_diffusion_solver.jl`
- `logistic_fit_result.csv`
- `supplementary_fig4-chemicalpotentialdata.csv`
- `supplementary_fig4-logisticfit.csv`
- `supplementary_fig5.csv`

---

### PhysiCell & Associated Files — BSD 3-Clause License

[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

The PhysiCell framework (v1.10.4) and all files within `PhysiCell_directory_poroelasticity.zip` — including `custom.cpp`, `PhysiCell_settings.xml`, and `cell_rules.csv` — are distributed under the **BSD 3-Clause License**. See the [PhysiCell repository](https://github.com/MathCancer/PhysiCell) for full license terms.
