# Pollutant Transport Model - 1D Advection Solver

## Overview

This project simulates pollutant transport in a river using the **1D advection equation**:

$$\frac{\partial\theta}{\partial t} = -U\frac{\partial\theta}{\partial x}$$

Where:
- θ = pollutant concentration (µg/m³)
- t = time (s)
- x = distance along the river (m)
- U = flow velocity (m/s)

The solver uses an **implicit finite difference scheme**, which is stable for a wide range of parameters.

---

## File Structure

```
├── README.md                       # This file - project overview and guide
├── data/
│   └── initial_conditions.csv      # Sample data for Test Case 2
├── notebooks/
│   └── CNM_CW_Group_9.ipynb        # Main notebook containing solver and all test cases
├── src/
│   └── advection_solver.py         # Standalone solver module
├── tests/
│   ├── test_unit.py                # Unit tests for individual functions
│   └── test_integration.py         # Integration tests for full workflows
└── results/                        # Output folder for figures and animations
```

---

## Quick Start

### 1. Open the Notebook

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/liz-lewis-manchester/CNM_2025_group_09/blob/main/notebooks/CNM_CW_Group_9.ipynb)

Or open `notebooks/CNM_CW_Group_9.ipynb` in Google Colab or Jupyter Notebook.

### 2. Run the Solver Cell
Run the first code cell (Main Code/Solver) to load the solver functions. You should see:
```
✓ Solver loaded successfully!
```

### 3. Run Any Test Case
Run any of the test case cells to see different scenarios.

---

## Notebook Structure

The notebook `CNM_CW_Group_9.ipynb` contains:

| Section | Contents |
|---------|----------|
| **Main Code (Solver)** | `advection_solver()`, `create_animation()`, `plot_snapshots()` functions |
| **Test Case 1** | Basic simulation with constant parameters |
| **Test Case 2** | Initial conditions from CSV file |
| **Test Case 3** | Sensitivity analysis (U, dx, dt) |
| **Test Case 4** | Exponentially decaying boundary condition |
| **Test Case 5** | Variable velocity profile |

---

## Solver Function

### `advection_solver()`

The main function that runs the simulation.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `end_time` | float | Total simulation duration (seconds) |
| `dt` | float | Time step size (seconds) |
| `length` | float | Domain length (metres) |
| `dx` | float | Spatial resolution (metres) |
| `U` | float or array | Flow velocity - constant or spatially varying |
| `C0` | float | Default concentration at inlet (µg/m³) |
| `initial_conditions` | array (optional) | Custom initial concentration profile |
| `boundary_condition` | float or function (optional) | Inlet concentration - constant or time-varying |

**Returns:** Dictionary containing:
- `distance` - spatial grid array
- `time` - time array  
- `concentrations` - list of concentration profiles at each time step
- `parameters` - simulation parameters (end_time, dt, dx, length, num_points)

**Example:**
```python
results = advection_solver(
    end_time=300,    # 5 minutes
    dt=10,           # 10s time step
    length=20,       # 20m domain
    dx=0.2,          # 0.2m resolution
    U=0.1,           # 0.1 m/s velocity
    C0=250           # 250 µg/m³
)
```

### `plot_snapshots()`

Creates a static plot showing concentration profiles at multiple times.

```python
plot_snapshots(results, num_snapshots=5, title="My Title")
```

### `create_animation()`

Creates an animated plot showing concentration evolution.

```python
create_animation(results, max_conc=250)
```

---

## Test Cases

### Test Case 1: Basic Simulation

**Description:** Test the case where the 1D model domain extends to 20m downstream (with a 20cm spatial resolution) of the point that the pollutant enters the river and model how the pollutant moves over the 5 minutes after it enters the river (with a temporal resolution of 10s). Assume that the initial concentration of the pollutant is 250 µg/m³ at x=0 and 0 elsewhere. Assume that U = 0.1ms⁻¹

**Where to make changes:** Modify the parameters passed to `advection_solver()`:
```python
results1 = advection_solver(
    end_time=300,    # ← Change simulation duration here
    dt=10,           # ← Change time step here
    length=20,       # ← Change domain length here
    dx=0.2,          # ← Change spatial resolution here
    U=0.1,           # ← Change velocity here
    C0=250           # ← Change inlet concentration here
)
```

---

### Test Case 2: Initial Conditions from CSV

**Description:** Test the case where, for the same model domain, the initial conditions are read in from a csv file 'initial_conditions.csv'. Note that the provided measurements are not aligned with the model grid. Your code should be written so that it can read in any initial conditions provided and to interpolate them onto the model grid.

**Where to make changes:**
```python
# ← Change filename and encoding here
df = pd.read_csv('initial_conditions.csv', encoding='cp1252')

# ← Change column selection here (0 = first column, 1 = second column)
data_dist = df.iloc[:, 0].values
data_conc = df.iloc[:, 1].values

# ← Change domain parameters here
length, dx = 20, 0.2
```

---

### Test Case 3: Sensitivity Analysis

**Description:** Test to see how sensitive your model results are to its parameters (U, spatial and temporal resolution)

**Where to make changes:**
```python
# ← Change the parameter values in these lists
for U in [0.001, 0.01, 0.1, 1, 10]:
    ...

for dx in [0.1, 0.2, 0.5, 1.0, 2.0]:
    ...

for dt in [1, 5, 10, 30, 60]:
    ...
```

---

### Test Case 4: Exponentially Decaying Boundary Condition

**Description:** Test to see how an exponentially decaying initial concentration of the pollutant in time alters your results.

**Where to make changes:**
```python
C0 = 250           # ← Change initial concentration here
decay_rate = 0.01  # ← Change decay rate here

results4 = advection_solver(
    ...
    # ← Change the boundary condition function here
    boundary_condition=lambda t: C0 * np.exp(-decay_rate * t)
)
```

---

### Test Case 5: Variable Velocity Profile

**Description:** Test to see how a variable stream velocity profile alters your results (for example, add a 10% random perturbation to the constant velocity profile).

**Where to make changes:**
```python
base_U = 0.1           # ← Change base velocity here
perturbation = 0.4     # ← Change perturbation fraction here (0.4 = ±40%)

np.random.seed(42)     # ← Change seed for different random pattern, or remove for new pattern each run

# ← Modify this line to change the velocity profile formula
variable_U = base_U * (1 + perturbation * (2*np.random.random(num_points) - 1))
```

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| `NameError: advection_solver is not defined` | Run the Main Code (Solver) cell first |
| `FileNotFoundError` for CSV | Ensure `initial_conditions.csv` is in the correct path |
| `UnicodeDecodeError` when reading CSV | Try `encoding='utf-8'` instead of `encoding='cp1252'` |
| Animation not displaying | Make sure you're in Google Colab or Jupyter Notebook |

---

## Dependencies

- `numpy` - numerical calculations
- `matplotlib` - plotting
- `pandas` - reading CSV files (Test Case 2 only)
- `IPython` - animation display in notebooks

All are pre-installed in Google Colab.

---

## Authors

Group 9: Ching Yau Chan, Hassan Alhamdani, Jiongjie Chen, Lucas So and Oyinmiebi Youdeowei
