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

**Code:**
```python
results1 = advection_solver(
    end_time=300,    # Simulation duration (s)
    dt=10,           # Time step (s)
    length=20,       # Domain length (m)
    dx=0.2,          # Spatial resolution (m)
    U=0.1,           # Flow velocity (m/s)
    C0=250           # Inlet concentration (µg/m³)
)

plot_snapshots(results1, title="Test Case 1: Basic Simulation")
create_animation(results1, max_conc=250)
```

**How to adapt:**

| Parameter | How to Change | Example Scenarios |
|-----------|---------------|-------------------|
| `end_time` | Increase for longer simulations | `3600` for 1 hour, `86400` for 1 day |
| `dt` | Decrease for more accuracy, increase for speed | `1` for detailed output, `60` for fast runs |
| `length` | Match your study area | `100` for 100m river section, `5000` for 5km |
| `dx` | Smaller = more detail, larger = faster | `0.1` for fine detail, `1.0` for coarse |
| `U` | Match your flow speed | `0.01` for slow groundwater, `2.0` for fast river |
| `C0` | Match your source concentration | Any positive value in your units |

---

### Test Case 2: Initial Conditions from CSV

**Description:** Test the case where, for the same model domain, the initial conditions are read in from a csv file 'initial_conditions.csv'. Note that the provided measurements are not aligned with the model grid. Your code should be written so that it can read in any initial conditions provided and to interpolate them onto the model grid.

**Code:**
```python
import pandas as pd

# Read CSV file
df = pd.read_csv('initial_conditions.csv', encoding='cp1252')
data_dist = df.iloc[:, 0].values  # Distance column
data_conc = df.iloc[:, 1].values  # Concentration column

# Create model grid and interpolate
length, dx = 20, 0.2
num_points = int(length / dx) + 1
model_distance = np.linspace(0, length, num_points)
initial_cond = np.interp(model_distance, data_dist, data_conc)

# Run solver
results2 = advection_solver(
    end_time=300, dt=10, length=20, dx=0.2, U=0.1,
    initial_conditions=initial_cond,
    boundary_condition=initial_cond[0]
)

plot_snapshots(results2, title="Test Case 2: CSV Initial Conditions")
create_animation(results2)
```

**How to adapt:**

| What to Change | How | Example |
|----------------|-----|---------|
| CSV filename | Change `'initial_conditions.csv'` | `'my_measurements.csv'` |
| File encoding | Change `encoding=` parameter | `'utf-8'` for most files |
| Column selection | Change `iloc[:, 0]` and `iloc[:, 1]` | `iloc[:, 2]` for third column |
| Column by name | Use column names instead of indices | `df['Distance']` and `df['Concentration']` |

---

### Test Case 3: Sensitivity Analysis

**Description:** Test to see how sensitive your model results are to its parameters (U, spatial and temporal resolution)

**Code:**
```python
# Sensitivity to Velocity (U)
fig, ax = plt.subplots(figsize=(10, 6))
for U in [0.001, 0.01, 0.1, 1, 10]:
    results = advection_solver(end_time=300, dt=10, length=20, dx=0.2, U=U, C0=250)
    ax.plot(results['distance'], results['concentrations'][-1],
            label=f'U = {U} m/s', linewidth=2)
ax.set_xlabel("Distance (m)")
ax.set_ylabel("Concentration (μg/m³)")
ax.set_title("Sensitivity to Velocity U (Final State)")
ax.legend()
ax.grid(True)
plt.show()

# Sensitivity to Spatial Resolution (dx)
fig, ax = plt.subplots(figsize=(10, 6))
for dx in [0.1, 0.2, 0.5, 1.0, 2.0]:
    results = advection_solver(end_time=300, dt=10, length=20, dx=dx, U=0.1, C0=250)
    ax.plot(results['distance'], results['concentrations'][-1],
            label=f'Δx = {dx} m', linewidth=2)

# Sensitivity to Time Step (dt)
fig, ax = plt.subplots(figsize=(10, 6))
for dt in [1, 5, 10, 30, 60]:
    results = advection_solver(end_time=300, dt=dt, length=20, dx=0.2, U=0.1, C0=250)
    ax.plot(results['distance'], results['concentrations'][-1],
            label=f'Δt = {dt} s', linewidth=2)
```

**How to adapt:**

| What to Change | How | Example |
|----------------|-----|---------|
| Parameter to vary | Change which parameter is in the loop | `for C0 in [100, 250, 500]` |
| Parameter values | Change the list of values | `[0.05, 0.1, 0.2]` for realistic range |
| Output to compare | Change `results['concentrations'][-1]` | `results['concentrations'][10]` for t=100s |

---

### Test Case 4: Exponentially Decaying Boundary Condition

**Description:** Test to see how an exponentially decaying initial concentration of the pollutant in time alters your results.

**Code:**
```python
C0 = 250           # Initial concentration (µg/m³)
decay_rate = 0.01  # Decay rate λ (1/s)

results4 = advection_solver(
    end_time=300, dt=10, length=20, dx=0.2, U=0.1, C0=C0,
    boundary_condition=lambda t: C0 * np.exp(-decay_rate * t)
)

plot_snapshots(results4, title=f"Test Case 4: Exponential Decay (λ={decay_rate})")
create_animation(results4, max_conc=C0)
```

**How to adapt:**

| Boundary Pattern | Code | Use Case |
|------------------|------|----------|
| Constant | `boundary_condition=250` | Continuous source |
| Exponential decay | `lambda t: C0 * np.exp(-k*t)` | Finite spill depleting |
| Linear decay | `lambda t: C0 * max(0, 1 - t/T)` | Source draining at constant rate |
| Step off | `lambda t: C0 if t < T else 0` | Source removed at time T |
| Sinusoidal | `lambda t: C0 * (1 + A*np.sin(2*np.pi*t/T))` | Tidal or cyclic variation |

---

### Test Case 5: Variable Velocity Profile

**Description:** Test to see how a variable stream velocity profile alters your results (for example, add a 10% random perturbation to the constant velocity profile).

**Code:**
```python
length, dx = 20, 0.2
num_points = int(length / dx) + 1
base_U = 0.1           # Base velocity (m/s)
perturbation = 0.4     # Perturbation fraction (±40%)

# Generate random velocity profile
np.random.seed(42)  # For reproducibility
variable_U = base_U * (1 + perturbation * (2*np.random.random(num_points) - 1))

# Plot velocity profile
fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(np.linspace(0, length, num_points), variable_U, 'b-', linewidth=2)
ax.axhline(base_U, color='r', linestyle='--', label='Base U')
ax.set_xlabel("Distance (m)")
ax.set_ylabel("Velocity (m/s)")
ax.set_title(f"Variable Velocity Profile (±{perturbation*100:.0f}%)")
ax.legend()
ax.grid(True)
plt.show()

# Run simulation with variable velocity
results5 = advection_solver(
    end_time=300, dt=10, length=20, dx=0.2,
    U=variable_U,
    C0=250
)

plot_snapshots(results5, title=f"Test Case 5: Variable Velocity (±{perturbation*100:.0f}%)")
create_animation(results5, max_conc=250)
```

**How to adapt:**

| Velocity Pattern | Code | Use Case |
|------------------|------|----------|
| Random perturbation | `base * (1 + p*(2*np.random.random(n)-1))` | Turbulent flow |
| Linear increase | `base * (1 + distance/length)` | Accelerating flow |
| Step change | `np.where(distance < L/2, U1, U2)` | Channel width change |
| From measurements | `np.interp(model_dist, meas_dist, meas_velocity)` | Real velocity data |

---

## Running Tests

To verify the solver works correctly, run the automated tests:

```python
# In Colab or terminal
!pip install pytest
!pytest tests/ -v
```

Expected output:
```
tests/test_unit.py::TestSpatialGrid::test_grid_size PASSED
tests/test_unit.py::TestTimeArray::test_time_frames PASSED
...
tests/test_integration.py::TestCase1BasicSimulation::test_simulation_runs PASSED
...
==================== 14 passed ====================
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
