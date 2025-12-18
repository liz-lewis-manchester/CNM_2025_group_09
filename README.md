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
├── README.md                  # This file - project overview and guide
├── advection_model.ipynb      # Main notebook containing solver and all test cases
├── initial_conditions.csv     # Sample data for Test Case 2
└── results/                   # Output folder for figures (optional)
```

---

## Quick Start

### 1. Open the Notebook
Open `advection_model.ipynb` in Google Colab or Jupyter Notebook.

### 2. Run Cell 1 (Solver)
Run the first cell to load the solver functions. You should see:
```
✓ Solver loaded successfully!
```

### 3. Run Any Test Case
Run any of the test case cells (Cells 2-6) to see different scenarios.

---

## Notebook Structure

The notebook `advection_model.ipynb` contains:

| Cell | Contents |
|------|----------|
| **Cell 1** | Solver function + visualisation functions |
| **Cell 2** | Test Case 1: Basic simulation |
| **Cell 3** | Test Case 2: Initial conditions from CSV |
| **Cell 4** | Test Case 3: Sensitivity analysis |
| **Cell 5** | Test Case 4: Exponentially decaying boundary |
| **Cell 6** | Test Case 5: Variable velocity profile |

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
- `parameters` - simulation parameters

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

## Test Cases and How to Adapt Them

### Test Case 1: Basic Simulation

**What it does:** Runs a baseline simulation with constant parameters.

**Original code:**
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

**How to adapt:**

| Parameter | How to Change | Example Scenarios |
|-----------|---------------|-------------------|
| `end_time` | Increase for longer simulations | `3600` for 1 hour, `86400` for 1 day |
| `dt` | Decrease for more accuracy, increase for speed | `1` for detailed output, `60` for fast runs |
| `length` | Match your study area | `100` for 100m river section, `5000` for 5km |
| `dx` | Smaller = more detail, larger = faster | `0.1` for fine detail, `1.0` for coarse |
| `U` | Match your flow speed | `0.01` for slow groundwater, `2.0` for fast river |
| `C0` | Match your source concentration | Any positive value in your units |

**Example adaptations:**

*Heat transfer in a 50m pipe over 1 hour:*
```python
results = advection_solver(
    end_time=3600,   # 1 hour
    dt=30,           # 30s time step
    length=50,       # 50m pipe
    dx=0.5,          # 0.5m resolution
    U=0.5,           # 0.5 m/s flow
    C0=80            # 80°C inlet temperature
)
```

*Sediment transport over 1km in 24 hours:*
```python
results = advection_solver(
    end_time=86400,  # 24 hours
    dt=600,          # 10 minute time step
    length=1000,     # 1km reach
    dx=10,           # 10m resolution
    U=0.3,           # 0.3 m/s flow
    C0=500           # 500 mg/L sediment
)
```

---

### Test Case 2: Initial Conditions from CSV

**What it does:** Loads measured data from a file instead of assuming concentration only at x=0.

**Original code:**
```python
df = pd.read_csv('initial_conditions.csv', encoding='cp1252')
data_dist = df.iloc[:, 0].values
data_conc = df.iloc[:, 1].values

length, dx = 20, 0.2
num_points = int(length / dx) + 1
model_distance = np.linspace(0, length, num_points)
initial_cond = np.interp(model_distance, data_dist, data_conc)

results = advection_solver(
    ...,
    initial_conditions=initial_cond,
    boundary_condition=initial_cond[0]
)
```

**How to adapt:**

| What to Change | How | Example |
|----------------|-----|---------|
| CSV filename | Change `'initial_conditions.csv'` | `'my_measurements.csv'` |
| File encoding | Change `encoding=` parameter | `'utf-8'` for most files, `'cp1252'` for Windows |
| Column selection | Change `iloc[:, 0]` and `iloc[:, 1]` | `iloc[:, 2]` for third column |
| Column by name | Use column names instead of indices | `df['Distance']` and `df['Concentration']` |
| Domain size | Change `length` and `dx` to cover your data | Must span your measurement range |

**Example adaptations:**

*Using a CSV with named columns:*
```python
df = pd.read_csv('field_data.csv')
data_dist = df['distance_m'].values      # Column named 'distance_m'
data_conc = df['concentration'].values   # Column named 'concentration'
```

*Using data that covers 100m:*
```python
length, dx = 100, 1.0  # 100m domain, 1m resolution
num_points = int(length / dx) + 1
model_distance = np.linspace(0, length, num_points)
initial_cond = np.interp(model_distance, data_dist, data_conc)
```

*Setting a different boundary condition (e.g., zero at inlet):*
```python
results = advection_solver(
    ...,
    initial_conditions=initial_cond,
    boundary_condition=0  # No new pollutant entering
)
```

---

### Test Case 3: Sensitivity Analysis

**What it does:** Runs multiple simulations to see how results change with different parameter values.

**Original code:**
```python
for U in [0.001, 0.01, 0.1, 1, 10]:
    results = advection_solver(end_time=300, dt=10, length=20, dx=0.2, U=U, C0=250)
    ax.plot(results['distance'], results['concentrations'][-1], label=f'U = {U} m/s')
```

**How to adapt:**

| What to Change | How | Example |
|----------------|-----|---------|
| Parameter to vary | Change which parameter is in the loop | `for dt in [...]`, `for C0 in [...]` |
| Parameter values | Change the list of values | `[0.05, 0.1, 0.2]` for realistic range |
| Number of values | Add or remove values from list | More values = smoother sensitivity curve |
| Output to compare | Change `results['concentrations'][-1]` | `results['concentrations'][10]` for t=100s |

**Example adaptations:**

*Sensitivity to initial concentration:*
```python
for C0 in [100, 250, 500, 1000]:
    results = advection_solver(end_time=300, dt=10, length=20, dx=0.2, U=0.1, C0=C0)
    ax.plot(results['distance'], results['concentrations'][-1], label=f'C0 = {C0} µg/m³')
```

*Sensitivity to domain length:*
```python
for length in [10, 20, 50, 100]:
    results = advection_solver(end_time=300, dt=10, length=length, dx=0.2, U=0.1, C0=250)
    ax.plot(results['distance'], results['concentrations'][-1], label=f'L = {length} m')
```

*Compare results at different times:*
```python
for U in [0.05, 0.1, 0.2]:
    results = advection_solver(end_time=300, dt=10, length=20, dx=0.2, U=U, C0=250)
    mid_time_index = len(results['concentrations']) // 2  # Middle of simulation
    ax.plot(results['distance'], results['concentrations'][mid_time_index], label=f'U = {U} m/s')
ax.set_title("Concentration at t = 150s")
```

*Sensitivity to multiple parameters (nested loops):*
```python
for U in [0.05, 0.1, 0.2]:
    for dx in [0.1, 0.2, 0.5]:
        results = advection_solver(end_time=300, dt=10, length=20, dx=dx, U=U, C0=250)
        ax.plot(results['distance'], results['concentrations'][-1], 
                label=f'U={U}, dx={dx}')
```

---

### Test Case 4: Exponentially Decaying Boundary Condition

**What it does:** Models a source that changes over time using a function.

**Original code:**
```python
C0 = 250
decay_rate = 0.01

results = advection_solver(
    ...,
    boundary_condition=lambda t: C0 * np.exp(-decay_rate * t)
)
```

**How to adapt:**

| Boundary Pattern | Code | Use Case |
|------------------|------|----------|
| Constant | `boundary_condition=250` | Continuous source |
| Exponential decay | `lambda t: C0 * np.exp(-k*t)` | Finite spill depleting |
| Linear decay | `lambda t: C0 * max(0, 1 - t/T)` | Source draining at constant rate |
| Step off | `lambda t: C0 if t < T else 0` | Source removed at time T |
| Step on | `lambda t: 0 if t < T else C0` | Source starts at time T |
| Pulse | `lambda t: C0 if t % period < duration else 0` | Intermittent releases |
| Sinusoidal | `lambda t: C0 * (1 + A*np.sin(2*np.pi*t/T))` | Tidal or cyclic variation |
| Linear increase | `lambda t: C0 * (1 + t/T)` | Growing source |

**Example adaptations:**

*Source turns off after 2 minutes:*
```python
C0 = 250
shutoff_time = 120  # seconds

results = advection_solver(
    ...,
    boundary_condition=lambda t: C0 if t < shutoff_time else 0
)
plot_snapshots(results, title=f"Source stops at t={shutoff_time}s")
```

*Tidal variation (sinusoidal):*
```python
C0 = 250
amplitude = 0.5      # 50% variation
period = 600         # 10 minute cycle

results = advection_solver(
    ...,
    boundary_condition=lambda t: C0 * (1 + amplitude * np.sin(2*np.pi*t/period))
)
plot_snapshots(results, title="Tidal Variation")
```

*Pulsed release every 60 seconds for 10 seconds:*
```python
C0 = 250
pulse_interval = 60  # seconds between pulses
pulse_duration = 10  # seconds each pulse lasts

results = advection_solver(
    ...,
    boundary_condition=lambda t: C0 if (t % pulse_interval) < pulse_duration else 0
)
plot_snapshots(results, title="Pulsed Release")
```

*Linear decay to zero over simulation time:*
```python
C0 = 250
end_time = 300

results = advection_solver(
    end_time=end_time,
    ...,
    boundary_condition=lambda t: C0 * max(0, 1 - t/end_time)
)
plot_snapshots(results, title="Linear Decay to Zero")
```

*Faster or slower exponential decay:*
```python
# Slower decay (half-life ≈ 693 seconds)
decay_rate = 0.001
results = advection_solver(..., boundary_condition=lambda t: C0 * np.exp(-decay_rate * t))

# Faster decay (half-life ≈ 7 seconds)  
decay_rate = 0.1
results = advection_solver(..., boundary_condition=lambda t: C0 * np.exp(-decay_rate * t))
```

---

### Test Case 5: Variable Velocity Profile

**What it does:** Uses spatially varying velocity instead of a constant value.

**Original code:**
```python
length, dx = 20, 0.2
num_points = int(length / dx) + 1
base_U = 0.1
perturbation = 0.4

np.random.seed(42)
variable_U = base_U * (1 + perturbation * (2*np.random.random(num_points) - 1))

results = advection_solver(..., U=variable_U)
```

**How to adapt:**

| Velocity Pattern | Code | Use Case |
|------------------|------|----------|
| Constant | `U = 0.1` | Uniform flow |
| Random perturbation | `U = base * (1 + p*(2*np.random.random(n)-1))` | Turbulent flow |
| Linear increase | `U = base * (1 + distance/length)` | Accelerating flow |
| Linear decrease | `U = base * (1 - 0.5*distance/length)` | Decelerating flow |
| Step change | `U = np.where(distance < L/2, U1, U2)` | Channel width change |
| Sinusoidal | `U = base * (1 + A*np.sin(2*np.pi*distance/wavelength))` | Periodic variations |
| From measurements | `U = np.interp(model_dist, meas_dist, meas_velocity)` | Real velocity data |
| Parabolic profile | `U = Umax * (1 - (distance - L/2)**2 / (L/2)**2)` | Pipe flow profile |

**Example adaptations:**

*Change perturbation amount:*
```python
perturbation = 0.1   # ±10% variation (calm flow)
perturbation = 0.3   # ±30% variation (moderate turbulence)
perturbation = 0.5   # ±50% variation (highly turbulent)
```

*Different random pattern each run:*
```python
# Remove or change the seed
np.random.seed()  # Uses system time - different each run
# or
np.random.seed(123)  # Different but reproducible pattern
```

*Linear velocity increase (e.g., narrowing channel):*
```python
distance = np.linspace(0, length, num_points)
variable_U = base_U * (1 + distance/length)  # Doubles by end of domain
```

*Step change in velocity (e.g., enters wider section):*
```python
distance = np.linspace(0, length, num_points)
U_upstream = 0.2    # Faster in narrow section
U_downstream = 0.05 # Slower in wide section
change_point = 10   # metres

variable_U = np.where(distance < change_point, U_upstream, U_downstream)
```

*Sinusoidal velocity (e.g., meandering river):*
```python
distance = np.linspace(0, length, num_points)
wavelength = 5  # metres
amplitude = 0.3 # 30% variation

variable_U = base_U * (1 + amplitude * np.sin(2*np.pi*distance/wavelength))
```

*Velocity from measurement data:*
```python
# Your measured velocities at specific points
measured_distance = np.array([0, 5, 10, 15, 20])  # metres
measured_velocity = np.array([0.08, 0.12, 0.15, 0.10, 0.09])  # m/s

# Interpolate to model grid
distance = np.linspace(0, length, num_points)
variable_U = np.interp(distance, measured_distance, measured_velocity)
```

---

## General Adaptation Guide

### Changing Units

The solver is unit-agnostic. Just be consistent:

| System | length | dx | dt | U | C0 |
|--------|--------|----|----|---|-----|
| SI (default) | metres | metres | seconds | m/s | µg/m³ |
| Imperial | feet | feet | seconds | ft/s | lb/ft³ |
| Large scale | km | km | hours | km/hr | mg/L |

Remember to update axis labels in the plotting functions.

### Combining Test Case Techniques

You can combine multiple adaptations:

*Real data + decaying source + variable velocity:*
```python
# Load initial conditions from CSV (Test 2)
initial_cond = np.interp(model_distance, data_dist, data_conc)

# Create variable velocity (Test 5)
variable_U = base_U * (1 + 0.2 * (2*np.random.random(num_points) - 1))

# Run with decaying boundary (Test 4)
results = advection_solver(
    end_time=300, dt=10, length=20, dx=0.2,
    U=variable_U,
    initial_conditions=initial_cond,
    boundary_condition=lambda t: initial_cond[0] * np.exp(-0.01*t)
)
```

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| `NameError: advection_solver is not defined` | Run Cell 1 (solver) first |
| `UnicodeDecodeError` when reading CSV | Try `encoding='cp1252'` or `encoding='utf-8'` |
| Animation not displaying | Make sure you're in Google Colab or Jupyter Notebook |
| Results look wrong | Check units are consistent; try smaller dt and dx |
| Concentration goes negative | Check velocity direction and boundary conditions |

---

## Dependencies

- `numpy` - numerical calculations
- `matplotlib` - plotting
- `pandas` - reading CSV files (Test Case 2 only)
- `IPython` - animation display in notebooks

All are pre-installed in Google Colab.

---

## Authors

[Group 9 / Ching Yau Chan, Hassan Alhamdani, Jiongjie Chen, Lucas So and Oyinmiebi Youdeowei]

