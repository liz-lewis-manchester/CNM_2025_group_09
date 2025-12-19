"""
Integration Tests for Advection Solver
======================================
Tests full simulation workflows matching the 5 test cases in the notebook.

Each test class corresponds to a test case in CNM_CW_Group_9.ipynb

Run with: python -m pytest tests/test_integration.py -v
Or in Colab: !pytest tests/test_integration.py -v
"""

import numpy as np
import sys
sys.path.append('../src')
from advection_solver import advection_solver



# TEST CASE 1: BASIC SIMULATION

class TestCase1BasicSimulation:
    """
    Integration tests for Test Case 1: Basic Simulation.
    
    From notebook:
        "Test the case where the 1D model domain extends to 20m downstream 
        (with a 20cm spatial resolution) of the point that the pollutant 
        enters the river and model how the pollutant moves over the 5 minutes 
        after it enters the river (with a temporal resolution of 10s). 
        Assume that the initial concentration of the pollutant is 250 µg/m³ 
        at x=0 and 0 elsewhere. Assume that U = 0.1ms-1"
    """
    
    def test_basic_simulation_runs(self):
        """Check basic simulation completes without error."""
        results = advection_solver(
            end_time=300,    # 5 minutes
            dt=10,           # 10s temporal resolution
            length=20,       # 20m downstream
            dx=0.2,          # 20cm spatial resolution
            U=0.1,           # U = 0.1 m/s
            C0=250           # 250 µg/m³ at x=0
        )
        
        assert results is not None
        assert 'concentrations' in results
    
    def test_correct_number_of_frames(self):
        """Check simulation produces 31 frames (300/10 + 1)."""
        results = advection_solver(
            end_time=300, dt=10, length=20, dx=0.2, U=0.1, C0=250
        )
        
        assert len(results['concentrations']) == 31
    
    def test_correct_number_of_points(self):
        """Check simulation produces 101 spatial points (20/0.2 + 1)."""
        results = advection_solver(
            end_time=300, dt=10, length=20, dx=0.2, U=0.1, C0=250
        )
        
        assert len(results['distance']) == 101
        assert len(results['concentrations'][0]) == 101



# TEST CASE 2: CSV INITIAL CONDITIONS

class TestCase2CSVInitialConditions:
    """
    Integration tests for Test Case 2: CSV Initial Conditions.
    
    From notebook:
        "Test the case where, for the same model domain, the initial conditions 
        are read in from a csv file 'initial_conditions.csv'. Note that the 
        provided measurements are not aligned with the model grid. Your code 
        should be written so that it can read in any initial conditions provided 
        and to interpolate them onto the model grid."
    """
    
    def test_interpolation_workflow(self):
        """
        Check full workflow: create data, interpolate with np.interp, simulate.
        
        From notebook:
            initial_cond = np.interp(model_distance, data_dist, data_conc)
        """
        # Simulate CSV data
        data_dist = np.array([0, 5, 10, 15, 20])
        data_conc = np.array([300, 200, 100, 50, 10])
        
        # Interpolate to model grid (same as notebook)
        length, dx = 20, 0.2
        num_points = int(length / dx) + 1
        model_distance = np.linspace(0, length, num_points)
        initial_cond = np.interp(model_distance, data_dist, data_conc)
        
        # Run simulation (same as notebook)
        results = advection_solver(
            end_time=300, dt=10, length=20, dx=0.2, U=0.1,
            initial_conditions=initial_cond,
            boundary_condition=initial_cond[0]
        )
        
        assert results is not None
        assert results['concentrations'][0][0] == 300
    
    def test_interpolation_produces_correct_values(self):
        """Check np.interp produces correct intermediate values."""
        data_dist = np.array([0, 10, 20])
        data_conc = np.array([100, 200, 100])
        
        model_distance = np.linspace(0, 20, 101)
        initial_cond = np.interp(model_distance, data_dist, data_conc)
        
        assert initial_cond[0] == 100    # At x=0
        assert initial_cond[50] == 200   # At x=10 (middle)
        assert initial_cond[100] == 100  # At x=20



# TEST CASE 3: SENSITIVITY ANALYSIS

class TestCase3SensitivityAnalysis:
    """
    Integration tests for Test Case 3: Sensitivity Analysis.
    
    From notebook:
        "Test to see how sensitive your model results are to its parameters 
        (U, spatial and temporal resolution)"
    """
    
    def test_sensitivity_to_velocity(self):
        """
        Check different velocities produce different results.
        
        From notebook: for U in [0.001, 0.01, 0.1, 1, 10]
        """
        results_slow = advection_solver(
            end_time=300, dt=10, length=20, dx=0.2, U=0.01, C0=250
        )
        results_fast = advection_solver(
            end_time=300, dt=10, length=20, dx=0.2, U=1, C0=250
        )
        
        # Faster velocity should spread pollutant further
        slow_spread = np.sum(results_slow['concentrations'][-1] > 1)
        fast_spread = np.sum(results_fast['concentrations'][-1] > 1)
        
        assert fast_spread > slow_spread
    
    def test_sensitivity_to_spatial_resolution(self):
        """
        Check different dx values all produce results.
        
        From notebook: for dx in [0.1, 0.2, 0.5, 1.0, 2.0]
        """
        for dx in [0.1, 0.2, 0.5, 1.0, 2.0]:
            results = advection_solver(
                end_time=300, dt=10, length=20, dx=dx, U=0.1, C0=250
            )
            assert results is not None
            assert len(results['concentrations']) > 0
    
    def test_sensitivity_to_time_step(self):
        """
        Check different dt values produce correct number of frames.
        
        From notebook: for dt in [1, 5, 10, 30, 60]
        """
        for dt in [1, 5, 10, 30, 60]:
            results = advection_solver(
                end_time=300, dt=dt, length=20, dx=0.2, U=0.1, C0=250
            )
            
            expected_frames = int(300 / dt) + 1
            assert len(results['concentrations']) == expected_frames



# TEST CASE 4: EXPONENTIAL DECAY

class TestCase4ExponentialDecay:
    """
    Integration tests for Test Case 4: Exponential Decay.
    
    From notebook:
        "Test to see how an exponentially decaying initial concentration 
        of the pollutant in time alters your results."
        
        C0 = 250
        decay_rate = 0.01
        boundary_condition=lambda t: C0 * np.exp(-decay_rate * t)
    """
    
    def test_exponential_decay_runs(self):
        """Check exponential decay simulation runs without error."""
        C0 = 250
        decay_rate = 0.01
        
        results = advection_solver(
            end_time=300, dt=10, length=20, dx=0.2, U=0.1, C0=C0,
            boundary_condition=lambda t: C0 * np.exp(-decay_rate * t)
        )
        
        assert results is not None
    
    def test_inlet_concentration_decreases(self):
        """
        Check inlet concentration decreases over time.
        
        At t=0: C = 250 * exp(0) = 250
        At t=300: C = 250 * exp(-3) ≈ 12.5
        """
        C0 = 250
        decay_rate = 0.01
        
        results = advection_solver(
            end_time=300, dt=10, length=20, dx=0.2, U=0.1, C0=C0,
            boundary_condition=lambda t: C0 * np.exp(-decay_rate * t)
        )
        
        inlet_start = results['concentrations'][0][0]
        inlet_end = results['concentrations'][-1][0]
        
        assert inlet_end < inlet_start
        assert inlet_end < 50  # Should be around 12.5
    
    def test_decay_rate_effect(self):
        """Check higher decay rate produces lower concentrations."""
        C0 = 250
        
        results_slow = advection_solver(
            end_time=300, dt=10, length=20, dx=0.2, U=0.1, C0=C0,
            boundary_condition=lambda t: C0 * np.exp(-0.001 * t)
        )
        results_fast = advection_solver(
            end_time=300, dt=10, length=20, dx=0.2, U=0.1, C0=C0,
            boundary_condition=lambda t: C0 * np.exp(-0.1 * t)
        )
        
        total_slow = np.sum(results_slow['concentrations'][-1])
        total_fast = np.sum(results_fast['concentrations'][-1])
        
        assert total_fast < total_slow



# TEST CASE 5: VARIABLE VELOCITY

class TestCase5VariableVelocity:
    """
    Integration tests for Test Case 5: Variable Velocity.
    
    From notebook:
        "Test to see how a variable stream velocity profile alters your results 
        (for example, add a 10% random perturbation to the constant velocity profile)."
        
        base_U = 0.1
        perturbation = 0.4  # ±40%
        variable_U = base_U * (1 + perturbation * (2*np.random.random(num_points) - 1))
    """
    
    def test_variable_velocity_runs(self):
        """Check variable velocity simulation runs without error."""
        length, dx = 20, 0.2
        num_points = int(length / dx) + 1
        base_U = 0.1
        perturbation = 0.4
        
        np.random.seed(42)
        variable_U = base_U * (1 + perturbation * (2*np.random.random(num_points) - 1))
        
        results = advection_solver(
            end_time=300, dt=10, length=20, dx=0.2,
            U=variable_U,
            C0=250
        )
        
        assert results is not None
        assert len(results['concentrations']) == 31
    
    def test_variable_differs_from_constant(self):
        """Check variable velocity produces different results than constant."""
        length, dx = 20, 0.2
        num_points = int(length / dx) + 1
        base_U = 0.1
        perturbation = 0.4
        
        # Constant velocity
        results_constant = advection_solver(
            end_time=300, dt=10, length=20, dx=0.2, U=0.1, C0=250
        )
        
        # Variable velocity (±40% as in notebook)
        np.random.seed(42)
        variable_U = base_U * (1 + perturbation * (2*np.random.random(num_points) - 1))
        results_variable = advection_solver(
            end_time=300, dt=10, length=20, dx=0.2, U=variable_U, C0=250
        )
        
        # Results should be different
        diff = np.sum(np.abs(
            results_constant['concentrations'][-1] - 
            results_variable['concentrations'][-1]
        ))
        
        assert diff > 0
    
    def test_different_perturbation_levels(self):
        """Check different perturbation levels all run successfully."""
        length, dx = 20, 0.2
        num_points = int(length / dx) + 1
        base_U = 0.1
        
        for perturbation in [0.1, 0.2, 0.4]:  # 10%, 20%, 40%
            np.random.seed(42)
            variable_U = base_U * (1 + perturbation * (2*np.random.random(num_points) - 1))
            
            results = advection_solver(
                end_time=300, dt=10, length=20, dx=0.2,
                U=variable_U, C0=250
            )
            
            assert results is not None


# FULL WORKFLOW TESTS

class TestFullWorkflow:
    """Integration tests for complete simulation workflow."""
    
    def test_results_dictionary_structure(self):
        """
        Check results dictionary has all required keys.
        
        From notebook return statement:
            return {
                'distance': distance,
                'time': time_array,
                'concentrations': all_concentrations,
                'parameters': {...}
            }
        """
        results = advection_solver(
            end_time=300, dt=10, length=20, dx=0.2, U=0.1, C0=250
        )
        
        assert 'distance' in results
        assert 'time' in results
        assert 'concentrations' in results
        assert 'parameters' in results
        
        assert 'end_time' in results['parameters']
        assert 'dt' in results['parameters']
        assert 'dx' in results['parameters']
        assert 'length' in results['parameters']
        assert 'num_points' in results['parameters']
    
    def test_reproducibility(self):
        """Check same inputs produce identical outputs."""
        results1 = advection_solver(
            end_time=300, dt=10, length=20, dx=0.2, U=0.1, C0=250
        )
        results2 = advection_solver(
            end_time=300, dt=10, length=20, dx=0.2, U=0.1, C0=250
        )
        
        assert np.allclose(
            results1['concentrations'][-1], 
            results2['concentrations'][-1]
        )
    
    def test_mass_increases_with_constant_source(self):
        """Check total mass increases when pollutant continuously enters."""
        results = advection_solver(
            end_time=300, dt=10, length=20, dx=0.2, U=0.1, C0=250
        )
        
        mass_initial = np.sum(results['concentrations'][0])
        mass_final = np.sum(results['concentrations'][-1])
        
        assert mass_final > mass_initial



# RUN TESTS

if __name__ == "__main__":
    import pytest
    pytest.main([__file__, "-v"])
