"""
Unit Tests for Advection Solver
===============================
Tests individual functions to ensure they work correctly.

These tests verify the advection_solver function from the main notebook.

Run with: python -m pytest tests/test_unit.py -v
Or in Colab: !pytest tests/test_unit.py -v
"""

import numpy as np
import sys
sys.path.append('../src')
from advection_solver import advection_solver



# TESTS FOR SPATIAL GRID

class TestSpatialGrid:
    """Tests for spatial grid creation."""
    
    def test_grid_size(self):
        """
        Check correct number of grid points.
        
        From notebook: num_points = int(length / dx) + 1
        With length=20, dx=0.2: num_points = 20/0.2 + 1 = 101
        """
        results = advection_solver(end_time=10, dt=10, length=20, dx=0.2, U=0.1, C0=250)
        assert len(results['distance']) == 101
    
    def test_grid_start_end(self):
        """
        Check grid starts at 0 and ends at length.
        
        From notebook: distance = np.linspace(0, length, num_points)
        """
        results = advection_solver(end_time=10, dt=10, length=20, dx=0.2, U=0.1, C0=250)
        assert results['distance'][0] == 0
        assert results['distance'][-1] == 20
    
    def test_grid_spacing(self):
        """Check grid has correct spacing between points."""
        results = advection_solver(end_time=10, dt=10, length=20, dx=0.5, U=0.1, C0=250)
        spacing = results['distance'][1] - results['distance'][0]
        assert abs(spacing - 0.5) < 0.001



# TESTS FOR TIME ARRAY

class TestTimeArray:
    """Tests for time array creation."""
    
    def test_time_frames(self):
        """
        Check correct number of time frames.
        
        From notebook: num_frames = int(end_time / dt) + 1
        With end_time=300, dt=10: num_frames = 300/10 + 1 = 31
        """
        results = advection_solver(end_time=300, dt=10, length=20, dx=0.2, U=0.1, C0=250)
        assert len(results['time']) == 31
    
    def test_time_start_end(self):
        """
        Check time starts at 0 and ends at end_time.
        
        From notebook: time_array = np.linspace(0, end_time, num_frames)
        """
        results = advection_solver(end_time=300, dt=10, length=20, dx=0.2, U=0.1, C0=250)
        assert results['time'][0] == 0
        assert results['time'][-1] == 300



# TESTS FOR INITIAL CONDITIONS


class TestInitialConditions:
    """Tests for initial conditions setup."""
    
    def test_default_initial_conditions(self):
        """
        Check default IC: C0 at x=0, zeros elsewhere.
        
        From notebook:
            initial_conditions = np.zeros(num_points)
            initial_conditions[0] = C0
        """
        results = advection_solver(end_time=10, dt=10, length=20, dx=0.2, U=0.1, C0=250)
        initial = results['concentrations'][0]
        
        assert initial[0] == 250      # C0 at inlet (x=0)
        assert all(initial[1:] == 0)  # Zero everywhere else
    
    def test_custom_initial_conditions(self):
        """
        Check custom initial conditions are used.
        
        From notebook Test Case 2: initial_conditions=initial_cond
        """
        custom_ic = np.ones(101) * 100  # 100 everywhere
        
        results = advection_solver(end_time=10, dt=10, length=20, dx=0.2, U=0.1, 
                                   initial_conditions=custom_ic)
        initial = results['concentrations'][0]
        
        assert all(initial == 100)



# TESTS FOR BOUNDARY CONDITIONS

class TestBoundaryConditions:
    """Tests for boundary condition handling."""
    
    def test_constant_boundary(self):
        """
        Check constant boundary condition stays at C0.
        
        From notebook: bc_func = lambda t: C0
        """
        results = advection_solver(end_time=100, dt=10, length=20, dx=0.2, U=0.1, C0=250)
        
        for conc in results['concentrations']:
            assert conc[0] == 250
    
    def test_function_boundary(self):
        """
        Check time-varying boundary condition works.
        
        From notebook: elif callable(boundary_condition): bc_func = boundary_condition
        """
        bc_func = lambda t: 100
        
        results = advection_solver(end_time=100, dt=10, length=20, dx=0.2, U=0.1, 
                                   C0=250, boundary_condition=bc_func)
        
        assert results['concentrations'][1][0] == 100
    
    def test_exponential_decay_boundary(self):
        """
        Check exponential decay boundary decreases over time.
        
        From notebook Test Case 4:
            boundary_condition=lambda t: C0 * np.exp(-decay_rate * t)
        """
        C0 = 250
        decay_rate = 0.01
        bc_func = lambda t: C0 * np.exp(-decay_rate * t)
        
        results = advection_solver(end_time=300, dt=10, length=20, dx=0.2, U=0.1,
                                   C0=C0, boundary_condition=bc_func)
        
        inlet_start = results['concentrations'][0][0]
        inlet_end = results['concentrations'][-1][0]
        
        assert inlet_end < inlet_start


# TESTS FOR CONCENTRATION VALUES

class TestConcentrationValues:
    """Tests for concentration calculations."""
    
    def test_concentration_non_negative(self):
        """Check concentration never goes negative."""
        results = advection_solver(end_time=300, dt=10, length=20, dx=0.2, U=0.1, C0=250)
        
        for conc in results['concentrations']:
            assert all(conc >= 0)
    
    def test_concentration_bounded(self):
        """Check concentration doesn't exceed C0."""
        results = advection_solver(end_time=300, dt=10, length=20, dx=0.2, U=0.1, C0=250)
        
        for conc in results['concentrations']:
            assert all(conc <= 250 * 1.01)  # Allow 1% tolerance
    
    def test_concentration_spreads_downstream(self):
        """Check pollutant moves downstream over time."""
        results = advection_solver(end_time=300, dt=10, length=20, dx=0.2, U=0.1, C0=250)
        
        initial = results['concentrations'][0]
        final = results['concentrations'][-1]
        
        initial_nonzero = np.sum(initial > 0)
        final_nonzero = np.sum(final > 0)
        
        assert final_nonzero > initial_nonzero



# TESTS FOR PARAMETERS STORAGE

class TestParameters:
    """Tests for parameter storage in results dictionary."""
    
    def test_parameters_stored(self):
        """
        Check all parameters are stored in results.
        
        From notebook return statement:
            'parameters': {'end_time': end_time, 'dt': dt, 'dx': dx,
                          'length': length, 'num_points': num_points}
        """
        results = advection_solver(end_time=300, dt=10, length=20, dx=0.2, U=0.1, C0=250)
        
        assert results['parameters']['end_time'] == 300
        assert results['parameters']['dt'] == 10
        assert results['parameters']['length'] == 20
        assert results['parameters']['dx'] == 0.2
    
    def test_num_points_stored(self):
        """Check num_points is correctly calculated and stored."""
        results = advection_solver(end_time=10, dt=10, length=20, dx=0.2, U=0.1, C0=250)
        
        assert results['parameters']['num_points'] == 101


# TESTS FOR VARIABLE VELOCITY

class TestVariableVelocity:
    """Tests for spatially varying velocity (Test Case 5)."""
    
    def test_array_velocity_accepted(self):
        """
        Check solver accepts velocity as an array.
        
        From notebook: if np.isscalar(A) converts to arrays
        """
        num_points = 101
        variable_U = np.ones(num_points) * 0.1
        
        results = advection_solver(end_time=100, dt=10, length=20, dx=0.2, 
                                   U=variable_U, C0=250)
        
        assert len(results['concentrations']) > 1
    
    def test_variable_velocity_with_perturbation(self):
        """
        Check solver runs with random perturbation velocity.
        
        From notebook Test Case 5:
            variable_U = base_U * (1 + perturbation * (2*np.random.random(num_points) - 1))
        """
        length, dx = 20, 0.2
        num_points = int(length / dx) + 1
        base_U = 0.1
        perturbation = 0.4
        
        np.random.seed(42)
        variable_U = base_U * (1 + perturbation * (2*np.random.random(num_points) - 1))
        
        results = advection_solver(end_time=300, dt=10, length=20, dx=0.2, 
                                   U=variable_U, C0=250)
        
        assert len(results['concentrations']) == 31



# RUN TESTS

if __name__ == "__main__":
    import pytest
    pytest.main([__file__, "-v"])
