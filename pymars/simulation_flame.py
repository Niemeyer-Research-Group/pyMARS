"""Flame Speed Simulation module

.. moduleauthor:: Cailin Moore
"""

# Standard libraries
import os
import logging

# Related modules
import numpy as np
import cantera as ct

ct.suppress_thermo_warnings()

class FlameSimulation(object):
    """Class for one dimensional flame simulations
    Parameters
    ----------
    idx : int
        Identifer index for case
    properties : InputLaminarFlame
        Object with initial conditions for simulation
    model : str
        Filename for Cantera-format model to be used
    phase_name : str, optional
        Optional name for phase to load from CTI file (e.g., 'gas'). 
    path : str, optional
        Path for location of output files
        
    """
    def __init__(self, idx, properties, model, phase_name='', path=''):
        self.idx = idx
        self.properties = properties
        self.model = model
        self.phase_name = phase_name
        self.path = path
    
    def setup_case(self):
        """Initialize simulation case.
        """
        self.gas = ct.Solution(self.model, self.phase_name)
    
        self.gas.TP = (
            self.properties.temperature, self.properties.pressure * ct.one_atm
            )
        # set initial composition using either equivalence ratio or general reactant composition
        if self.properties.equivalence_ratio:
            self.gas.set_equivalence_ratio(
                self.properties.equivalence_ratio,
                self.properties.fuel,
                self.properties.oxidizer
                )
        else:
            if self.properties.composition_type == 'mole':
                self.gas.TPX = (
                    self.properties.temperature, self.properties.pressure * ct.one_atm, 
                    self.properties.reactants
                    )
            else:
                self.gas.TPY = (
                    self.properties.temperature, self.properties.pressure * ct.one_atm, 
                    self.properties.reactants
                    )

        if self.properties.kind == 'constant pressure':
            self.reac = ct.IdealGasConstPressureReactor(self.gas)
        else:
            self.reac = ct.IdealGasReactor(self.gas)

        self.sim = ct.ReactorNet([self.reac])
        
        # implementation of allowing option for width to be input - would this be in properties 
        #if not self.properties.width:
            #default flame width to 25 cm 
        #self.properties.width = 0.014

        # Create the flame object
        self.flame = ct.FreeFlame(self.gas, width=self.properties.width)

        # Define tolerances for the solver
        self.flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1)
        
        # Set file for later data file 
        self.save_file = os.path.join(self.path, str(self.idx) + '.h5')
        self.sample_points = []

        self.Su0 = 0.0
     
    def run_case(self, restart=False):
        """Run simulation case set up ``setup_case``.
        Runs canteras flame.solve funciton on given parameters. 
        Parameters
        ----------
        stop_at_flame : remove? no gridpoint calculation needed 
        restart : bool
            If ``True``, skip if results file exists.
        
        Returns
        -------
        self. : float
            Computed ignition delay in seconds
        """

        if restart and os.path.isfile(self.save_file):
            print('Skipped existing case ', self.idx)
            return
        
        try:
            self.flame.solve(loglevel=0, refine_grid=True, auto=True)
            self.Su0 = self.flame.u[0]
        except:
            logging.error(f'No flame detected for laminar flame case {self.idx}')
            raise RuntimeError(f'No flame detected for laminar flame case {self.idx}')

        return self.Su0

    def calculate_flamespeed(self):
        """Run simulation case set up ``setup_case``, just for flame speed.
        """    
        try:
            self.run_case()
        except:
            logging.error(f'No flame detected for laminar flame case {self.idx}')
            raise RuntimeError(f'No flame detected for laminar flame case {self.idx}')

        return self.Su0
    
    def process_results(self, skip_data=False):
        """Process integration results to sample data
        Parameters
        ----------
        skip_data : bool
            Flag to skip sampling thermochemical data
        Returns
        -------
        tuple of float, numpy.ndarray or float
            Flame speed, or flame speed and sampled data
        """

        self.Su0 = self.calculate_flamespeed()
        if skip_data == True:
            return self.Su0

        delta = 0.05
        deltas = np.arange(delta, 1 + delta, delta)
        
        temperatures = self.flame.T
        pressures = self.flame.P
        mass_fractions = self.flame.Y
        flame_speeds = self.Su0

        temperature_initial = temperatures[0]
        temperature_max = temperatures[len(temperatures)-1]
        temperature_diff = temperature_max - temperature_initial 

        sampled_data = np.zeros((len(deltas), 4 + mass_fractions.shape[1]))

        #processing to get the 20 data points here
        idx = 0
        for temp, pres, mass in zip(
            temperatures, pressures, mass_fractions
            ):
            
            if temp >= temperature_initial + (deltas[idx] * temperature_diff):
                sampled_data[idx, 0:2] = [temp, pres]
                sampled_data[idx, 2:] = mass

                idx += 1
                if idx == 20:
                    self.sampled_data = sampled_data
                    return flame_speeds, sampled_data
        #processing to get 20 data points for additional parameters with different sizes? 

    def clean(self):
        """Delete HDF5 file with full integration data.
        """
        try:
            os.remove(self.save_file)
        except OSError:
            pass        
    