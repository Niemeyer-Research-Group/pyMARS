"""Autoignition simulation module

.. moduleauthor:: Kyle Niemeyer <kyle.niemeyer@gmail.com>
"""

# Standard libraries
import os
import logging

# Related modules
import numpy as np
import tables
import cantera as ct


ct.suppress_thermo_warnings()

class Simulation(object):
    """Class for general simulations

    Parameters
    ----------
    idx : int
        Identifer index for case
    conditions : InputIgnition or InputLaminarFlame
        Object with initial conditions for simulation
    model : str
        Filename for Cantera-format model to be used
    phase_name : str, optional
        Optional name for phase (e.g., 'gas'). 
    path : str, optional
        Path for location of output files
        
    """
    def __init__(self, sim_type, idx, conditions, model, phase_name='', path=''):
        self.idx = idx
        self.conditions = conditions
        self.model = model
        self.phase_name = phase_name
        self.path = path
        self.sim_type = sim_type

    def setup_case(self):
        """Initialize simulation case.
        """
        self.gas = ct.Solution(os.path.join(self.path,self.model), self.phase_name)

        if self.sim_type == 'ignition':
            # Default maximum number of steps
            self.max_steps = 10000
            if self.conditions.max_steps:
                self.max_steps = self.conditions.max_steps

            # By default, simulations will run to steady state, with the maximum number of steps
            # given by ``self.max_steps``. Alternatively, an end time (in seconds) can be
            # given in cases where something specific is needed (e.g., longer than normal)
            self.time_end = 0.0
            if self.conditions.end_time:
                self.time_end = self.conditions.end_time
        
            self.gas.TP = (
                self.conditions.temperature, self.conditions.pressure * ct.one_atm
                )
            # set initial composition using either equivalence ratio or general reactant composition
            if self.conditions.equivalence_ratio:
                self.gas.set_equivalence_ratio(
                    self.conditions.equivalence_ratio,
                    self.conditions.fuel,
                    self.conditions.oxidizer
                    )
            else:
                if self.conditions.composition_type == 'mole':
                    self.gas.TPX = (
                        self.conditions.temperature, self.conditions.pressure * ct.one_atm, 
                        self.conditions.reactants
                        )
                else:
                    self.gas.TPY = (
                        self.conditions.temperature, self.conditions.pressure * ct.one_atm, 
                        self.conditions.reactants
                        )

            if self.conditions.kind == 'constant pressure':
                self.reac = ct.IdealGasConstPressureReactor(self.gas)
            else:
                self.reac = ct.IdealGasReactor(self.gas)

            # Create ``ReactorNet``
            self.sim = ct.ReactorNet([self.reac])

            # Set file for later data file
            self.save_file = os.path.join(self.path, self.sim_type + '_' + str(self.idx) + '.h5')
            self.sample_points = []

            self.ignition_delay = 0.0

        elif self.sim_type == 'flame':
            self.gas.TP = (
                self.conditions.temperature, self.conditions.pressure * ct.one_atm
                )
            # set initial composition using either equivalence ratio or general reactant composition
            if self.conditions.equivalence_ratio:
                self.gas.set_equivalence_ratio(
                    self.conditions.equivalence_ratio,
                    self.conditions.fuel,
                    self.conditions.oxidizer
                    )
            else:
                if self.conditions.composition_type == 'mole':
                    self.gas.TPX = (
                        self.conditions.temperature, self.conditions.pressure * ct.one_atm, 
                        self.conditions.reactants
                        )
                else:
                    self.gas.TPY = (
                        self.conditions.temperature, self.conditions.pressure * ct.one_atm, 
                        self.conditions.reactants
                        )

            if self.conditions.kind == 'premixed':
                self.gas.transport_model = self.conditions.transport
                self.sim = ct.FreeFlame(self.gas, width=self.conditions.width)
                self.sim.set_refine_criteria(ratio=self.conditions.refine_ratio,slope=self.conditions.refine_slope,
                                             curve=self.conditions.refine_curve,prune=self.conditions.refine_prune)

            # Set file for later data file
            self.save_file = os.path.join(self.path, self.sim_type + '_' + str(self.idx) + '.h5')
            self.sample_points = []

            self.ignition_delay = 0.0


    def run_case(self, stop_at_ignition=False, restart=False):
        """Run simulation case set up ``setup_case``.

        If no end time is specified for the integration, the function integrates
        to steady state (or a maximum of 10,000 steps, by default). This is done
        by checking whether the system state changes below a certain threshold,
        with the residual computed using feature checking. This is blatantly stolen
        from Cantera's :meth:`cantera.ReactorNet.advance_to_steady_state` method.

        Parameters
        ----------
        stop_at_ignition : bool
            If ``True``, stop integration at ignition point, don't save data.
        restart : bool
            If ``True``, skip if results file exists.
        
        Returns
        -------
        self.ignition_delay : float
            Computed ignition delay in seconds

        """

        if restart and os.path.isfile(self.save_file):
            print('Skipped existing case ', self.idx)
            return

        # Save simulation results in hdf5 table format.
        table_def = {'metric': tables.Float64Col(pos=0),
                     'temperature': tables.Float64Col(pos=1),
                     'pressure': tables.Float64Col(pos=2),
                     'mass_fractions': tables.Float64Col(
                          shape=(self.gas.n_species), pos=3
                          ),
                     }

        with tables.open_file(self.save_file, mode='w',
                              title=str(self.idx)
                              ) as h5file:

            table = h5file.create_table(where=h5file.root,
                                        name='simulation',
                                        description=table_def
                                        )
                                        
            # Row instance to save timestep information to
            fluid_parcel = table.row

            if self.sim_type == 'ignition':
                # Save initial conditions
                fluid_parcel['metric'] = self.sim.time
                fluid_parcel['temperature'] = self.reac.T
                fluid_parcel['pressure'] = self.reac.thermo.P
                fluid_parcel['mass_fractions'] = self.reac.Y
                # Add ``fluid_parcel`` to table
                fluid_parcel.append()

                ignition_flag = False

                # Main time integration loop
                if self.time_end:
                    # if end time specified, continue integration until reaching that time
                    while self.sim.time < self.time_end:
                        self.sim.step()

                        # Save new fluid_parcel information
                        fluid_parcel['metric'] = self.sim.time
                        fluid_parcel['temperature'] = self.reac.T
                        fluid_parcel['pressure'] = self.reac.thermo.P
                        fluid_parcel['mass_fractions'] = self.reac.Y

                        if self.reac.T >= self.conditions.temperature + 400.0 and not ignition_flag:
                            self.ignition_delay = self.sim.time
                            ignition_flag = True

                            if stop_at_ignition:
                                break

                        # Add ``fluid_parcel`` to table
                        fluid_parcel.append()
                    
                else:
                    # otherwise, integrate until steady state, or maximum number of steps reached
                    self.sim.reinitialize()
                    max_state_values = self.sim.get_state()
                    residual_threshold = 10. * self.sim.rtol
                    absolute_tolerance = self.sim.atol

                    for step in range(self.max_steps):
                        previous_state = self.sim.get_state()

                        self.sim.step()

                        # Save new fluid_parcel information
                        fluid_parcel['metric'] = self.sim.time
                        fluid_parcel['temperature'] = self.reac.T
                        fluid_parcel['pressure'] = self.reac.thermo.P
                        fluid_parcel['mass_fractions'] = self.reac.Y

                        if self.reac.T >= self.conditions.temperature + 400.0 and not ignition_flag:
                            self.ignition_delay = self.sim.time
                            ignition_flag = True

                            if stop_at_ignition:
                                break

                        # Add ``fluid_parcel`` to table
                        fluid_parcel.append()

                        state = self.sim.get_state()
                        max_state_values = np.maximum(max_state_values, state)
                        residual = np.linalg.norm(
                            (state - previous_state) / (max_state_values + absolute_tolerance)
                            ) / np.sqrt(self.sim.n_vars)

                        if residual < residual_threshold:
                            break
                    
                    if step == self.max_steps - 1:
                        logging.warning(
                            'Warning: Maximum number of steps reached before '
                            f'convergence for ignition case {self.idx}. '
                            f'Do NOT use the final mechanism if this warning '
                            f'appears at the end right before generating the final mechanism.'
                            )
                    
                if not ignition_flag:
                    logging.error(f'No ignition detected for ignition case {self.idx}')
                    raise RuntimeError(f'No ignition detected for ignition case {self.idx}')
                
                out_metric = self.ignition_delay
                    
            elif self.sim_type == 'flame':
                self.sim.solve(loglevel=0)

                # Extract only points from Tu + 0.05*(Tb-Tu) < T < Tb - 0.05*(Tb-Tu)
                T_low = self.sim.T[0] + 0.05*(self.sim.T[-1] - self.sim.T[0])
                T_high = self.sim.T[-1] - 0.05*(self.sim.T[-1] - self.sim.T[0])

                for i in range(len(self.sim.grid)):
                    if (self.sim.T[i] > T_low) and (self.sim.T[i] < T_high):
                        fluid_parcel['metric'] = self.sim.velocity[0] # save the laminar flame speed. A bit redundant as this saves it for every parcel, but fits well with current architecture so we don't need to refactor
                        fluid_parcel['temperature'] = self.sim.T[i]
                        fluid_parcel['pressure'] = self.sim.P
                        fluid_parcel['mass_fractions'] = self.sim.Y[:,i]

                        fluid_parcel.append()
                        
                out_metric = self.sim.velocity[0]

            # Write ``table`` to disk
            table.flush()

        return out_metric

    def process_results(self, skip_data=False):
        """Process integration results to sample data

        Parameters
        ----------
        skip_data : bool
            Flag to skip sampling thermochemical data

        Returns
        -------
        tuple of float, numpy.ndarray or float
            Ignition delay, or ignition delay and sampled data

        """
        delta = 0.05
        deltas = np.arange(delta, 1 + delta, delta)

        # Load saved integration results
        self.save_file = os.path.join(self.path, self.sim_type + '_' + str(self.idx) + '.h5')
        if self.sim_type == 'ignition':
            with tables.open_file(self.save_file, 'r') as h5file:
                # Load Table with Group name simulation
                table = h5file.root.simulation

                times = table.col('metric')
                temperatures = table.col('temperature')
                pressures = table.col('pressure')
                mass_fractions = table.col('mass_fractions')

            temperature_initial = temperatures[0]
            temperature_max = temperatures[len(temperatures)-1]
            temperature_diff = temperature_max - temperature_initial 

            sampled_data = np.zeros((len(deltas), 2 + mass_fractions.shape[1]))

            # need to add processing to get the 20 data points here
            self.ignition_delay = 0.0
            ignition_flag = False
            idx = 0
            for time, temp, pres, mass in zip(
                times, temperatures, pressures, mass_fractions
                ):
                if temp >= temperature_initial + 400.0 and not ignition_flag:
                    self.ignition_delay = time
                    ignition_flag = True
                    if skip_data:
                        return self.ignition_delay, None, self.conditions.target_weights
                
                if temp >= temperature_initial + (deltas[idx] * temperature_diff):
                    sampled_data[idx, 0:2] = [temp, pres]
                    sampled_data[idx, 2:] = mass

                    idx += 1
                    if idx == 20:
                        self.sampled_data = sampled_data
                        return self.ignition_delay, sampled_data, np.tile(self.conditions.target_weights, (idx, 1))
        elif self.sim_type == 'flame':
            with tables.open_file(self.save_file, 'r') as h5file:
                # Load Table with Group name simulation
                table = h5file.root.simulation

                velocities = table.col('metric')
                temperatures = table.col('temperature')
                pressures = table.col('pressure')
                mass_fractions = table.col('mass_fractions')

            flame_speed = velocities[0]
            temperature_initial = temperatures[0]
            temperature_max = temperatures[len(temperatures)-1]
            temperature_diff = temperature_max - temperature_initial 

            sampled_data = np.zeros((len(deltas), 2 + mass_fractions.shape[1]))

            # need to add processing to get the 20 data points here
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
                        return flame_speed, sampled_data, np.tile(self.conditions.target_weights, (idx,1))

    def clean(self):
        """Delete HDF5 save file
        """
        try:
            os.remove(self.save_file)
        except OSError:
            pass
