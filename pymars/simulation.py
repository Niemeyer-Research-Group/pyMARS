"""Autoignition simulation module

.. moduleauthor:: Kyle Niemeyer <kyle.niemeyer@gmail.com>
"""

# Standard libraries
import os
import numpy as np

# Related modules
try:
    import cantera as ct
    ct.suppress_thermo_warnings()
except ImportError:
    print("Error: Cantera must be installed.")
    raise

try:
    import tables
except ImportError:
    print('PyTables must be installed')
    raise

class Simulation(object):
    """Class for ignition delay simulations

    Parameters
    ----------
    idx : int
        Identifer index for case
    properties : dict
        Dictionary with all properties needed
    """
    def __init__(self, idx, properties, model, path=''):
        self.idx = idx
        self.properties = properties
        self.model = model
        self.path = path

    def setup_case(self):
        """Initialize simulation case.
        """
        self.gas = ct.Solution(self.model)

        self.time_end = 10.0
        if 'end-time' in self.properties:
            self.time_end = self.properties['end-time']

        self.gas.TP = (
            self.properties['temperature'], self.properties['pressure'] * ct.one_atm
            )
        self.gas.set_equivalence_ratio(self.properties['equivalence-ratio'],
                                       self.properties['fuel'],
                                       self.properties['oxidizer']
                                      )

        if self.properties['kind'] == 'constant pressure':
            self.reac = ct.IdealGasConstPressureReactor(self.gas)
        else:
            self.reac = ct.IdealGasReactor(self.gas)

        # Create ``ReactorNet`` newtork
        self.reac_net = ct.ReactorNet([self.reac])

        # Set file for later data file
        self.save_file = os.path.join(self.path, str(self.idx) + '.h5')
        self.sample_points = []

        self.ignition_delay = 0.0

    def run_case(self, stop_at_ignition=False, restart=False):
        """Run simulation case set up ``setup_case``.

        :param bool stop_at_ignition: If ``True``, stop integration at ignition point.
        :param bool restart: If ``True``, skip if results file exists.
        """

        if restart and os.path.isfile(self.save_file):
            print('Skipped existing case ', self.idx)
            return

        # Save simulation results in hdf5 table format.
        table_def = {'time': tables.Float64Col(pos=0),
                     'temperature': tables.Float64Col(pos=1),
                     'pressure': tables.Float64Col(pos=2),
                     'mass_fractions': tables.Float64Col(
                          shape=(self.reac.thermo.n_species), pos=3
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
            timestep = table.row
            # Save initial conditions
            timestep['time'] = self.reac_net.time
            timestep['temperature'] = self.reac.T
            timestep['pressure'] = self.reac.thermo.P
            timestep['mass_fractions'] = self.reac.Y
            # Add ``timestep`` to table
            timestep.append()

            ignition_flag = False

            # Main time integration loop; continue integration while time of
            # the ``ReactorNet`` is less than specified end time.
            while self.reac_net.time < self.time_end:
                self.reac_net.step()

                # Save new timestep information
                timestep['time'] = self.reac_net.time
                timestep['temperature'] = self.reac.T
                timestep['pressure'] = self.reac.thermo.P
                timestep['mass_fractions'] = self.reac.Y

                if self.reac.T > self.properties['temperature'] + 400.0 and not ignition_flag:
                    self.ignition_delay = self.reac_net.time
                    ignition_flag = True

                    if stop_at_ignition:
                        continue

                # Add ``timestep`` to table
                timestep.append()

            # Write ``table`` to disk
            table.flush()

        return self.ignition_delay

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
        self.save_file = os.path.join(self.path, str(self.idx) + '.h5')
        with tables.open_file(self.save_file, 'r') as h5file:
            # Load Table with Group name simulation
            table = h5file.root.simulation

            times = table.col('time')
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
            if temp > temperature_initial + 400.0 and not ignition_flag:
                    self.ignition_delay = time
                    ignition_flag = True
                    if skip_data:
                        return self.ignition_delay
            
            if temp >= temperature_initial + (deltas[idx] * temperature_diff):
                sampled_data[idx, 0:2] = [temp, pres]
                sampled_data[idx, 2:] = mass

                idx += 1
                if idx == 20:
                    self.sampled_data = sampled_data
                    return self.ignition_delay, sampled_data

    def clean(self):
        """Delete HDF5 file with full integration data.
        """
        try:
            os.remove(self.save_file)
        except OSError:
            pass
