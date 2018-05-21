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

# Local imports
#from .utils import units


class Simulation(object):
    """Class for ignition delay simulations."""

    def __init__(self, idx, properties, model, path=''):
        """Initialize simulation case.

        :param idx: identifier number for this case
        :type idx: int
        :param properties: set of properties for this case
        :type properties: dict
        :param str model_file: Filename for Cantera-format model
        :param str path: Path for data file
        """
        self.idx = idx
        self.properties = properties

        self.gas = model


        self.time_end = 10 #This is just a filler idk how end time should actually be determined yet

        self.gas.TP = (self.properties['temperature'], self.properties['pressure']*float(ct.one_atm))
        self.gas.set_equivalence_ratio(self.properties['equivalence_ratio'],
                                       self.properties['fuel'],
                                       self.properties['oxidizer']
                                      )

        # Create non-interacting ``Reservoir`` on other side of ``Wall``
        env = ct.Reservoir(ct.Solution('air.xml'))

        # All reactors are ``IdealGasReactor`` objects
        self.reac = ct.IdealGasReactor(self.gas)
        self.wall = ct.Wall(self.reac, env, A=1.0, velocity=0)

        # Create ``ReactorNet`` newtork
        self.reac_net = ct.ReactorNet([self.reac])

        # Set file for later data file
        file_path = os.path.join(path, str(self.idx) + '.h5')
        self.save_file = file_path
        self.sample_points = []



        self.ignition_delay = 0.0

    def run_case(self, stop_at_ignition=False, restart=False):
        """Run simulation case set up ``setup_case``.

        :param bool stop_at_ignition: If ``True``, stop integration at ignition point.
        :param bool restart: If ``True``, skip if results file exists.
        """

        if restart and os.path.isfile(self.meta['save-file']):
            print('Skipped existing case ', self.meta['id'])
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

        print('Done with case ', self.idx)
        return self.ignition_delay

    def process_results(self):
        """Process integration results to sample data

        Returns arrays of sampled data.
        """

        number_sampled_points = 20

        # Load saved integration results
        with tables.open_file(self.save_file, 'r') as h5file:
            # Load Table with Group name simulation
            table = h5file.root.simulation

            temperatures = table.col('temperature')
            pressures = table.col('pressure')
            mass_fractions = table.col('mass_fractions')

        temperature_initial = temperatures[0]
        temperature_max = temperatures[len(temperatures)-1]
        delta = temperature_max - temperature_initial 

        multiplier = .05

        sampled_data = []

        # need to add processing to get the 20 data points here
        for temperature, pressure, mass_fraction in zip(temperatures, pressures, mass_fractions):
            # get temperature, pressure, mass fractions

            if temperature >= temperature_initial + (multiplier * delta):
                point_data = []
                point_data.append(temperature) 
                point_data.append(pressure) 
                point_data.append(mass_fraction)
                sampled_data.append(point_data)

                multiplier += .05
                if multiplier >= 1:
                    self.sample_points = sampled_data
                    return sampled_data

    def clean(self):
        """Delete HDF5 file with full integration data.
        """
        try:
            os.remove(self.save_file)
        except OSError:
            pass
