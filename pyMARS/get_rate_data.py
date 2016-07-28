import cantera as ct

def get_rates(mech_file, state_data):

    times=state_data.time
    temps=state_data.temp
    species_data = state_data.sp_data

    for i, t in enumerate(times):
        solution=ct.Solution(mech_file)
        solution.T = temps[i]












"for testing"
from autoignition_module import run_sim
A=run_sim('gri301.cti')

get_rates('gri301.cti', A)
