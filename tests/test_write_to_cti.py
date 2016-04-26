#tests the cti output file
import cantera as ct

def test_solution_object():
    Solution_object=ct.Solution('Output_Data_Files/gri30_converted.cti')
    assert type(Solution_object).__name__ == 'Solution'
