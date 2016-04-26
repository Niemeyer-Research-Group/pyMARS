#tests the cti output file
import cantera as ct
import os

def test_solution_object():
    input_file=os.path.abspath('Input_Data_Files/' + 'gri30_converted.cti')
    output_file_name=os.path.abspath('Output_Data_Files/'+'gri30_converted.cti')
    Solution_object=ct.Solution(output_file_name)
    assert type(Solution_object).__name__ == 'Solution'
