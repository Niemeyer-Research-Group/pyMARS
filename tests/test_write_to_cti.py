#tests the cti output file
import cantera as ct
import os

def test_solution_object():
    input_file=os.path.abspath('Input_Data_Files/' + 'gri30.cti')
    output_file_name=os.path.abspath('Output_Data_Files/'+'trimmed_gri30.cti')
    Solution_object=ct.Solution(output_file_name)
    assert type(Solution_object).__name__ == 'Solution'
