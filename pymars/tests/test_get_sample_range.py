import numpy as np
import matplotlib.pyplot as plt
from pyMARS import get_sample_range as gsr


def log_data(midpoint):
    x = np.linspace(-2,8.0,10000)
    y = []
    x_0 = midpoint
    k = 1.0
    L = 1.0
    for xi in x:
        y.append( L/(1+np.exp(-k*(xi-x_0))))

    middle_index = np.argmin(np.abs(x-x_0))
    sdata = np.zeros((len(x), len(x)))
    production_data = np.zeros((len(x), len(x)))

    #plt.plot(x,y)
    #plt.show()
    return [x,y, sdata, production_data, middle_index, x_0]

#generate test data
sample_1 = log_data(1.0)
sample_2 = log_data(3.0)

#call get sample range function for testing
range_1 = gsr.get_range(sample_1[0], sample_1[1], sample_1[2], sample_1[3])
range_2 = gsr.get_range(sample_2[0], sample_2[1], sample_2[2], sample_2[3])

#check ign delay  (tau of 1.0 and 3.0 expected, and .1% error allowed)
assert abs(range_1.tau - sample_1[0][sample_1[4]])/abs(sample_1[0][sample_1[4]])*100 < .1
assert abs(range_2.tau - sample_2[0][sample_2[4]])/abs(sample_2[0][sample_2[4]])*100 < .1

#check sample range size
assert len(range_1.times) and len(range_1.temps) == 40
assert len(range_2.times) and len(range_2.temps) == 40
