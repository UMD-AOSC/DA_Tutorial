from class_state_vector import state_vector
from class_da_system import da_system
import numpy as np
from sys import argv
import matplotlib.pyplot as plt

nature_file ='x_nature.pkl'
nature = state_vector()
nature = nature.load(nature_file)

freerun_file = 'x_freerun.pkl'
freerun = state_vector()
freerun = freerun.load(freerun_file)

method = argv[1]
analysis_file = 'x_analysis_{method}.pkl'.format(method=method)
das = da_system()
das = das.load(analysis_file)
analysis = das.getStateVector()

plt.plot(nature.getTimes(),
         np.linalg.norm(freerun.getTrajectory() - nature.getTrajectory(),
                        axis=1), label='Free run')
plt.plot(nature.getTimes(),
         np.linalg.norm(analysis.getTrajectory() - nature.getTrajectory(),
                        axis=1), label='Analysis ({method})'.format(method=method))
plt.legend()
plt.xlabel('Time')
plt.ylabel('Error', rotation='horizontal', labelpad=20)
plt.tight_layout()
plt.savefig("img/tmp.png")
