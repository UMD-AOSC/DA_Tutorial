from class_state_vector import state_vector
from class_da_system import da_system
import numpy as np
from sys import argv
import matplotlib.pyplot as plt

def main():
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
    plot_rmse_all(nature, freerun, analysis, method, np.s_[:, :],
        "img/%s/rmse_all.pdf" % method)
    plot_rmse_all(nature, freerun, analysis, method, np.s_[:, 0:10],
        "img/%s/rmse_atmos_psi.pdf" % method)
    plot_rmse_all(nature, freerun, analysis, method, np.s_[:, 10:20],
        "img/%s/rmse_atmos_temp.pdf" % method)
    plot_rmse_all(nature, freerun, analysis, method, np.s_[:, 20:28],
        "img/%s/rmse_ocean_psi.pdf" % method)
    plot_rmse_all(nature, freerun, analysis, method, np.s_[:, 28:36],
        "img/%s/rmse_ocean_temp.pdf" % method)

def plot_rmse_all(nature, freerun, analysis, method, slice, img_name):
    plt.plot(nature.getTimes(),
             np.linalg.norm(freerun.getTrajectory()[slice] - nature.getTrajectory()[slice],
                            axis=1), label='Free run')
    plt.plot(nature.getTimes(),
             np.linalg.norm(analysis.getTrajectory()[slice] - nature.getTrajectory()[slice],
                            axis=1), label='Analysis ({method})'.format(method=method))
    if analysis.getEnsembleTrajectory() is not None:
        ptb = analysis.getEnsembleTrajectory()[:, :, :] - analysis.getTrajectory()[:, :, np.newaxis]
        edim = ptb.shape[2]
        sprd = (np.sum(ptb ** 2, axis=2) / (edim - 1.0)) ** 0.5
        plt.plot(nature.getTimes(), np.linalg.norm(sprd[slice], axis=1),
            label='Analysis spread ({method})'.format(method=method))
    plt.legend()
    plt.yscale("log")
    plt.xlabel('Time')
    plt.ylabel('Error', rotation='horizontal', labelpad=20)
    plt.title(img_name)
    plt.tight_layout()
    plt.savefig(img_name)
    plt.close()

main()
