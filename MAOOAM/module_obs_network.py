import numpy as np
import pickle
import matplotlib
matplotlib.use("TKAgg")
import matplotlib.pyplot as plt

NDIM = 36

def get_grid_val(waves, x, y, elem):
    def get_atm(is_temp):
        types = ["A", "K", "L", "A", "K", "L", "K", "L", "K", "L"]
        hs = [0, 1, 1, 0, 1, 1, 2, 2, 2, 2]
        ps = [1, 1, 1, 2, 2, 2, 1, 1, 2, 2]
        gridval = 0.0
        for j in range(na):
            j_all = j + na if is_temp else j
            if types[j] == "A":
                gridval = gridval + waves[j_all] * np.sqrt(2.0) * np.cos(ps[j] * y)
            elif types[j] == "K":
                gridval = gridval + waves[j_all] * 2.0 * np.cos(hs[j] * n * x) * np.sin(ps[j] * y)
            else:
                gridval = gridval + waves[j_all] * 2.0 * np.sin(hs[j] * n * x) * np.sin(ps[j] * y)
        if is_temp:
            gridval *= (f0 ** 2 * L ** 2) / R
        else:
            gridval *= L ** 2 * f0
        return gridval

    def get_ocn(is_temp):
        hos = [1, 1, 1, 1, 2, 2, 2, 2]
        pos = [1, 2, 3, 4, 1, 2, 3, 4]
        gridval = 0.0
        for j in range(no):
            j_all = j + (na * 2 + no) if is_temp else j + na * 2
            gridval = gridval + waves[j_all] * 2.0 * np.sin(0.5 * hos[j] * n * x) * np.sin(pos[j] * y)
        if is_temp:
            gridval *= (f0 ** 2 * L ** 2) / R
        else:
            gridval *= L ** 2 * f0
        return gridval

    assert waves.__class__ == np.ndarray
    assert waves.shape == (NDIM,)
    assert y.__class__ in [float, np.float32, np.float64]
    assert x.__class__ in [float, np.float32, np.float64]

    f0 = 1.032e-4
    n = 1.5
    na = 10
    no = 8
    R = 287.0
    L = 5000000.0 / np.pi

    if elem == "a_psi":
        return get_atm(False)
    elif elem == "a_tmp":
        return get_atm(True)
    elif elem == "o_psi":
        return get_ocn(False)
    elif elem == "o_tmp":
        return get_ocn(True)
    else:
        raise Exception("get_grid_val overflow. Element %s not found." % elem)

def __test_get_grid_val():
    state = model_state_exsample()
    x = 1.2 * np.pi / n  # 0.0 <= x <= 2.0 * pi / n
    y = 0.8 * np.pi      # 0.0 <= y <= pi

    # get_grid_val() returns one of four variables
    # {atmosphere|ocean} x {streamfunction|temperature} at point (x, y)
    # unit: [m^2/s] for streamfunction and [K] for temperature
    a_psi = get_grid_val(state, x, y, "a_psi")
    a_tmp = get_grid_val(state, x, y, "a_tmp")
    o_psi = get_grid_val(state, x, y, "o_psi")
    o_tmp = get_grid_val(state, x, y, "o_tmp")
    print(a_psi, a_tmp, o_psi, o_tmp)

def __get_obs_grid_atmos():
    n = 1.5
    xmax = 2.0 * np.pi / n
    ymax = np.pi
    nxobs = 5
    nyobs = 2
    x1d = np.linspace(0, xmax, nxobs, endpoint=False)
    y1d = np.linspace(0, ymax, nyobs, endpoint=False) + ymax / nyobs * 0.25
    x2d, y2d = np.meshgrid(x1d, y1d)
    return x2d, y2d

def __get_obs_grid_ocean():
    n = 1.5
    xmax = 2.0 * np.pi / n
    ymax = np.pi
    nxobs = 2
    nyobs = 4
    x1d = np.linspace(0, xmax, nxobs, endpoint=False) + xmax / nxobs * 0.5
    y1d = np.linspace(0, ymax, nyobs, endpoint=False) + ymax / nyobs * 0.5
    x2d, y2d = np.meshgrid(x1d, y1d)
    return x2d, y2d

def get_h_full_coverage():
    nobs = NDIM
    h_mat = np.empty((nobs, NDIM))
    xgrid_atm, ygrid_atm = __get_obs_grid_atmos()
    xgrid_ocn, ygrid_ocn = __get_obs_grid_ocean()
    nobs_atm = xgrid_atm.size
    nobs_ocn = xgrid_ocn.size
    assert nobs_atm * 2 + nobs_ocn * 2 == nobs
    for i in range(NDIM):
        state_unit = np.zeros(NDIM)
        state_unit[i] = 1.0
        for j in range(nobs):
            if j < nobs_atm:
                k = j
                elem = "a_psi"
                xgrid = xgrid_atm
                ygrid = ygrid_atm
            elif j < nobs_atm * 2:
                k = j - nobs_atm
                elem = "a_tmp"
                xgrid = xgrid_atm
                ygrid = ygrid_atm
            elif j < nobs_atm * 2 + nobs_ocn:
                k = j - nobs_atm * 2
                elem = "o_psi"
                xgrid = xgrid_ocn
                ygrid = ygrid_ocn
            elif j < nobs_atm * 2 + nobs_ocn * 2:
                k = j - (nobs_atm * 2 + nobs_ocn)
                elem = "o_tmp"
                xgrid = xgrid_ocn
                ygrid = ygrid_ocn
            else:
                raise Exception("__get_h_full_coverage overflow")
            h_mat[j, i] = get_grid_val(state_unit, xgrid.flatten()[k], ygrid.flatten()[k], elem)
    return h_mat

def plot_mat(mat):
    plt.imshow(mat, norm=matplotlib.colors.SymLogNorm(linthresh=0.1),
               cmap="RdBu_r")
    plt.colorbar()
    plt.show()

if __name__ == "__main__":
    mat = get_h_full_coverage()
    plot_mat(mat)

