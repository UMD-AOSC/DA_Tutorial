import sys
import numpy as np
import pickle
import matplotlib
matplotlib.use("TKAgg")
import matplotlib.pyplot as plt

NDIM = 36

def get_grid_val(waves, x, y, is_atm, elem):
    assert waves.__class__ == np.ndarray
    assert waves.shape == (NDIM,)
    assert y.__class__ in [float, np.float32, np.float64]
    assert x.__class__ in [float, np.float32, np.float64]
    assert elem in ["psi", "tmp", "u", "v"]

    f0 = 1.032e-4
    n = 1.5
    na = 10
    no = 8
    R = 287.0
    L = 5000000.0 / np.pi

    gridval = 0.0
    if is_atm:
        types = ["A", "K", "L", "A", "K", "L", "K", "L", "K", "L"]
        hs = [0, 1, 1, 0, 1, 1, 2, 2, 2, 2]
        ps = [1, 1, 1, 2, 2, 2, 1, 1, 2, 2]
        for j in range(na):
            j_all = j + na if elem == "tmp" else j
            if elem == "u":
                if types[j] == "A":
                    gridval += waves[j_all] * np.sqrt(2.0) * ps[j] * np.sin(ps[j] * y)
                elif types[j] == "K":
                    gridval += waves[j_all] * 2.0 * np.cos(hs[j] * n * x) \
                               * ps[j] * (-1) * np.cos(ps[j] * y)
                else:
                    gridval += waves[j_all] * 2.0 * np.sin(hs[j] * n * x) \
                               * ps[j] * (-1) * np.cos(ps[j] * y)
            elif elem == "v":
                if types[j] == "A":
                    gridval += 0.0
                elif types[j] == "K":
                    gridval += waves[j_all] * 2.0 * hs[j] * n * (-1) \
                               * np.sin(hs[j] * n * x) * np.sin(ps[j] * y)
                else:
                    gridval += waves[j_all] * 2.0 * hs[j] * n \
                               * np.cos(hs[j] * n * x) * np.sin(ps[j] * y)
            else:
                if types[j] == "A":
                    gridval += waves[j_all] * np.sqrt(2.0) * np.cos(ps[j] * y)
                elif types[j] == "K":
                    gridval += waves[j_all] * 2.0 * np.cos(hs[j] * n * x) * np.sin(ps[j] * y)
                else:
                    gridval += waves[j_all] * 2.0 * np.sin(hs[j] * n * x) * np.sin(ps[j] * y)
    else:
        hos = [1, 1, 1, 1, 2, 2, 2, 2]
        pos = [1, 2, 3, 4, 1, 2, 3, 4]
        for j in range(no):
            j_all = j + (na * 2 + no) if elem == "tmp" else j + na * 2
            if elem == "u":
                gridval += waves[j_all] * 2.0 * np.sin(0.5 * hos[j] * n * x) \
                           / (-1) * pos[j] * np.cos(pos[j] * y)
            elif elem == "v":
                gridval += waves[j_all] * 2.0 * 0.5 * hos[j] * n \
                           * np.cos(0.5 * hos[j] * n * x) * np.sin(pos[j] * y)
            else:
                gridval += waves[j_all] * 2.0 * np.sin(0.5 * hos[j] * n * x) * np.sin(pos[j] * y)
    if elem == "tmp":
        gridval *= (f0 ** 2 * L ** 2) / R
    elif elem == "psi":
        gridval *= L ** 2 * f0
    else:
        gridval *= L * f0
    return gridval

def __test_get_grid_val():
    state = __model_state_exsample()
    n = 1.5
    x = 1.2 * np.pi / n  # 0.0 <= x <= 2.0 * pi / n
    y = 0.8 * np.pi      # 0.0 <= y <= pi

    # get_grid_val() returns one of four variables
    # {atmosphere|ocean} x {streamfunction|temperature} at point (x, y)
    # unit: [m^2/s] for streamfunction and [K] for temperature
    a_psi = get_grid_val(state, x, y, True, "psi")
    a_tmp = get_grid_val(state, x, y, True, "tmp")
    o_psi = get_grid_val(state, x, y, False, "psi")
    o_tmp = get_grid_val(state, x, y, False, "tmp")
    print(a_psi, a_tmp, o_psi, o_tmp)
    assert np.isclose(-28390877.979826435, a_psi)
    assert np.isclose(-6.915992899895785, a_tmp)
    assert np.isclose(-16019.881464394632, o_psi)
    assert np.isclose(-39.234272164275836, o_tmp)
    __test_difference_u_v(n, state)

def __test_difference_u_v(n, state):
    eps = 1.0e-8
    for is_atm in [True, False]:
        for i in range(1000):
            L = 5000000.0 / np.pi
            pivot_x = np.random.uniform(0, 2.0 * np.pi / n)
            pivot_y = np.random.uniform(0, np.pi)
            ptb_x = pivot_x + eps
            ptb_y = pivot_y + eps
            if ptb_x < 0.0 or ptb_x > 2.0 * np.pi / n or ptb_y < 0.0 or ptb_y > np.pi:
                print("out of range. skip")
                continue
            psi_pivot = get_grid_val(state, pivot_x, pivot_y, is_atm, "psi")
            psi_ptb_x = get_grid_val(state, ptb_x, pivot_y, is_atm, "psi")
            psi_ptb_y = get_grid_val(state, pivot_x, ptb_y, is_atm, "psi")
            u = get_grid_val(state, pivot_x, pivot_y, is_atm, "u")
            v = get_grid_val(state, pivot_x, pivot_y, is_atm, "v")
            try:
                assert np.isclose((psi_ptb_x - psi_pivot) / eps / L, v, atol=1e-6)
                assert np.isclose(- (psi_ptb_y - psi_pivot) / eps / L, u, atol=1e-6)
            except:
                cmp = "atm" if is_atm else "ocn"
                print("assertion error at %s, i = %d" % (cmp, i))

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
    y1d = np.linspace(0, ymax, nyobs, endpoint=False) + ymax / nyobs * 0.25
    x2d, y2d = np.meshgrid(x1d, y1d)
    return x2d, y2d

def __model_state_exsample():
    xini = np.array([
        4.695340259215241E-002,
        2.795833230987369E-002,
        -2.471191763590483E-002,
        -7.877635082773315E-003,
        -4.448292568544942E-003,
        -2.756238610924190E-002,
        -4.224400051368891E-003,
        5.914241112882518E-003,
        -1.779437742222920E-004,
        5.224450720394076E-003,
        4.697982667229096E-002,
        5.149282577209392E-003,
        -1.949084549066326E-002,
        4.224006062949761E-004,
        -1.247786759371923E-002,
        -9.825952138046594E-003,
        -2.610941795170075E-005,
        2.239286581216401E-003,
        -7.891896725509534E-004,
        7.470171905055880E-004,
        -9.315932162526787E-007,
        3.650179005106874E-005,
        1.064122403269511E-006,
        3.937836448211443E-008,
        -2.208288760403859E-007,
        -3.753762121228048E-006,
        -7.105126469908465E-006,
        1.518110190916469E-008,
        -5.773178576933025E-004,
        0.187369278208256,
        1.369868543156558E-003,
        7.023608700166264E-002,
        -4.539810680860224E-004,
        -1.882650440363933E-003,
        -3.900412687995408E-005,
        -1.753655087903711E-007])
    return xini

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
                is_atm = True
                elem = "u"
                xgrid = xgrid_atm
                ygrid = ygrid_atm
            elif j < nobs_atm * 2:
                k = j - nobs_atm
                is_atm = True
                elem = "tmp"
                xgrid = xgrid_atm
                ygrid = ygrid_atm
            elif j < nobs_atm * 2 + nobs_ocn:
                k = j - nobs_atm * 2
                is_atm = False
                elem = "u"
                xgrid = xgrid_ocn
                ygrid = ygrid_ocn
            elif j < nobs_atm * 2 + nobs_ocn * 2:
                k = j - (nobs_atm * 2 + nobs_ocn)
                is_atm = False
                elem = "tmp"
                xgrid = xgrid_ocn
                ygrid = ygrid_ocn
            else:
                raise Exception("__get_h_full_coverage overflow")
            h_mat[j, i] = get_grid_val(state_unit, xgrid.flatten()[k], ygrid.flatten()[k], is_atm, elem)
    return h_mat

def __test_h_matrix_conditional_number():
    h = get_h_full_coverage()
    print(np.linalg.svd(h)[1])
    print(np.linalg.cond(h))
    plot_mat(h)

def plot_mat(mat):
    plt.imshow(mat, cmap="RdBu_r")
    # plt.imshow(mat, norm=matplotlib.colors.SymLogNorm(linthresh=0.1),
    #            cmap="RdBu_r")
    plt.colorbar()
    plt.show()

if __name__ == "__main__":
    __test_h_matrix_conditional_number()

