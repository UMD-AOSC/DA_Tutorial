#!/usr/bin/env python

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import subprocess

N = 36
# N = 228   # atm 6x6 ocn 6x6
# N = 318   # atm 6x6 ocn 9x9
# N = 414   # atm 9x9 ocn 6x6
DT = 10.0  # write interval in [timeunit]
NT = 100  # number of write. filesize = 8 * (N ** 2 + N) * NT [bytes]
ONEDAY = 8.64  # [timeunit/day] a46p51
FNAME = "dat/evol_field_tlm_2.00.dat"
GINELLI = True
NT_ABORT = 0
T_VERIF_LIST = [0]

def main():
    np.random.seed(10 ** 8 + 7)
    np.set_printoptions(formatter={'float': '{: 0.4f}'.format}, threshold=np.inf, linewidth=400)
    mkdirs()
    if GINELLI:
        trajs, gs, ms, rs = integ_forward_ginelli()
        cs = integ_backward_ginelli(rs)
        vs = obtain_clvs_ginelli(cs, gs)
        # print_verif(trajs, ms, gs, rs, cs, vs)
    else:
        trajs, gs, ms = integ_forward()
        fs = integ_backward(ms)
        vs = calc_clv(gs, fs)
    trajs, gs, ms, vs = trim_spinup(trajs, gs, ms, vs)
    test_growth_rate(ms, vs)
    # test_growth_rate_long(ms, vs, 100, [0, 1, 4, 9, 14, 19, 24, 29, 35])

def trim_spinup(trajs, gs, ms, vs):
    st = NT_ABORT
    ed = NT - NT_ABORT
    return trajs[st:ed, :], gs[st:ed, :, :], ms[st:ed, :, :], vs[st:ed, :, :]

def integ_forward():
    all_traj = np.empty((NT, N))
    all_blv = np.empty((NT, N, N))
    all_tlm = np.empty((NT, N, N))
    blv = np.random.normal(0.0, 1.0, (N, N))
    blv, ble = orth_norm_vectors(blv)
    for i in range(NT):
        traj, tlm = read_file_part(FNAME, i)
        all_traj[i, :] = traj
        all_blv[i, :, :] = blv
        all_tlm[i, :, :] = tlm
        blv = np.dot(tlm, blv)
        blv, ble = orth_norm_vectors(blv)
    return all_traj, all_blv, all_tlm

def integ_backward(all_tlm):
    all_flv = np.empty((NT, N, N))
    flv = np.random.normal(0.0, 1.0, (N, N))
    flv, fle = orth_norm_vectors(flv)
    for i in reversed(range(NT)):
        traj, tlm = read_file_part(FNAME, i)
        flv = np.dot(tlm.T, flv)
        flv, fle = orth_norm_vectors(flv)
        all_flv[i, :, :] = flv
    return all_flv

def calc_clv(all_blv, all_flv):
    all_clv = np.empty((NT, N, N))
    for i in range(NT):
        for k in range(N):
            all_clv[i, :, k] = vector_common(all_blv[i, :, :k+1], all_flv[i, :, k:], k)
    return all_clv

def orth_norm_vectors(lv):
    lv_oath, r = np.linalg.qr(lv)
    le = np.zeros(N)
    eigvals = np.abs(np.diag(r))
    le = np.log(eigvals)
    return lv_orth, le

def vector_common(blv, flv, k):
    ab = np.empty((N, N+1))
    ab[:, :k+1] = blv[:, :]
    ab[:, k+1:] = - flv[:, :]
    coefs = nullspace(ab)
    clv = np.dot(ab[:, :k+1], coefs[:k+1, np.newaxis])[:,0]
    clv /= np.linalg.norm(clv)
    return clv

def nullspace(a):
    # refer to 658Ep19
    u, s, vh = np.linalg.svd(a)
    return vh.T[:,-1]

def integ_forward_ginelli():
    # trajs[i, :] and gs[i, :, :] are at time i
    # ms[i, :, :] and rs[i, :, :] are times i -> i + 1
    # g = np.random.normal(0.0, 1.0, (N, N))
    g = np.identity(N)
    g, r = np.linalg.qr(g)
    trajs = np.empty((NT, N))
    gs = np.empty((NT, N, N))
    ms = np.empty((NT, N, N))
    rs = np.empty((NT, N, N))
    for i in range(NT):
        traj, m = read_file_part(FNAME, i)
        trajs[i, :] = traj
        gs[i, :, :] = g
        g = np.dot(m, g)
        g, r = np.linalg.qr(g)
        ms[i, :, :] = m
        rs[i, :, :] = r
    return trajs, gs, ms, rs

def integ_backward_ginelli(rs):
    # rs[i, :, :] is times i -> i + 1
    # cs[i, :, :] is at time i
    cs = np.empty((NT, N, N))
    # c = np.random.normal(0.0, 1.0, (N, N))
    c = np.identity(N)
    dummy, c = np.linalg.qr(c)
    c = normalize_column(c)
    for i in reversed(range(NT)):
        r = rs[i, :, :]
        ri = np.linalg.inv(r)
        c = np.dot(ri, c)
        c = normalize_column(c)
        cs[i, :, :] = c
    return cs

def obtain_clvs_ginelli(cs, gs):
    # cs[i, :, :], gs[i, :, :], and vs[i, :, :] are at time i
    vs = np.empty((NT, N, N))
    for i in range(NT):
        c = cs[i, :, :]
        g = gs[i, :, :]
        v = np.dot(g, c)
        vs[i, :, :] = v
    return vs

def orth_norm_vectors(lv):
    # lv     <- np.array[N,N] : Lyapunov vectors (column)
    # return -> np.array[N,N] : orthonormalized LVs in descending order
    # return -> np.array[N]   : ordered Lyapunov Exponents
    q, r = np.linalg.qr(lv)
    le = np.zeros(N)
    eigvals = np.abs(np.diag(r))
    lv_orth = q
    le = np.log(eigvals)
    return lv_orth, le

def normalize_column(m):
    for i in range(N):
        c = m[:, i]
        lc = np.dot(c, c) ** 0.5
        m[:, i] /= lc
    return m

def test_growth_rate(ms, vs):
    nt = ms.shape[0]
    rate = np.zeros(N)
    for i in range(nt):
        v = vs[i, :, :]
        m = ms[i, :, :]
        vn = np.dot(m, v)
        for j in range(N):
            r = np.log(np.linalg.norm(vn[:, j]) / np.linalg.norm(v[:, j]))
            rate[j] += r / (DT / ONEDAY) / nt
    print("Growth rate:")
    print(rate)

def test_growth_rate_long(ms, vs, ntg, ilist):
    nt = ms.shape[0]
    lengths = np.zeros((nt - ntg, ntg + 1, len(ilist)))
    for k, i in enumerate(ilist):
        for it in range(nt - ntg):
            vi = vs[it, :, i]
            lengths[it, 0, k] = np.linalg.norm(vi)
            for jt in range(ntg):
                vi = np.dot(ms[it + jt, :, :], vi[:])
                lengths[it, jt + 1, k] = np.linalg.norm(vi)
    length_log_mean = np.mean(np.log(lengths), axis=0)
    length_mean_log = np.log(np.mean(lengths, axis=0))
    growth_log_mean = np.copy(length_log_mean)
    growth_mean_log = np.copy(length_mean_log)
    for jt in reversed(range(ntg)):
        growth_log_mean[jt + 1, :] -= growth_log_mean[jt, :]
        growth_mean_log[jt + 1, :] -= growth_mean_log[jt, :]
    plot_growth_rate(growth_log_mean, ntg, ilist, "lognorm")
    plot_growth_rate(growth_mean_log, ntg, ilist, "L2_norm")

def plot_growth_rate(growth_mean, ntg, ilist, suff):
    # growth_mean.shape = (ntg + 1, len(ilist))
    x = (np.arange(float(ntg)) + 0.5) * (DT / ONEDAY)
    for k, i in enumerate(ilist):
        plt.plot(x, growth_mean[1:, k], label="CLV %02d" % (i + 1))
    plt.xlabel("days")
    plt.legend()
    plt.savefig("img/fig_8_%s.png" % suff)
    plt.close()

def print_verif(trajs, ms, gs, rs, cs, vs):
    # trajs[i, :], gs[i, :, :], cs[i, :, :], and vs[i, :, :] are at time i
    # ms[i, :, :] and rs[i, :, :] are times i -> i + 1
    for i in T_VERIF_LIST:
        print(i)
        if i < NT:
            print("Traj:")
            print(trajs[i, :])
            print("Q:")
            print(gs[i, :, :])
            print("C:")
            print(cs[i, :, :])
            print("CLVs:")
            print(vs[i, :, :])
        elif i == NT:
            print("Traj: None")
            g = np.dot(ms[i - 1, :, :], gs[i - 1, :, :])
            g, r = np.linalg.qr(g)
            print("Q:")
            print(g)
            assert np.allclose(r, rs[i - 1, :, :])
            print("C: Identity")
            print("CLVs: same as Q")
        if i > 0:
            print("R:")
            print(rs[i - 1, :, :])

def mkdirs():
    subprocess.run("rm -rf img", check=True, shell=True)
    subprocess.run("mkdir -p img", check=True, shell=True)

def read_file_part(fname, it):
    dbyte = it * (N ** 2 + N) * 8
    x = np.array(np.memmap(fname, offset=dbyte, dtype=np.float64, mode='r', shape=(N ** 2 + N)))
    assert x.shape == (N ** 2 + N, )
    traj = x[:N]
    tlm = x[N:].reshape((N, N))
    return traj, tlm

if __name__ == "__main__":
    main()