import numpy as np
import uproot
import matplotlib.pyplot as plt
import matplotlib as mpl
import boost_histogram as bh
from functions import *
import lnl


f = uproot.open(
    "/sps/km3net/users/spenamar/oscillations/chi2_tool/results/v7.10/SelectedEvents_mc_nu.root"
)  # merged DST of events passing selection
# E = f['sel:energy_true'].array().to_numpy() # True energy
E = f["sel:energy_recoJEnergy"].array().to_numpy()
# dir_z = -1.*f['sel:cos_zenith_true'].array().to_numpy() # True cos_zenith
dir_z = -1.0 * f["sel:cos_zenith_recoJGandalf"].array().to_numpy()
nu_type = f["sel:type"].array().to_numpy()
is_cc = f["sel:is_cc"].array().to_numpy()
print(is_cc)
w2 = f["sel:w2"].array().to_numpy()
livetime = f["sel:run_duration"].array().to_numpy()
ngen = f["sel:ngen"].array().to_numpy()


osci_w = compute_evt_weight(nu_type, E, dir_z, is_cc, w2, livetime, ngen)
nu_params = {
    "dm_21": 7.42e-5,
    "dm_31": 2.210e-3,
    # "dm_31": 0,
    "theta_12": 40.0 * np.pi / 180,
    "theta_23": 70.0 * np.pi / 180,
    # "theta_23": 0. * np.pi/180,
    "theta_13": 8.62 * np.pi / 180,
    "dcp": 230 * np.pi / 180,
}
osci_w_bad = compute_evt_weight(
    nu_type, E, dir_z, is_cc, w2, livetime, ngen, nu_params=nu_params
)
print(osci_w)
print(osci_w_bad)


# plothist = lambda h: plt.bar(*h.axes.centers, h.values(), width=h.axes.widths[0], alpha=0.6);

h1 = bh.Histogram(
    bh.axis.Regular(50, 1, 100, transform=bh.axis.transform.log),
)
h2 = bh.Histogram(
    bh.axis.Regular(50, 1, 100, transform=bh.axis.transform.log),
)

h1.fill(E, weight=osci_w)
h2.fill(E, weight=osci_w_bad)
fig, ax = plt.subplots()
plothist(h1, ax)
plothist(h2, ax)
plt.xscale("log")
plt.show()


Ebins = 30
zbins = 20

h1 = bh.Histogram(
    bh.axis.Regular(Ebins, 1, 100, transform=bh.axis.transform.log),
    bh.axis.Regular(zbins, -1, 0),
)
h2 = bh.Histogram(
    bh.axis.Regular(Ebins, 1, 100, transform=bh.axis.transform.log),
    bh.axis.Regular(zbins, -1, 0),
)
h11 = bh.Histogram(
    bh.axis.Regular(Ebins, 1, 100, transform=bh.axis.transform.log),
    bh.axis.Regular(zbins, -1, 0),
)
h22 = bh.Histogram(
    bh.axis.Regular(Ebins, 1, 100, transform=bh.axis.transform.log),
    bh.axis.Regular(zbins, -1, 0),
)
h1.fill(E, -1.0 * dir_z, weight=osci_w)
h2.fill(E, -1.0 * dir_z, weight=osci_w_bad)
h11.fill(E, -1.0 * dir_z, weight=osci_w**2)
h22.fill(E, -1.0 * dir_z, weight=osci_w_bad**2)
s = np.sqrt(h22.values()) / h2.values()
fig, ax = plt.subplots()
plothist2d(h1, ax)
# plothist(h2)
plt.xscale("log")
plt.show()
fig, ax = plt.subplots()
plothist2d(h2, ax)
plt.xscale("log")
plt.show()


print("#" * 20)
print("#" * 20)
print("#" * 20)
print("#" * 20)
vLnL = np.vectorize(lnl.LnL)
print("h1")
print(h1.values())
print(np.shape(h1.values()))
print("h2")
print(h2.values())
print(np.shape(h2.values()))
print("s")
print(s)
print(np.shape(s))
chi2 = vLnL(h1.values(), h2.values(), s)
print("chi2")
print(chi2)
print(np.shape(chi2))
print("#" * 20)
print("#" * 20)
print("#" * 20)
print("#" * 20)

fig, ax = plt.subplots()
h = compute_chi2_map(E, dir_z, osci_w, osci_w_bad)

plothist2d(h, ax)
plt.xscale("log")
plt.show()
