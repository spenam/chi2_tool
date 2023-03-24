import numpy as np
import uproot
import matplotlib.pyplot as plt
import matplotlib as mpl
import boost_histogram as bh
from functions import *
import lnl


f = uproot.open("/sps/km3net/users/alflazo/dstProd/v7.1_v7.2_jsh/mcv7.2.gsg_muon-CC_1-100GeV.km3sim.jorcarec.jsh.aanet.dst_merged.root") # nu mu CC
#f = uproot.open("/sps/km3net/users/alflazo/dstProd/v7.1_v7.2_jsh/mcv7.2.gsg_elec-CC_1-100GeV.km3sim.jorcarec.jsh.aanet.dst_merged.root") # nu e CC
#f = uproot.open("/sps/km3net/users/alflazo/dstProd/v7.1_v7.2_jsh/mcv7.2.gsg_muon-NC_1-100GeV.km3sim.jorcarec.jsh.aanet.dst_merged.root") # nu mu NC
#f = uproot.open("/sps/km3net/users/spenamar/oscillations/chi2_tool/results/v7.10/SelectedEvents_mc_nu.root") # merged DST of events passing selection
#E = f['E:Evt/mc_trks/mc_trks.E'].array()[:,0].to_numpy()
#dir_z = f['E:Evt/mc_trks/mc_trks.dir.z'].array()[:,0].to_numpy()
E = f['E:Evt/trks/trks.E'].array()[:,0].to_numpy()
dir_z = f['E:Evt/trks/trks.dir.z'].array()[:,0].to_numpy()
nu_type = f['E:Evt/mc_trks/mc_trks.type'].array()[:,0].to_numpy()
is_cc = f['T:sum_mc_nu/cc'].array().to_numpy()
print(is_cc)
w2 = f['E:Evt/w'].array()[:,1].to_numpy()
w_noOsc = f['T:sum_mc_evt/weight_noOsc'].array().to_numpy()
w_Osc = f['T:sum_mc_evt/weight'].array().to_numpy()
livetime = f['T:sum_mc_evt/livetime_DAQ'].array().to_numpy()
ngen = f['T:sum_mc_evt/n_gen'].array().to_numpy()


osci_w = compute_evt_weight(nu_type,E,dir_z,is_cc,w2,livetime,ngen)
nu_params = {
                        "dm_21": 7.42e-5,
                        "dm_31": 2.210e-3,
                        #"dm_31": 0,
                        "theta_12": 40. * np.pi/180,
                        "theta_23": 70. * np.pi/180,
                        #"theta_23": 0. * np.pi/180,
                        "theta_13": 8.62 * np.pi/180,
                        "dcp": 230 * np.pi/180,
                        }
osci_w_bad = compute_evt_weight(nu_type,E,dir_z,is_cc,w2,livetime,ngen, nu_params=nu_params)
print(osci_w)
print(osci_w_bad)


#plothist = lambda h: plt.bar(*h.axes.centers, h.values(), width=h.axes.widths[0], alpha=0.6);

h1 = bh.Histogram(bh.axis.Regular(50,1,100, transform=bh.axis.transform.log),
                )
h2 = bh.Histogram(bh.axis.Regular(50,1,100, transform=bh.axis.transform.log),
                 )
h3 = bh.Histogram(bh.axis.Regular(50,1,100, transform=bh.axis.transform.log),
                 )

h1.fill(E, weight=osci_w)
h2.fill(E, weight=osci_w_bad)
fig, ax = plt.subplots()
plothist(h1,ax)
plothist(h2,ax)
plt.xscale("log")
plt.show()

h1 = bh.Histogram(bh.axis.Regular(50,1,100, transform=bh.axis.transform.log),
                )
h2 = bh.Histogram(bh.axis.Regular(50,1,100, transform=bh.axis.transform.log),
                 )
h3 = bh.Histogram(bh.axis.Regular(50,1,100, transform=bh.axis.transform.log),
                 )
h1.fill(E, weight=osci_w)
h2.fill(E, weight=w_noOsc)
h3.fill(E, weight=w_Osc)
fig, ax = plt.subplots()
plothist(h1,ax)
#plothist(h2)
plothist(h3, ax)
plt.xscale("log")
plt.show()

h1 = bh.Histogram(bh.axis.Regular(30,1,100, transform=bh.axis.transform.log),
                )
h2 = bh.Histogram(bh.axis.Regular(30,1,100, transform=bh.axis.transform.log),
                 )
h11 = bh.Histogram(bh.axis.Regular(30,1,100, transform=bh.axis.transform.log),
                )
h22 = bh.Histogram(bh.axis.Regular(30,1,100, transform=bh.axis.transform.log),
                 )
h1.fill(E, weight=osci_w)
h2.fill(E, weight=w_noOsc)
h11.fill(E, weight=osci_w**2)
h22.fill(E, weight=w_noOsc**2)
fig, ax = plt.subplots()
plothist(h1, ax)
plothist(h2, ax)
plt.xscale("log")
plt.show()
fig, ax = plt.subplots()
plothist(h11, ax)
plothist(h22, ax)
plt.xscale("log")
plt.show()

#def plothist2d(h):
#        return plt.pcolormesh(*h.axes.edges.T, h.values().T)



Ebins = 30
zbins = 20

h1 = bh.Histogram(bh.axis.Regular(Ebins,1,100, transform=bh.axis.transform.log),
                  bh.axis.Regular(zbins,-1,0),
                )
h2 = bh.Histogram(bh.axis.Regular(Ebins,1,100, transform=bh.axis.transform.log),
                  bh.axis.Regular(zbins,-1,0),
                 )
h11 = bh.Histogram(bh.axis.Regular(Ebins,1,100, transform=bh.axis.transform.log),
                  bh.axis.Regular(zbins,-1,0),
                )
h22 = bh.Histogram(bh.axis.Regular(Ebins,1,100, transform=bh.axis.transform.log),
                  bh.axis.Regular(zbins,-1,0),
                 )
h1.fill(E, -1.*dir_z, weight=osci_w)
h2.fill(E, -1.*dir_z, weight=osci_w_bad)
h11.fill(E, -1.*dir_z, weight=osci_w**2)
h22.fill(E, -1.*dir_z, weight=osci_w_bad**2)
s = np.sqrt(h22.values())/h2.values()
fig, ax = plt.subplots()
plothist2d(h1,ax)
#plothist(h2)
plt.xscale("log")
plt.show()
fig, ax = plt.subplots()
plothist2d(h2,ax)
plt.xscale("log")
plt.show()


print("#"*20)
print("#"*20)
print("#"*20)
print("#"*20)
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
chi2 = vLnL(h1.values(),h2.values(),s)
print("chi2")
print(chi2)
print(np.shape(chi2))
print("#"*20)
print("#"*20)
print("#"*20)
print("#"*20)

fig, ax = plt.subplots()
h = compute_chi2_map(E, dir_z, osci_w, osci_w_bad)

plothist2d(h,ax)
plt.xscale("log")
plt.show()


