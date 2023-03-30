import boost_histogram as bh
import uproot 
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from plot_functions import *

def define_binning_function(fpath: str, is_pid: bool = False, pid: str = "low_track"):
    # Define the binning for the histograms in
    # Energy and cos theta
    default_theta = 20
    default_E = 30
    lowlim = 1
    uplim = 100
    f = uproot.open(fpath)
    E = f["sel:energy_recoJEnergy"].array().to_numpy()
    if is_pid:
        pid_proba_track = f["sel:pid_proba_track"].array().to_numpy()
        antimu_proba_bkg = f["sel:antimu_proba_bkg"].array().to_numpy()
        classes ={
           "low_track":(E>lowlim) & (E<uplim) & (pid_proba_track > 0.85) & (antimu_proba_bkg <= 1e-4), 
           "low_shower":(E>lowlim) & (E<uplim) & (pid_proba_track <= 0.85) & (antimu_proba_bkg <= 1e-4), 
           "high_track":(E>lowlim) & (E<uplim) & (pid_proba_track > 0.85) & (antimu_proba_bkg > 1e-4) & (antimu_proba_bkg <= 2e-3),
           "high_shower":(E>lowlim) & (E<uplim) & (pid_proba_track <= 0.85) & (antimu_proba_bkg > 1e-4) & (antimu_proba_bkg <= 2e-3),
        }
        condition = classes[pid]
        
    else:
        condition = (E>lowlim) & (E<uplim)
    E = E[condition]
    dir_z = -1.0 * f["sel:cos_zenith_recoJGandalf"].array().to_numpy()[condition]
    #dir_z = dir_z[dir_z>0.0]
    #dir_z = dir_z[dir_z<1.0]
    #E = f["sel:energy_true"].array().to_numpy()
    #dir_z = -1.0 * f["sel:cos_zenith_true"].array().to_numpy()
    w2 = f["sel:w2"].array().to_numpy()/1e17
    w = f["sel:wOsc"].array().to_numpy()[condition]
    ngen = f["sel:ngen"].array().to_numpy()


    hcostheta = bh.Histogram(
        bh.axis.Regular(default_theta, -1, 0)
    )
    hE = bh.Histogram(
        bh.axis.Regular(default_E, lowlim, uplim, transform=bh.axis.transform.log)
    )
    h = bh.Histogram(
        bh.axis.Regular(default_E, lowlim, uplim, transform=bh.axis.transform.log),
        bh.axis.Regular(default_theta, -1, 0),
    )
    hcostheta.fill(-1.0*dir_z, weight=w)
    hE.fill(E, weight=w)
    h.fill(E, -1.0*dir_z, weight=w)
    costhetacum = np.cumsum(hcostheta.values())
    costhetacum = np.insert(costhetacum, 0, 0.0)
    fcum = interp1d(hcostheta.axes.edges[0], costhetacum, kind = "linear")
    Ecum = np.cumsum(hE.values())
    Ecum = np.insert(Ecum, 0, 0.0)
    fcumE = interp1d(hE.axes.edges[0], Ecum, kind = "linear")
    #plt.plot(hE.axes.edges[0], fcumE(hE.axes.edges[0]))
    #plt.show()
    hcostheta_flat = bh.Histogram(
        bh.axis.Regular(20, 0, np.amax(costhetacum)+0.0001)
    )
    hnew = bh.Histogram(
        #bh.axis.Regular(default_E, 1, 100, transform=bh.axis.transform.log)
        bh.axis.Regular(default_E, 0, np.amax(Ecum)+0.0001),
        bh.axis.Regular(default_theta, 0, np.amax(costhetacum)+0.0001),
    )
    #hcostheta_flat.fill(fcum(-1.0*dir_z), weight=w2*E**(-2)/ngen[0])
    hcostheta_flat.fill(fcum(-1.0*dir_z), weight = w)
    #hnew.fill(fcum(-1.0*dir_z), E, weight = w)
    print(Ecum)
    print(fcumE(E))
    hnew.fill(fcumE(E), fcum(-1.0*dir_z), weight = w)
    return h, hnew


H, Hnew = define_binning_function("/home/spm/Documents/PhD/KM3NeT/local_work/chi2_tool/chi2_tool/SelectedEvents_mc_nu.root")
mii = 1
scale = 4
fig, ax = plt.subplots(1,2,figsize = (scale*3,scale*1))
im0 = plothist2d(H,ax=ax[0], vmin = mii)
im1 = plothist2d(Hnew,ax=ax[1], vmin = mii)
ax[0].set_xscale("log")
im0.cmap.set_under('white')
im1.cmap.set_under('white')
fig.colorbar(im0, ax = ax[0])
fig.colorbar(im1, ax = ax[1])
ax[0].set_title("Original binning")
ax[1].set_title("New binning")
ax[0].set_xlabel(r"E [GeV]")
ax[1].set_xlabel(r"CDF(E [GeV])")
ax[0].set_ylabel(r"cos ($\theta$)")
ax[1].set_ylabel(r"CDF(cos ($\theta$))")
plt.tight_layout()
plt.show()
for item in ["low_track", "low_shower", "high_track", "high_shower"]:
    H, Hnew = define_binning_function("/home/spm/Documents/PhD/KM3NeT/local_work/chi2_tool/chi2_tool/SelectedEvents_mc_nu.root",is_pid=True, pid = item)
    fig, ax = plt.subplots(1,2,figsize = (scale*3,scale*1))
    im0 = plothist2d(H,ax=ax[0], vmin = mii)
    im1 = plothist2d(Hnew,ax=ax[1], vmin = mii)
    ax[0].set_xscale("log")
    im0.cmap.set_under('white')
    im1.cmap.set_under('white')
    fig.colorbar(im0, ax = ax[0])
    fig.colorbar(im1, ax = ax[1])
    ax[0].set_title("Original binning")
    ax[1].set_title("New binning")
    ax[0].set_xlabel(r"E [GeV]")
    ax[1].set_xlabel(r"CDF(E [GeV])")
    ax[0].set_ylabel(r"cos ($\theta$)")
    ax[1].set_ylabel(r"CDF(cos ($\theta$))")
    fig.suptitle(item)
    plt.tight_layout()
    plt.show()