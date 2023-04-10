import boost_histogram as bh
import uproot 
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from plot_functions import *


def cumulative_1d(hist: bh.Histogram, **kwargs):
    """
    Returns a function of the cumulative distribution of a variable based on the provided boost_histogram object.

    Parameters
    ----------
    hist : boost_histogram.Histogram
        A 1D histogram object from the boost_histogram library.

    Returns
    -------
    cdf : callable
        A function that takes a single numerical argument x and returns the cumulative distribution
        of the variable based on the provided histogram object. If the histogram is 2D, the function
        will calculate the marginal cumulative distribution along the x-axis.

    Raises
    ------
    ValueError
        If the input histogram is not 1D .

    Example
    -------
    >>> import boost_histogram as bh
    >>> import numpy as np
    >>> hist = bh.Histogram(bh.axis.Regular(10, 0, 1))
    >>> hist.fill(np.random.rand(1000))
    >>> cdf_func = cumulative_1d(hist)
    >>> cdf_func(0.5)
    0.494
    """
    if hist.ndim not in [1]:
        raise ValueError("Input histogram must be 1D.")
    hcum = np.cumsum(hist.values())
    hcum = hcum/np.amax(hcum)
    hcum = np.insert(hcum, 0, 0.0)
    cdf = interp1d(hist.axes.edges[0], hcum, kind = "linear")
    return cdf


    

def define_binning_function(fpath: str, is_pid: bool = False, pid: str = "low_track", bin_theta: int = 20, bin_E: int = 30, lowlimE: int = 1, uplimE: int = 100):
    # Define the binning for the histograms in
    # Energy and cos theta
    f = uproot.open(fpath)
    E = f["sel:energy_recoJEnergy"].array().to_numpy()
    if is_pid:
        pid_proba_track = f["sel:pid_proba_track"].array().to_numpy()
        antimu_proba_bkg = f["sel:antimu_proba_bkg"].array().to_numpy()
        classes ={
           "low_track":(E>lowlimE) & (E<uplimE) & (pid_proba_track > 0.85) & (antimu_proba_bkg <= 1e-4), 
           "low_shower":(E>lowlimE) & (E<uplimE) & (pid_proba_track <= 0.85) & (antimu_proba_bkg <= 1e-4), 
           "high_track":(E>lowlimE) & (E<uplimE) & (pid_proba_track > 0.85) & (antimu_proba_bkg > 1e-4) & (antimu_proba_bkg <= 2e-3),
           "high_shower":(E>lowlimE) & (E<uplimE) & (pid_proba_track <= 0.85) & (antimu_proba_bkg > 1e-4) & (antimu_proba_bkg <= 2e-3),
        }
        condition = classes[pid]
        
    else:
        condition = (E>=lowlimE) & (E<=uplimE)
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
        bh.axis.Regular(bin_theta, -1, 0)
    )
    hcostheta.fill(-1.0*dir_z, weight=w)
    h = bh.Histogram(
        bh.axis.Regular(bin_theta, lowlimE, uplimE, transform=bh.axis.transform.log),
        bh.axis.Regular(bin_theta, -1, 0)
    )
    h.fill(E, -1.0*dir_z, weight=w)
    fcum = cumulative_1d(hcostheta)
    hcostheta_flat = bh.Histogram(
        bh.axis.Regular(bin_theta, 0, 1)
    )
    hcostheta_flat.fill(fcum(-1.0*dir_z), weight=w)
    thetac = fcum(-1*dir_z)
    indices = np.argsort(thetac)
    thetac = thetac[indices]
    E = E[indices]
    w = w[indices]
    CDFs = []
    thetaind = []
    for i in range(bin_theta):
        mid_ind = (thetac < hcostheta_flat.axes.edges[0][i+1]) & (thetac > hcostheta_flat.axes.edges[0][i])
        Eis = E[mid_ind]
        wis = w[mid_ind]
        hmid = bh.Histogram(
            bh.axis.Regular(bin_E, lowlimE, uplimE, transform=bh.axis.transform.log)
        )
        hmid.fill(Eis, weight=wis)
        CDFs.append(cumulative_1d(hmid))
        thetaind.append(mid_ind)
    Ecs = []
    for j,item in enumerate(thetaind):
        Ec = CDFs[j](E[item])
        Ecs.append(Ec)
    Ecs = np.concatenate(Ecs,axis=None)

    hnew = bh.Histogram(
        bh.axis.Regular(bin_theta, 0, 1),
        bh.axis.Regular(bin_theta, 0, 1)
    )
    hnew.fill(thetac, Ecs, weight=w)

    return h, hnew


#h, hnew = define_binning_function("/media/spm/a33240d0-5a68-477a-9a84-f5f40d57b19c/Documents/PhD/local_work/dev/chi2_tool_project/SelectedEvents_mc_nu.root")


scale = 7
mii = 0.0001


for item in ["low_track", "low_shower", "high_track", "high_shower"]:
    H, Hnew = define_binning_function("/media/spm/a33240d0-5a68-477a-9a84-f5f40d57b19c/Documents/PhD/local_work/dev/chi2_tool_project/SelectedEvents_mc_nu.root",is_pid=True, pid = item, bin_theta = 20, bin_E = 30)
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