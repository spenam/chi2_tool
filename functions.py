import ROOT
import numpy as np
import h5py
import sys,os
import km3flux
import matplotlib.pyplot as plt
import boost_histogram as bh
import lnl


OscProbDir="/pbs/throng/km3net/software/oscprob/master-root_6.18.04"
osc_prob_lib = OscProbDir+"/libOscProb.so"
ROOT.gSystem.Load(osc_prob_lib)
op = ROOT.OscProb


def _infinite(scalar): # From km3services/_tools.py
    """Creates an infinite generator from a scalar.

    Useful when zipping finite length iterators with single scalars.
    """
    while True:
        yield scalar


def _zipequalize(*iterables): # From km3services/_tools.py
    """Creates equal length iterables for a mix of single and n-length arrays"""

    dims = set(map(len, iterables))
    if len(dims) == 1:
        return iterables

    if 1 not in dims or len(dims) > 2:
        raise ValueError("Input arrays dimensions mismatch.")

    out = []
    for it in iterables:
        if len(it) == 1:
            out.append(_infinite(it[0]))
        else:
            out.append(it)

    return out


def _pdgid2flavor(pdgid):
    
    """Converts PDG ID to OscProb flavor"""
    
    if abs(pdgid) == 12:
        return 0
    if abs(pdgid) == 14:
        return 1
    if abs(pdgid) == 16:
        return 2
    raise ValueError("Unsupported PDG ID, please use neutrinos")

    
def _osc_prob(pmns, prem, flav_in, flav_out, energies, cos_zenith):
    params = _zipequalize(flav_in, flav_out, energies, cos_zenith)


    # use a PremModel to make the paths through the earth
    # with the class PremModel
    # chose an angle for the neutrino and fill the paths with cosTheta,
    # e.g. cosTheta = -1 (vertical up-going)
    P = []
    for fl_in, fl_out, E, cosZ in zip(*params):

        if fl_in < 0:
            pmns.SetIsNuBar(True)
        else:
            pmns.SetIsNuBar(False)
        prem.FillPath(cosZ)
        pmns.SetPath(prem.GetNuPath())
        P.append(pmns.Prob(_pdgid2flavor(fl_in), _pdgid2flavor(fl_out), E))
    #print("Here I am")
    #print(P)
    return P

def use_oscprob(flavour_in,flavour_out,energies,cos_zeniths,nu_params=None):
    
    """
    Use oscprob from Zineb's scripts
    """

    #oscprob config
    pmns = op.PMNS_Fast()
    prem = op.PremModel()

    ORCA_DEPTH_KM = 2.4 # km
    prem.SetDetPos(6371.-ORCA_DEPTH_KM) # values from Zineb

    # nufit parameters with SK atmospheric data. values from Nufit2022 in http://www.nu-fit.org/?q=node/256 
    print(nu_params)
    if nu_params is None:
        
        nu_params = {
                "dm_21": 7.41e-5,
                "dm_31": 2.507e-3, # if not inv_hierarchy else dm_21 - 2.465e-3
                "theta_12": 33.41 * np.pi/180, 
                "theta_23": 42.2 * np.pi/180, # if inverted ordering 49.0 * np.pi/180
                "theta_13": 8.58 * np.pi/180, # if inverted ordering 8.57 * np.pi/180
                "dcp": 232 * np.pi/180  # if inverted ordering 276 * np.pi / 180
            }

    pmns.SetDm(2, nu_params["dm_21"])  # set delta_m21 in eV^2
    pmns.SetDm(3, nu_params["dm_31"])  # set delta_m31 in eV^2
    pmns.SetAngle(1, 2, nu_params["theta_12"])  # set Theta12 in radians
    pmns.SetAngle(1, 3, nu_params["theta_13"])  # set Theta13 in radians
    pmns.SetAngle(2, 3, nu_params["theta_23"])  # set Theta23 in radians
    pmns.SetDelta(1, 3, nu_params["dcp"])  # set Delta13 in radians
    print(pmns.GetDm(3))
    print(pmns.GetAngle(2, 3))
    
    
    p = _osc_prob(
    pmns,
    prem,
    flavour_in,
    flavour_out,
    energies,
    cos_zeniths,
    )
    
    return p


def compute_osc_weight(nu_type, energy, dir_z, is_cc,oscillation=True,nu_params=None):

    '''
    Computes oscillation weight based on event properties

    Parameters
    ----------
    nu_type : array
        The flavors (pdgid) of the particles for which to get the weights.
    energy : array
        The energies of the interactions.
    dir_z : array
        The z-dir, cos(theta), of the interaction.
    is_cc : array
        The description of the current of the interaction: 2 is cc, 3 is nc.
    oscillation : bool
        True for consider oscillations, false for not 
    nu_params : dict
        A dict to give specific oscillation parameters. Default are the current best nufit values.

    Returns
    -------
    weight : array
        The oscillation weights for each event.
    '''
    honda = km3flux.flux.Honda()

    #make sure these are all arrays (there were problems with data frames)
    nu_type = np.array(nu_type)
    energy = np.array(energy)
    dir_z = np.array(dir_z)
    is_cc = np.array(is_cc)
 
    
    nu_dict = {
        12: honda.flux(2014, "Frejus", solar="min", averaged="azimuth")["nue"],
        14: honda.flux(2014, "Frejus", solar="min", averaged="azimuth")["numu"],
        -12: honda.flux(2014, "Frejus", solar="min", averaged="azimuth")["anue"],
        -14: honda.flux(2014, "Frejus", solar="min", averaged="azimuth")["anumu"],
    }

    # nufit parameters with SK atmospheric data. values from Nufit2022 in http://www.nu-fit.org/?q=node/256 
    if nu_params is None:

        nu_params = {
                    "dm_21": 7.42e-5,
                    "dm_31": 2.510e-3,
                    "theta_12": 33.45 * np.pi/180,
                    "theta_23": 42.1 * np.pi/180,
                    "theta_13": 8.62 * np.pi/180,
                    "dcp": 230 * np.pi/180,
                    }

    weight = np.zeros(len(nu_type))

    #make sure the correct is_cc convention is used
    #print(is_cc)
    cc_mask = is_cc==1 ######################## SOMETHING WRONG HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #print(cc_mask)
    nc_mask = np.invert(cc_mask)
    #print(cc_mask)
    #print(nu_type[cc_mask])

    for flav in [12, 14]:
        print(flav)

        #get oscillation probs for CC events
        if oscillation:
            #convention: use - dir_z for this
            #flavour_in,flavour_out,...
            #print(flav*np.sign(nu_type[cc_mask]))
            #print(nu_type[cc_mask])

            osc_prob_cc = use_oscprob(flav*np.sign(nu_type[cc_mask]),nu_type[cc_mask],
                                                    energy[cc_mask],-dir_z[cc_mask],nu_params)
            #print(osc_prob_cc)

        else:
            osc_prob_cc = np.zeros(np.count_nonzero(cc_mask))
            flavor_out_is_same_as_flavor_in_mask = flav*np.sign(nu_type[cc_mask]) == nu_type[cc_mask]
            osc_prob_cc[flavor_out_is_same_as_flavor_in_mask] = np.ones(np.count_nonzero(flavor_out_is_same_as_flavor_in_mask))

        #get flux in a loop as it only takes single numbers (as of yet)
        flux  = []
        names = np.vectorize(nu_dict.get)(flav*np.sign(nu_type))
        flux = [nu_dict[flav*np.sign(nu_type[i])](energy[i], dir_z[i]) for i in range(len(nu_type))]
        
        flux = np.asarray(flux)

        #set all the CC weights: osc_prob times flux
        weight[cc_mask] += osc_prob_cc*flux[cc_mask]

        #set all NC weights; total flux, consists of e & mu, no oscillations
        weight[nc_mask] += flux[nc_mask]
        #print(nc_mask)
        #print(cc_mask)

    return weight

def compute_evt_weight(nu_type, energy, dir_z, is_cc, w2, livetime, ngen, oscillation=True,nu_params=None):

    '''
    Computes event weight with oscillations

    Parameters
    ----------
    nu_type : array
        The flavors (pdgid) of the particles for which to get the weights.
    energy : array
        The energies of the interactions.
    dir_z : array
        The z-dir, cos(theta), of the interaction.
    is_cc : array
        The description of the current of the interaction: 2 is cc, 3 is nc.
    oscillation : bool
        True for consider oscillations, false for not 
    nu_params : dict
        A dict to give specific oscillation parameters. Default are the current best nufit values.

    Returns
    -------
    weight : array
        The weights for each event including the oscillations.
    '''
    
    #make sure these are all arrays (there were problems with data frames)
    nu_type = np.array(nu_type)
    energy = np.array(energy)
    dir_z = np.array(dir_z)
    is_cc = np.array(is_cc)
    w2 = np.array(w2)
    livetime = np.array(livetime)
    ngen = np.array(ngen)
    
    oscilation_w = compute_osc_weight(nu_type, energy, dir_z, is_cc,oscillation,nu_params)
    weight = livetime *  w2  / ngen * oscilation_w
    
    return weight

def plothist(h,ax):
    return ax.bar(*h.axes.centers, h.values(), width=h.axes.widths[0], alpha=0.6)

def plothist2d(h,ax):
    return ax.pcolormesh(*h.axes.edges.T, h.values().T)


def compute_chi2_map(E, dir_z, data, mc, Ebins = 30, dir_zbins = 20):

    '''
    Computes the chi2 map from a 2d histogram of data and mc values

    Parameters
    ----------
    E : array
        Energies of the events.
    dir_z : array
        Direction of the events as they come from the DSTs.
    data : array
        Weight of each individual data event.
    mc : array
        Weight of each individual mc event.
    Ebins : int
        Number of bins in the E axis.
    dir_zbins : int
        Number of bins in the dir_z axis.
    

    Returns
    -------
    h : boost_histogram object
        Boost histogram object containing the 2d histogram of dir_z vs E filled with the corresponding Chi2 value for each bin.
    '''
    
    hdata = bh.Histogram(bh.axis.Regular(Ebins,1,100, transform=bh.axis.transform.log),
                  bh.axis.Regular(dir_zbins,-1,0),
                )
    hmc = bh.Histogram(bh.axis.Regular(Ebins,1,100, transform=bh.axis.transform.log),
                      bh.axis.Regular(dir_zbins,-1,0),
                     )
    hdata2 = bh.Histogram(bh.axis.Regular(Ebins,1,100, transform=bh.axis.transform.log),
                      bh.axis.Regular(dir_zbins,-1,0),
                    )
    hmc2 = bh.Histogram(bh.axis.Regular(Ebins,1,100, transform=bh.axis.transform.log),
                      bh.axis.Regular(dir_zbins,-1,0),
                     )
    hdata.fill(E, -1.*dir_z, weight=data)
    hmc.fill(E, -1.*dir_z, weight=mc)
    hdata2.fill(E, -1.*dir_z, weight=data**2)
    hmc2.fill(E, -1.*dir_z, weight=mc**2)
    s = np.sqrt(hmc2.values())/hmc.values()
    
    #vectorize LnL function
    vLnL = np.vectorize(lnl.LnL)
    chi2 = vLnL(hdata.values(),hmc.values(),s)
    
    h = bh.Histogram(bh.axis.Regular(Ebins,1,100, transform=bh.axis.transform.log),
                  bh.axis.Regular(dir_zbins,-1,0),
                )
    h[...] = np.stack(chi2.T, axis=-1)
    
    return h