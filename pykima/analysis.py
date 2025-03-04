import numpy as np

from .utils import get_planet_mass, get_planet_semimajor_axis

def most_probable_np(results):
    """ 
    Return the value of Np with the highest posterior probability.
    Arguments:
        results: a KimaResults instance. 
    """
    Nplanets = results.posterior_sample[:, results.index_component].astype(int)
    values, counts = np.unique(Nplanets, return_counts=True)
    return values[counts.argmax()]


def passes_threshold_np(results, threshold=150):
    """ 
    Return the value of Np supported by the data, considering a posterior ratio
    threshold (default 150).
    Arguments:
        results: a KimaResults instance. 
        threshold: posterior ratio threshold, default 150.
    """
    Nplanets = results.posterior_sample[:, results.index_component].astype(int)
    values, counts = np.unique(Nplanets, return_counts=True)
    ratios = counts[1:] / counts[:-1]
    if np.any(ratios > threshold):
        i = np.argmax(np.where(ratios > 150)[0]) + 1
        return values[i]
    else:
        return values[0]



def planet_parameters(results, star_mass=1.0, sample=None, printit=True):
    if sample is None:
        sample = results.maximum_likelihood_sample(printit=printit)

    if printit:
        print()
        print('Calculating planet masses with Mstar = %.3f Msun' % star_mass)
    indices = results.indices
    mc = results.max_components

    # planet parameters
    pars = sample[indices['planets']].copy()
    # number of planets in this sample
    nplanets = (pars[:mc] != 0).sum()

    if printit:
        print(20* ' ' + ('%12s' % 'Mp [Mearth]') + ('%12s' % 'Mp [Mjup]'))

    masses = []
    for j in range(int(nplanets)):
        P = pars[j + 0 * mc]
        if P == 0.0:
            continue
        K = pars[j + 1 * mc]
        # phi = pars[j + 2 * mc]
        # t0 = t[0] - (P * phi) / (2. * np.pi)
        ecc = pars[j + 3 * mc]
        # w = pars[j + 4 * mc]
        m_mjup, m_mearth = get_planet_mass(P, K, ecc, star_mass=star_mass)

        if printit:
            s = 20 * ' '
            s += '%12.5f' % m_mearth
            s += '%12.5f' % m_mjup
            print(s)
            masses.append(m_mearth)

    return np.array(masses)



def column_dynamic_ranges(results):
    """ Return the range of each column in the posterior file """
    return results.posterior_sample.ptp(axis=0)

def columns_with_dynamic_range(results):
    """ Return the columns in the posterior file which vary """
    dr = column_dynamic_ranges(results)
    return np.nonzero(dr)[0]