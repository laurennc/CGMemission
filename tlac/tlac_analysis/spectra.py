"""Functions to:
   - obtain spectra
   - compare spectra

A `spectrum` is generally a [x,y,yerr]-list.
"""

import numpy as np
import random
from tlac_plot import plot_spectrum
from physics import lambda0, c


def get_spectrum(dat, weights = None, return_error = False,
                 error_method = 'auto', **kwargs):
    """Returns a binned spectrum ( [x,y,err] ) given the raw frequency (and
    weights).
    
    Keyword Arguments:
    dat           -- List of photons frequencies
    weights       -- List of weights (default: None)
    return_error  -- If True (default: False), will return error as well
    error_method  -- Can be 'auto' (default) or 'poisson'
    **kwargs      -- Will be passed to `tlac_analysis.tlac_plot`,
                     i.e., includes "binning", "smooth", "rng" ...

    Returns:
    x and y values of spectrum (if return_error is True, also the sd on y)
    """
    
    v, y =  plot_spectrum(dat, weights, plot = False, density = False,
                          **kwargs)

    
    normc = (v[1] - v[0]) * np.sum(y)
    if return_error:
        if (weights is None and error_method == 'auto')\
           or error_method == 'poisson':
            ## Poisson error
            err = np.sqrt(y) / normc
        elif error_method == 'bootstrap':
            ## Bootstrap --> To be optimized
            ss_size =  1000. # subsample size
            repetitions = 150
            ypart = []
            for i in range(repetitions):
                # or use np.random.choice with replacement
                vtmp, ytmp = plot_spectrum(random.sample(dat,ss_size),
                                           weights, plot = False,
                                           density = False, **kwargs)
                normctmp = (vtmp[1] - vtmp[0]) * np.sum(ytmp)
                ytmp /= normctmp
                ypart.append(ytmp)
            err = np.std(ypart, axis=0)
        else:
            # "properly" calculated error --> No difference though
            v, wsq = plot_spectrum(dat, weights**2, plot = False,
                                   density = False, **kwargs)
            err = np.sqrt(wsq - 1. / len(dat) * y**2) / normc

    
    y /= normc
    
    if return_error:
        return v, y, err
    else:
        return v, y


    

def compare_spectra(spec_dat, spec_model, check_spectrum_format = False):
    """Compares two spectra.
    
    Keyword Arguments:
    spec_dat              -- [x,y,yerr] of spectrum one (data)
    spec_model            -- [x,y] of spectrum two (model)
    check_spectrum_format -- If true, will check if spectra are compatible
                             (default False)

    Returns:
    log(likelihood)
    """
    v, y = spec_model[:2]
    
    if check_spectrum_format:
        if np.sum( ((v - spec_dat[0])/v)**2 ) > 1e-5:
            raise Exception("Uncomparable formats...")

    dat    = spec_dat[1]
    daterr = spec_dat[2]
    m = daterr == 0
    daterr[m] = np.min(daterr[~m])
    invvar = 1. / (daterr**2)
    lnp = -0.5 * np.sum( (dat - y)**2 * invvar - np.log(invvar))
    
    return lnp

    

    
        
        

    
def spectrum_estimate_error(spec, method = 'continuum', method_args = None):
    """Estimates error from spectrum.
    
    Keyword Arguments:
    spec        -- [x,y]-of (binned) spectrum
    method      -- Method how to estimate the uncertainty (default 'continuum').
                   Available are:
                      + 'contimuum'   -- From spread of contimuum data points.
                                         Excludes `method_args[0 - 1]`


    Returns:
    List of uncertainties for each bin.
    """
    x, y = spec[:2]
    
    if method == 'continuum':
        lims = method_args
        m = ~((x > lims[0]) & (x < lims[1]))
        if np.sum(m) < 25:
            raise Exception("Too few points to calculate uncertainty!")

        return np.repeat(np.std(y[m]), len(x))

    else:
        raise ValueError("Unknown method.")
    
    



def extract_spectrum(spec_in, v_rng, z_rng):
    """Extracts Lya spectrum.
    
    Keyword Arguments:
    spec_in -- Input spectrum (x,y,yerr)
    v_rng   -- Supported velocity range [vmin,vmax] in km/s
    z_rng   -- Chosen z range (zmin,zmax)

    Returns:
    [x,y,yerr] of extracted spectrum
    """
    if len(spec_in) != 3 or len(v_rng) != 2 or len(z_rng) != 2:
        raise ValueError("Unexpected input.")

    l_rng = [ (v_rng[0] * 1e5 / c + 1) * lambda0 * (z_rng[1] + 1),
              (v_rng[1] * 1e5 / c + 1) * lambda0 * (z_rng[0] + 1) ]

    lmax = spec_in[0][ spec_in[1] == np.max(spec_in[1]) ]

    if l_rng[0] > lmax or l_rng[1] < lmax:
        print "[WARNING] Range (%s) does not include lmax (%.2f)" %(str(l_rng),
                                                                    lmax)

    m = (spec_in[0] > l_rng[0]) & (spec_in[0] < l_rng[1])
        
    return [ spec_in[0][m], spec_in[1][m], spec_in[2][m] ]


