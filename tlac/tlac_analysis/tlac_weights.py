import numpy as np
from scipy.stats import norm 
from physics import c as c_cgs
from physics import thermal_velocity, x_to_v

from TlacDatSet import TlacDatSet

def tlac_weights(cds, sigmap, taup, EWi = None, T = 'auto', cHI = 'auto'):
    """
    Returns weights of photons to mimic another intrinsic spectrum and/or dust.
    The intrinsic spectrum it will model is a gaussian with standard deviation
    of sigmap / thermal velocity.
    If EWi is given, also the continuum will be modelled.
    
    Keyword Arguments:
    cds    -- TlacDatSet instance of data
    sigmap -- intrinsic sigma to model (in km/s)
    taup   -- dust optical depth to model. If None, no dust is added.
    EWi    -- Intrinsic equivalent width in Angstrom.
    T/cHI  -- can be 'auto' (default) for automatic reading from header
              (only works for some filling modes). If float is given, will
              overwrite this.

    Returns:
    For EWi = None (default): List of weight for Lya photons.
    Otherwise: [weights(Lya), weights(UV) ]

    Details:
    The weight of a photon is given by
    \begin{equation}
    w = \frac{f_p(x_i)}{f_r(x_i)}
    \end{equation}
    where $f_p$ is the pdf (the intrinsic spectrum) one wants to model and
    $f_r$ is the one you actually ran.
    For dust the weight of each photon is given by
    \begin{equation}
    w = \e^{-\tau_d^{(d)}}
    \end{equation}
    with
    \begin{equation}
    \tau_d^{(d)} = c_d^{(d)}\sigma_d = c_{HI}^{(d)}
    \frac{c_d^{(c)}}{c_{HI}^{(c)}}\sigma_d = c_{HI}^{(c)}
    \frac{\tau_d^{c}}{c_{HI}^{(c)}}.
    \end{equation}
    Here, we used ${}^{(d)}$ for symbols that are from the data and
    ${}^{(c)}$ for quantities from the config file.

    See arXiv:1506.03836 for details.
    """
    # Some backward compatibility
    if EWi is 0 or EWi is False:
        EWi = None
    
    if(cds.header['emission', 'frequency_mode' ] != 2):
         # Comment: More intrinsic emission PDFs can be added.
        raise Exception("Needs Gaussian emission (frequency_mode = 2) to weight"
                        " photons.")
    
    ## Intrinsic weights
    if EWi is not None:
        weights_lya, weights_uv = _intrinsic_weights_combined(cds, sigmap, EWi,
                                                              T = T)
    else:
        weights_lya = _intrinsic_weights_lya(cds, sigmap, T = T)
        
    ## Dust weighting
    if taup is not None:
        w_dust = _dust_weights(cds, taup, EWi is not None, cHI = cHI)
        
    assert np.sum(weights_lya) > 0, "Zero weights assigned!"


    if EWi is not None:
        return weights_lya * w_dust[0], weights_uv * w_dust[1]
    else:
        return weights_lya * w_dust
    
    
        
        
def compute_fesc(cds, sigmap, taud, EWi, T = 'auto', cHI = 'auto'):
    """Computes theoretical escape fraction of Lya photons through
    weighting.

    NOT YET IMPLEMENTED ONLY FOR A TLAC RUN ONLY WITH LYA PHOTONS
    (BUT OF COURSE COMPUTES f_esc_lya). 
    
    Keyword Arguments:
    cds    -- TlacDatSet to compute
    sigmap -- sigma_i to use
    taud   -- dust optical depth to use
    T, cHI -- Temperature, hydrogen column density to use
              (default 'auto' = try to read from config)
    """
    f_post, f_run = _intrinsic_weights_combined(cds, sigmap, EWi,
                                                T = T,
                                                return_pdf_values = True)

    if T == 'auto':
        T = cds.header['grid', 'temperature' ]

    vi  = x_to_v(np.concatenate((cds['Lya','x_i'], cds['UV','x_i'])), T) / 1e5
    f_lya = norm.pdf(vi, 0, sigmap)

    N = np.sum(f_lya / f_run) # normalization
    w_dust = np.concatenate(_dust_weights(cds, taud, True, cHI = cHI))
    fesc = np.sum( f_lya / f_run * w_dust) / N

    return fesc
    
    

    
        
        

def _intrinsic_weights_lya(cds, sigmap, T = 'auto'):
    """
    Internal function called from `tlac_weights` to mimic intrinsic spectrum
    of the Lyman alpha photons.
    
    Keyword Arguments:
    cds    -- TlacDatSet instance of data
    sigmap -- intrinsic sigma to model (in km/s)
    T      -- can be 'auto' (default) for automatic reading from header
              (only works for global filling mode). If float is given, will
              overwrite this.

    Returns:
    weights(Lya)
    """

    sigmar = float(cds.header['emission', 'frequency_param' ])

    if T == 'auto':
        T = cds.header['grid', 'temperature' ]
    #b = 0.128444 * np.sqrt(T) # thermal motion in km/s
    #b = cds.get_v_thermal() / 1e3
    b = thermal_velocity(T) / 1e5

    xi = cds['Lya','x_i']
    weights_lya = sigmar / sigmap * np.exp(-0.5 * (b * xi)**2 * (sigmap**(-2) -\
                                                                 sigmar**(-2)) )

    return weights_lya
    



def _intrinsic_weights_combined(cds, sigmap, EWi, T = 'auto',
                                return_pdf_values = False):
    """
    Internal function called from `tlac_weights` to mimic intrinsic spectrum
    including the continuum.
    
    Keyword Arguments:
    cds    -- TlacDatSet instance of data
    sigmap -- intrinsic sigma to model (in km/s)
    EWi    -- Intrinsic equivalent width in Angstrom.
    T      -- can be 'auto' (default) for automatic reading from header
              (only works for global filling mode). If float is given, will
              overwrite this.
    return_pdf_values -- If True (default: False) returns instead two lists
                         with f_post and f_run values

    Returns:
    [weights(Lya), weights(UV) ]
    """
    if(cds.header['emission', 'uv_frequency_mode'] != 3):
        raise Exception("uv_frequency_mode = 3 expected to calculate "
                        "continuum weights!")

    deltav = 2. * cds.header['emission', 'uv_frequency_param' ]
    nlya = cds.get_nphot('Lya')
    nuv = cds.get_nphot('UV')
    eta = nuv / float(nlya)
    if T == 'auto':
        v_th = cds.get_v_thermal() / 1e3 # in km/s
    else:
        v_th = thermal_velocity(T) / 1e5
    c = c_cgs / 1e5 # in km/s
    xi = np.concatenate((cds['Lya','x_i'], cds['UV','x_i']))
    vi = c / (1. - c / (v_th * xi))
    exparg = -0.5 * vi**2
    
    sigmar = float(cds.header['emission', 'frequency_param' ])

    etap = deltav * 1215.67 / (c * EWi)

    f_post   = ((eta + 1) * (etap / deltav + np.exp(exparg / sigmap**2) /
                             (np.sqrt(2 * np.pi) * sigmap)))
    f_run    = ((etap + 1) * (eta / deltav + np.exp(exparg / sigmar**2) /
                           (np.sqrt(2 * np.pi) * sigmar)))
    w = f_post / f_run

    if return_pdf_values:
        return f_post, f_run
    else:
        return [ w[:nlya], w[nlya:] ]
    
    
    

def _dust_weights(cds, taup, uv_weights = False,
                  cHI = 'auto'):
    """
    Returns weighting due to dust.
    
    Keyword Arguments:
    cds  -- 
    taup -- Dust optical depth
    cHI  -- Can give a default column density (default 'auto' = will try to
            read from config)
    """
    if cHI == 'auto':
        hydrogen = cds.header['grid','hydrogen']
        if(hydrogen[0] != 'c'):
            raise Exception("Needs hydrogen given as column density to compute "
                            "weights!")
        cHIc = float(hydrogen[2:])
    else:
        cHIc = cHI
    cHId = cds['Lya','c_HI_total']

    weights_lya = np.exp(-cHId / cHIc * taup)

    if uv_weights:
        cHId_uv = cds['UV','c_HI_total']
        weights_uv = np.exp(-cHId_uv / cHIc * taup)
        return weights_lya, weights_uv
            
    else:
        return weights_lya
        
