"""A bunch of functions helping the handling of a TlacDB for the tlac_web
interface.
"""
import numpy as np

from bokeh.properties import Instance
from bokeh.models.widgets import (HBox, Slider, TextInput, VBoxForm, Select, PreText,
                                  Paragraph, TableColumn, VBox, RadioButtonGroup)

import tlac_analysis as ta

def hdr_val_to_float(i, v, db = None, raw = True):
    """If raw if True, comes straight from hdr, otherwise just round
    """
    prec = 4
    if not raw:
        return np.round(float(v), prec)
    if(db.choice_keys[i] == 'hydrogen'):
        return np.round(np.log10(float(v[2:])), prec)
    elif(db.choice_keys[i] == 'velocity_params'):
        return np.round(float(v), prec)
    else:
        return np.round(np.log10(float(v)), prec)

def hdr_vals_to_floats(vals, return_tuple = False, **kwargs):
    r = []
    for i, v in enumerate(vals):
        r.append(hdr_val_to_float(i,v, **kwargs))
    if(return_tuple):
        r = tuple(r)
    return r


def db_description_text(db):
    r = "# of grid points:\n"
    for i, ck in enumerate(db.choice_keys):
        r += ck + "\t-->" + str(len(db.choice_lst[i])) + "\n"
    r += "Total = " + str(len(db.db)) + "\n"
    r += "----------------\n"
    r += "# of photons:\n"
    r += "Lya-->" + str(db.db[0].get_nphot('Lya')) + "\n"
    r += "UV -->" + str(db.db[0].get_nphot('UV')) + "\n"
    r += "----------------\n"
    r += "Simulated values:\n"
    r += "sigma_i\t= " + \
         str(db.db[0].header['emission', 'frequency_param']) +\
         "km/s" + "\n"
    r += "tau_d\t= " + \
         str(db.db[0].header['grid', 'dust']) + "\n"
    r += "EW_i\t= " + \
         str(round(db.db[0].get_EWi(),2)) + " AA\n"

    return r


def find_dsid(db, keys, values):
    if '' in values:
        return -1
    
    r = db.lookup(hdr_vals_to_floats(values, return_tuple = True,
                                     raw = False, db = db))
    if r is False:
        return -1
    if r >= 0:
        return r
    else:
        return -1

    
def get_data(db, dsid, freq_type, sigma_i, tau_d, EW_i, nbins, nsmooth):
    if(dsid < 0):
        return [0, [0], [0]]
        
    cds = db.db[dsid]
    if(EW_i == 0):
        EW_i = False
        dat = cds.get_frequencies(freq_type, phot_type = 'Lya')
        weights = ta.tlac_weights(cds, sigma_i, tau_d, EWi = EW_i)
    else:
        w, w_uv = ta.tlac_weights(cds, sigma_i, tau_d, EWi = EW_i)
        dat = cds.get_frequencies(freq_type, phot_type = 'both')
        weights = np.concatenate((w,w_uv))
    if(nbins == 0):
        nbins = 'FD'
    elif nbins > 0:
        nbins = int(nbins)

    if nsmooth == 0:
        x,y = ta.plot_spectrum(dat, weights = weights, plot = False,
                               binning = nbins)
    else:
        x,y = ta.plot_spectrum(dat, weights = weights, plot = False,
                               smooth = {'width' : nsmooth},
                               binning = nbins)

    return [np.sum(weights)] + [x,y]


def get_sigma_i_range(db, frng = [0.1, 1.5]):
    """Returns sigma_i range corresponding to
    frng * sigma_i
    """
    cds = db.db[0]
    si = float(cds.header['emission','frequency_param'])
    return [ i * si for i in frng ]

    

def create_input_instance(db, i):
    """Returns instance of `Select` or `Slider` depending on value of `i`
    """
    o = db.choice_type[i]
    if(o == 0):
        return Instance(Select)
    elif(o == 1):
        return Instance(Slider)
    else:
        raise ValueError("i must be 0 or 1!")

    
    

    
