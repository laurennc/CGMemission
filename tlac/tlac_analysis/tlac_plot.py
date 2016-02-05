from weighted_kde import Covariator, gaussian_kde
from tlac_weights import tlac_weights
from TlacDatSet import TlacDatSet

from matplotlib.widgets import Slider

import matplotlib
import matplotlib.pyplot as plt

import numpy as np
from scipy.signal import gaussian

import brewer2mpl

from physics import pc

def plot_init(usetex = False):
    """Initializes plotting. Changes style etc.
    """
    plt.clf()
    
    plot_set_default_palette()
    
    matplotlib.rcParams['lines.linewidth'] = 1.2

    if(usetex):
        matplotlib.rc('text', usetex=True) 
        matplotlib.rc('font', family='serif')

    #matplotlib.rcParams['lines.linewidth'] = 2

    plt.grid(True)


def plot_get_colors(mode = ['brewer','qualitative','Set1'],
                    ncolors = 8):
    """Gets default colors
    """
    if(mode[0] == 'brewer'):
        bmap = brewer2mpl.get_map(mode[2], mode[1], ncolors)
        colors = bmap.mpl_colors
        colors.pop(5) # Remove yellow
    else:
        raise ValueError("Wrong value for 'mode'!")

    return colors

def plot_set_default_palette(**kwargs):
    """Sets default pallette accordingly.
    """
    matplotlib.rcParams['axes.color_cycle'] = plot_get_colors(**kwargs)


def plot_spectrum(dat, weights = None, rng = None, binning = 'FD', kde = False,
                  plot = True, uv_dat = None, uv_weights = None, smooth = None,
                  density = True, scale_frequencies = 1.,
                  **kwargs):
    """
    Keyword Arguments:
    dat        -- array of photon frequencies
    weights    -- None (default) for no weighting, otherwise array of weights
    rng        -- None (default) for min to max of data, otherwise (min,max)
    binning    -- how many bins to use. Can be
                   * 'FD' for Friedmann-Diaconis rule (default)
                   * integer for fixed number of bins
                   * float > 0 for average number of photons per bin
                   * float < 0 for binwidth
    kde        -- if true, kde will be overplotted. (default: False)
    plot       -- If false, instead of plotting [xlst, ylst] are returned.
    uv_dat     -- If UV photon frequencies are given, also continuum is plotted.
    uv_weights -- uv weights can be given
    smooth     -- None (default) for no smoothing. Otherwise dictionary which
                  will be passed to `ta.smooth_spectrum` (not x,y). Note that
                  `binning` will be overwritten if not given explicitly.
    density    -- `density` argument passed to numpy.histogram (default: True)
    scale_frequencies -- Multiplies frequencies (x-axis) by this number
    **kwargs   -- will be passed to matplotlib.plot
    """
    if(weights is None):
        weights = np.ones(len(dat))

    if(rng is None):
        rng = [np.min(dat[weights > 0]),
               np.max(dat[weights > 0]) ]

    #n = np.sum(weights)
    n = len(dat)
    if(binning == 'FD'):
        if smooth is None:
            p25, p75 = np.percentile(dat, [25, 75])
            iqr = p75 - p25
            h = iqr / n**(1.0/3.0)
            nbins = np.round( (rng[1] - rng[0]) / h )
        else:
            # choose bins when smoothing
            nbins = np.round( (rng[1] - rng[0]) / smooth['width'] ) * 25
    elif(isinstance(binning, int)):
        nbins = binning
    elif(isinstance(binning, float)):
        if(binning > 0):
            nbins = np.round( n / binning )
        else:
            nbins = np.round((-1.0) * (rng[1] - rng[0]) / binning)
    else:
        raise ValueError("Unknown value for 'binning'. Can be 'FD', "
                         "integer or float.")

    hist = np.histogram(dat, bins = nbins, range = rng, weights = weights,
                        density = density)
    
    xlst = 0.5 * (hist[1][:-1] + hist[1][1:])
    ylst = hist[0]

    if smooth is not None:
        xlst, ylst = smooth_spectrum(xlst, ylst, **smooth)
        
    if(plot):
        # save next color in cycle so can be re-used for KDE
        if('color' in kwargs):
            color = kwargs['color']
        else:
            cc = plt.gca()._get_lines.color_cycle
            color = cc.next()
            kwargs['color'] = color

        if density:
            yscale = scale_frequencies
        else:
            yscale = 1.0
        plt.plot(xlst / scale_frequencies, ylst * yscale, **kwargs)
    else:
        return [xlst, ylst]

    if(kde):
        plot_kde(dat, weights, rng = rng, style = 'shady', color = color)



def plot_kde(dat, weights = None, rng = None, resolution = 10, style = 'normal',
             bw_method = 'scott', plot = True, **kwargs):
    """
    Calculates & plots kernel density estimator
    
    Keyword Arguments:
    dat                 -- Data array
    weights             -- Weights of data (default None)
    rng                 -- Range to plot (default None = from min to max)
    resolution          -- Plot points per unit of rng (default 10)
    bw_method           -- Bandwith selection method. Can be 'scott' (default),
                           'silverman' or a constant
    style               -- Can be 'normal' for normal lines or 'shady' for
                           thicker, more transparent lines
    plot                -- If true (default) kde will be plotted. Otherwise
                           [x,y] is returned.
    **kwargs            -- will be passed to matplotlib.plot
    """
    if(weights is None):
        weights = np.ones(len(dat))

    if(rng is None):
        rng = [np.min(dat), np.max(dat) ]
        
    if(style == 'shady'):
        if(not 'alpha' in kwargs):
            kwargs['alpha'] = 0.5
        if(not 'linewidth' in kwargs):
            kwargs['linewidth'] = 1.8

    cov = Covariator([dat], weights)
    inv_cov, normf = cov(bw_method)
    kde = gaussian_kde([dat], weights, inv_cov, normf)

    x = np.linspace(rng[0],rng[1], np.round( resolution * (rng[1] - rng[0])))
    y = kde.evaluate([x])
    
    if(plot):
        plt.plot(x,y,**kwargs)
    else:
        return [x, y]


    

def iplot_spectrum_ds(ds, init = True,
                      sigmarng = [0.05, 1.5], taurng = [0,4], **kwargs):
    """
    Plots spectrum from TlacDatSet interactively using matplotlib.
    NOT FINISHED >> see tlac_web instead
    
    Keyword Arguments:
    ds         -- TlacDatSet instance to plot
    init       --- If true (default) plot will be initialized first
    sigmarng   -- Range for intrinsic spectrum width. Measured in
                  units of the actual sigma_i (default [0,1.5])
    taurng     -- Range for dust optical depth tau_d (default [0, 4])
    **kwargs   -- Will be passes to `tlac_analysis.plot_spectrum`
    """
    if(not isinstance(ds, TlacDatSet)):
        raise ValueError("ds has to be a TlacDatSet instance")

    dat = ds['Lya','x']
    sigmar = ds.header['emission', 'frequency_param' ]

    if(init):
        plot_init()
        plt.gcf().subplots_adjust(top = 0.95, bottom = 0.05)
    
    plot_spectrum(dat, np.ones(len(dat)), **kwargs)

    fig   = plt.gcf()
    ax    = fig.axes[0]
    line  = ax.lines[-1]
    nbins = len(line.get_xdata())

    bm = fig.subplotpars.bottom
    fig.subplots_adjust(bottom = bm + 0.04)
    
    axsigma = plt.axes([0.2, bm - 0.03, 0.65, 0.03])
    ssigma  = Slider(axsigma, r"$\sigma_i$", sigmarng[0] * sigmar,
                     sigmarng[1] * sigmar, valinit=sigmar)
    
    def update(val):
        sigma = ssigma.val
        w = tlac_weights(ds, sigma, 0)
        x,y = plot_spectrum(dat, w, plot = False, **kwargs)
        line.set_xdata(x)
        line.set_ydata(y)
        ymax = ax.get_ylim()[1]
        if(np.max(y) > 1.1 * ymax):
            ax.set_ylim(ymax = 1.5 * ymax)
        elif(np.max(y) < 0.5 * ymax):
            ax.set_ylim(ymax = 0.5 * ymax)
        fig.canvas.draw_idle()
        

    ssigma.on_changed(update)

    # change current axis back
    plt.sca(ax)
    plt.show()
    
    # have to return ref to slider. otherwise unresponsive!
    return fig, ax, ssigma
    



def smooth_spectrum(x, y, width, window_type = 'gaussian',
                    convolve_mode = 'valid', window_len = 4.5):
    """
    Smooths spectrum with window function.
    
    Keyword Arguments:
    x             -- x-values of spectrum in chosen units
    y             -- y-values of spectrin in whatever units
    width         -- width of smoothing in x-units
    window_type   -- Shape of smoothing/window function (default 'gaussian').
                     Can be:
                        * 'gaussian' - Gaussian window. width is standard dev
    convolve_mode -- mode of `np.convolve`. Default: 'same'
    window_len    -- Length of window in units of `width` (will be modified
                     to be uneven)
    
    Returns:
    [x,y] of smoothed spectrum
    """
    Deltax = np.abs(x[1] - x[0])
    width_conv = float(width) / Deltax
    
    if window_type == 'gaussian':
        window_len = np.ceil(window_len * width_conv)
        window_len += window_len % 2 + 1  # make window length uneven
        window = gaussian(window_len, width_conv)
    else:
        raise ValueError("'window_type' can be 'gaussian'.")
        
    window /= np.sum(window)

    ycon = np.convolve(window, y, mode=convolve_mode)

    if convolve_mode == 'valid':
        s = 0.5 * (window_len - 1)
        xcon = x[s:-s]
    elif convolve_mode == 'same':
        xcon = x
    else:
        raise ValueError("convolve_mode other than 'valid' and 'same' "\
                         "not yet implemented")

    return xcon, ycon
        


def plot_projection(ds,  cmap = 'RdYlBu_r', labels = True,
                    background_color = 'white',
                    percentile_range = [2.5, 97.5],
                    colorbar = True):
    """Plots spacial projection of origin of photons.
    
    Keyword Arguments:
    ds                 -- TlacDat Set to plot
    labels             -- Whether or not to show axis labels (default: True)
    cmap               -- Color map to use (default 'RdYlBu_r')
    background_color   -- Background color (default 'white')
    percentile_range   -- Range of frequencies in percentile (default: [2.5,97.5])
    colorbar           -- Whether to show a frequency colorbar (default: True)
    """
    ax = plt.gca()
    
    # Get data
    try:
        box   = np.array([ float(i) * pc
                           for i in ds.header['grid','size'][1:-1].split(",") ])
    except: # spherical?
        box = np.array([0,0,0])
    freqs = ds.get_frequencies(phot_type='Lya', freq_type = 'v')
    k     = np.vstack( [ ds['Lya','dir0'], ds['Lya','dir1'],
                         ds['Lya','dir2'] ]).T
    pos   = np.vstack( [ ds['Lya','pos0'], ds['Lya','pos1'],
                         ds['Lya','pos2'] ]).T - 0.5 * box

    fmin, fmax = np.percentile(freqs,percentile_range)
    flim = np.max(np.abs([fmin,fmax]))

    xp = -np.sum(k * pos,axis=1)
    yp = np.sqrt( np.linalg.norm(pos,axis=1)**2 - xp**2) * np.sign(pos[:,2])

    # Color
    plt.scatter(xp / pc, yp / pc, c = freqs, marker='o', alpha=0.5,
                edgecolor='black', linewidth=0.15, cmap=cmap,
                vmin=-flim,vmax=flim)

    if labels:
        plt.xlabel(r"$-\vec{k}\cdot\vec{r}$ (pc)")
        plt.ylabel(r"$\sqrt{|\vec{r}|^2+"\
                   r"(\vec{k}\cdot\vec{r})^2}\frac{r_y}{|r_y|}$ (pc)")

    ax.patch.set_facecolor(background_color)

    if colorbar:
        plt.colorbar(label=r"$v$ (km/s)")
    
    
    
