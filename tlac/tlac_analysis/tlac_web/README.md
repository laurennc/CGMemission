# tlac_web
Version 0.3 (date 23.3.2015)

## Usage
Simply play with the sliders on the LHS to see how the Lya spectrum changes.
You can push the buttons "0"-"4" to add new lines as comparison. Below the plot
you will find information on what is shown currently and on the right hand side
there is general information about the database as well as an option to change
the frequency units (where x = (nu - nu_Lya)/Delta nu and v is the velocity
shift in km/s).
Direclty above the plot, there is also a toolbox containing zoom etc.
provided directly by [bokeh](http://bokeh.pydata.org) (more help on this if you
click the "?" icon).

### Upload/Download data
When tlac_web is run embedded in a [flask](http://flask.pocoo.org) application,
it is possible to download the shown data & upload own spectra to compare with.

To download data, press the "save" button to save as many spectra you like
(each time, all the shown spectra are saved and the file counter should go up).
Then, press the "Download generated files" link on top of the page to download
all the created spectra in a tar archive.

In order to upload own spectra, press the "Upload data" link on the top of the
page. Then select a file (allowed extensions are "dat", "tsv" and "txt") which
contains "x" and "y" column or "x" and "y" row and press "upload".
Afterwards, the file should appear in the drop-down menu in the right column. If
not, simply change the parameters (i.e. a slider) a little bit to update the
list.

## Physics behind it
The setup is a simple "shell model", i.e., a central source surrounded
by an outflowing shell of hydrogen.
We ran the radiative transfer simulation for several hydrogen
densities, gas velocities & temperatures.
The other parameters (dust content, width of intrinsic spectrum,
intrinsic equivalent width) are post-processed.


## Parameters
* tau_d  -- optical depth of (all-absorbing) dust
* sigma_i -- width of intrinsic spectrum in km/s
* EW_i -- intrinsic equivalent width in Angstrom. Zero means no continuum.
* hydrogen -- log10 of column density in cm^-2
* velocity_params -- outflowing velocity in km/s
* temperature -- log10 of temperature in K
* nbins -- number of bins, auto is chosen by the Friedmann-Diaconis rule
* nsmooth -- Possibility to smooth the resulting spectra through a convolution
with a Gaussian where the standard deviation is `nsmooth` in the current x-axis
units

## About
Author: Max Gronke (maxbg -sign-used-in-emails- astro.uio.no)

Any feedback/comments are welcome!

* Plotting is done using [bokeh](http://bokeh.pydata.org) (Copyright (c) 2013,
Continuum Analytics, Inc. All rights reserved.)
* The web-frontend is powered by [flask](http://flask.pocoo.org/) (Copyright (c)
2014 by Armin Ronacher)

