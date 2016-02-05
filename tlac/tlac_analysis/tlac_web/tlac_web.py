"""
Using bokeh (http://bokeh.pydata.org/) to visualize a TlacDB on the web

There is two ways to run tlac_web:
- via the bokeh-server    -- see `bokeh_server.py`
- embedded in a flask app -- see `flask_app.py`

"""

################################################################################
## Imports
try:
    import tlac_analysis as ta
except:
    import os
    import sys
    myp = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                 '../../'))
    sys.path.append(myp)
    import tlac_analysis as ta

from TlacDB import TlacDB
from TlacHeader import TlacHeader

import sys, copy, time, os, glob
from prettytable import PrettyTable
import logging
import numpy as np

from bokeh.plotting import figure
from bokeh.models import Plot, ColumnDataSource, Range1d, LinearAxis
from bokeh.properties import Instance, String
from bokeh.models.widgets import (HBox, Slider, TextInput, VBoxForm, Select,
                                  PreText, Paragraph, TableColumn, VBox,
                                  RadioButtonGroup, Toggle)
from bokeh.models.glyphs import Line

from db_helpers import *

################################################################################
## Config
LOG_FORMAT = '[tlac_web] %(asctime)-15s %(message)s'
logging.basicConfig(format=LOG_FORMAT, level=logging.INFO)

db = TlacDB(".", verbose=True, add_all = True)

## Replace this line by automatic picking up of differences
db.choice_keys = ['hydrogen', 'velocity_params', 'temperature' ]
db.choice_type = [ -5, -5, -5] # 0 = Select, 1 = Slider, <0 = dynamic

## Keys of post-processed values
db.post_keys = [ "sigma_i", "tau_d", "EW_i", 'nbins', 'nsmooth' ]

## how many lines to choose from
nlines = 4
colorlst = ["blue", "red", "green", "orange"]

## Version
version = "0.03"


################################################################################
## Global stuff

# Make lookup table
k = [ ["grid", i] for i in db.choice_keys ]
db.create_lookup_table(k, hdr_vals_to_floats, store_ids = True,
                       return_tuple = True, db = db)

## Assemble lists to choose from
db.choice_lst = []
for i, ck in enumerate(db.choice_keys):
    db.choice_lst.append(np.sort(list(set([ j[i] for j in db.lookup_table ]))))

## Assign correct type if "dynamic" was chosen
for i, cl in enumerate(db.choice_lst):
    if db.choice_type[i] < 0:
        if len(cl) < -db.choice_type[i]:
            db.choice_type[i] = 0
        else:
            db.choice_type[i] = 1
    

db.description_text = db_description_text(db)


################################################################################
### Class stuff

class TlacWebApp(VBox):
    extra_generated_classes = [["TlacWebApp", "TlacWebApp", "VBox"]]

    inputs = Instance(VBoxForm)
    firstrow = Instance(HBox)
    
    cdsid = Instance(PreText) # stores ID of current data set
    values = Instance(PreText) # stores values selected
    dbdesc = Instance(PreText) # info about DB
    
    tau_d   = Instance(Slider)
    sigma_i = Instance(Slider)
    EW_i = Instance(Slider)

    nbins = Instance(Slider)
    nsmooth = Instance(Slider)

    lineselect = Instance(RadioButtonGroup)
    xaxisselect = Instance(RadioButtonGroup)

    lbl_xaxisselect = Instance(Paragraph)

    rightcol = Instance(VBoxForm)

    dl_button = Instance(Toggle)
    lbl_dl_button = Instance(Paragraph)

    ul_select = Instance(Select)
    lbl_ul_select = Instance(Paragraph)
    
    # shit solution but tried many other things:
    # dictionary, list, attribute of other functions ....
    selects0 = create_input_instance(db, 0)
    selects1 = create_input_instance(db, 1)
    selects2 = create_input_instance(db, 2)
    selects3 = Instance(Select)
    selects4 = Instance(Select)
    selects5 = Instance(Select)

    steps = []
    
    plot = Instance(Plot)
    source = Instance(ColumnDataSource)

    def __init__(self, *args, **kwargs):
        super(TlacWebApp, self).__init__(*args, **kwargs)

        # Initial data for data source object
        d = {}
        for i in range(nlines):
            d['x' + str(i)] = [] # x-values
            d['y' + str(i)] = [] # y-values
            d['p' + str(i)] = [] # parameters
        self._init_data = d

        # List of all keys
        self._parameter_keys = [ "#", "cdsid", "color"] + \
                               db.choice_keys + db.post_keys


    @classmethod
    def create(cls, data_prefix = None, upload_prefix = None):
        """
        This function is called once, and is responsible for
        creating all objects (plots, datasources, etc)

        Keywords:
        If `data_prefix` is given (default None), a button to generate data is
        shown. Data will then be placed in "`data_prefix`_X.dat" where X is
        counting upwards.

        If `upload_prefix` is given (default not), lets the user show uploaded
        data with this prefix.
        """
        obj = cls()

        obj.cdsid = PreText(width=30,height=6,text="-1")

        obj.values = PreText(width=1000,height=150,text="(please select values)")

        obj.dbdesc =  PreText(width=200,height=300,text=db.description_text)

        obj.tau_d = Slider(
            title="tau_d", name='tau_d',
            value=0.0, start=0.0, end=4
        )

        a,b = np.round(get_sigma_i_range(db))
        obj.sigma_i = Slider(
            title="sigma_i", name='sigma_i',
            value= 0.5 * (a+b), start=a, end=b
        )

        obj.EW_i = Slider(
            title="EW_i (0 = cont. off)", name='EW_i',
            value=0, start=0, end=100
        )

        obj.nbins = Slider(
            title="nbins (0=auto.)", name='nbins',
            value=0, start=0, end=200
        )

        obj.nsmooth = Slider(
            title="smooth (0=none)", name='nsmooth',
            value=0, start=0, end=100
        )

        s = [ str(i) for i in range(nlines) ]
        obj.lineselect = RadioButtonGroup(labels=s, active=0)

        obj.xaxisselect =  RadioButtonGroup(labels= ["x", "v"], name="x-axis",
                                            active=0)

        obj.lbl_xaxisselect = Paragraph(text="x-axis:")
        
        # Create dropdown menus
        #obj.selects = {}
        for i, ck in enumerate(db.choice_keys):
            lst = db.choice_lst[i]
            if (db.choice_type[i] == 0):
                lst = [ str(j) for j in lst ]
                setattr(obj,"selects"+str(i), Select(title=ck,
                                                     name = ck, value = lst[0],
                                                     options = ["(none)"] + lst))
            elif db.choice_type[i] == 1:
                if(len(lst) < 5):
                    raise Warning("Does a slider make sense with < 5 values?")
                    
                #steps = np.round(lst[1:] - lst[:-1],3)
                #if(len(set(steps)) != 1):
                #    print "Unique Steps:",set(steps)
                #    raise Exception("Unequal spacings between parameters."
                #                    "Not supported yet.")
                #print "Allsteps", steps
                setattr(obj,"selects"+str(i), Slider(title=ck,
                                                     name = ck,
                                                     start = lst[0] * 0.999,
                                                     end = lst[-1] * 1.001,
                                                     # Problem here???
                                                     # not taking float steps...?
                                                     # TODO: Hack it!
                                                     # step = steps[0],
                                                     value = lst[0]))
        # Create bokeh source object
        obj.source = ColumnDataSource(data=obj._init_data)

        toolset = "pan,reset,resize,save,wheel_zoom,box_zoom"
        obj.plot = figure(title_text_font_size="12pt",
                          plot_height=600,
                          plot_width=600,
                          tools=toolset,
                          title="tlac_web (first choose from all the menus)"
                          #x_mapper_type='Photon frequency (units see on the right)',
                          #y_mapper_type='Intensity'
        )
        obj.plot.line('x0', 'y0', source=obj.source,
                      line_width=3, line_alpha=0.6,
                  )
        
        ## Set axis labels --> Doesn't work since x/y order depends on machine
        obj.plot.xaxis.axis_label = "Frequency (x)"
        obj.plot.xaxis.name = "x_axis"
        obj.plot.yaxis.axis_label = "Intensity" 

        for iline in range(1,nlines):
            l = Line(x = "x" + str(iline), y = "y" + str(iline),
                     line_color = colorlst[iline], line_width=3,
                     line_alpha = 0.6)
            obj.plot.add_glyph(obj.source, l)


        # Download possibility
        if data_prefix is not None:
            obj.dl_button = Toggle(
                name="dl_button",
                label="Save",
                type = "primary"
            )
            obj.source.data['data_prefix'] = data_prefix
            obj.lbl_dl_button = Paragraph(text="0 files created.")
        else:
            obj.dl_button = None

        # Upload possibility
        if upload_prefix is not None:
            obj.ul_select = Select(
                name = "ul_select",
                title = "Uploaded file",
                options = ["(none)"]
            )
            obj.source.data['upload_prefix'] = upload_prefix
            l = Line(x = "x_ul", y = "y_ul",
                     line_color = 'black', line_width=3,
                     line_alpha = 0.6)
            obj.plot.add_glyph(obj.source, l)

        else:
            obj.ul_select = None
            
        # Really stupid....
        objselectlst = [ obj.selects0, obj.selects1, obj.selects2,\
                         obj.selects3, obj.selects4, obj.selects5 ]
        
        obj.inputs = VBoxForm(
            children=[
                obj.cdsid, obj.tau_d, obj.sigma_i, obj.EW_i
            ] + objselectlst[:len(db.choice_keys)] + [ obj.nbins, obj.nsmooth,
                                                       obj.lineselect ]
        )

        l = [ obj.dbdesc, obj.lbl_xaxisselect, obj.xaxisselect ]
        if obj.dl_button is not None:
            l += [ obj.lbl_dl_button, obj.dl_button ]
        if obj.ul_select is not None:
            l += [ obj.ul_select ]

        obj.rightcol = VBoxForm(
            children=l
        )

        obj.firstrow = HBox(children=[obj.inputs,obj.plot, obj.rightcol])

        obj.children.append(obj.firstrow)
        obj.children.append(obj.values)

        logging.info("New form created.")

        return obj

    def setup_events(self):
        super(TlacWebApp, self).setup_events()
        if not self.sigma_i: # to check if up an running?
            return
            
        # Register function called on change
        allwidgets = [ 'sigma_i', 'tau_d', 'EW_i', 'nbins', 'nsmooth' ] +\
                     [ 'selects' + str(i) for i in range(len(db.choice_keys)) ]
        for w in allwidgets:
            getattr(self, w).on_change('value', self, 'input_change')

        self.xaxisselect.on_click(self.update_xaxis)

        if self.dl_button is not None:
            self.dl_button.on_click(self.create_data_file)

        if self.ul_select is not None:
            self.ul_select.on_change('value', self, 'overplot_uploaded')
        
        self.update_dataset()
        self.update_parameter_values()

        
    def input_change(self, obj, attrname, old, new):
        """
        This callback is executed whenever the input form changes. It is
        responsible for updating the plot, or anything else you want.
        The signature is:

        Args:
            obj : the object that changed
            attrname : the attr that changed
            old : old value of attr
            new : new value of attr

        """
        self.values.text = "(...loading...)"
        self.update_parameter_values()
        
        msg = "Changed %s from %g to %g" %(obj.name, old, new)
        logging.info(msg)
        if( obj in [self.sigma_i, self.tau_d, self.nbins, self.nsmooth]):
            if obj == self.nsmooth:
                if new < 0.1 * self.nsmooth.end:
                    self.nsmooth.value = 0
            self.update_data(int(self.cdsid.text))
        else:
            if obj.name in db.choice_keys:
                # Step hack
                i = db.choice_keys.index(obj.name)
                if(db.choice_type[i] == 1):
                    l = db.choice_lst[i]
                    v = new
                    #print "Try to find closest to ", v
                    m = (np.abs(l - v) == np.min(np.abs(l - v)))
                    #print "Locked in ", l[m][0]
                    obj.value = l[m][0]
            self.update_dataset()
            self.update_data(int(self.cdsid.text))

        if self.ul_select is not None:
            self.update_uploaded_data_list()
        
        self.update_values_text()

        
    def update_parameter_values(self):
        """Updates values in the source.data object (prefix 'p')
        """
        nselect = 0
        v = []
        for j,ckey in enumerate(self._parameter_keys):
            if ckey == "cdsid":
                v.append(int(self.cdsid.text))
            elif ckey in db.choice_keys:
                v.append(getattr(self,"selects"+str(nselect)).value)
                nselect += 1
            elif ckey in db.post_keys:
                v.append(getattr(self,ckey).value)
            else:
                if ckey == "#":
                    v.append(self.lineselect.active)
                elif ckey == "color":
                    v.append(colorlst[self.lineselect.active])
                else:
                    raise ValueError("Weird key...?" + ckey)
                
        self.source.data['p' + str(self.lineselect.active) ] = v
        
        
    def update_values_text(self):
        """ Updates info table in the text box below the plot.
        """
        p = PrettyTable(self._parameter_keys)
        for i in range(nlines):
            v = self.source.data['p' + str(i)]
            if(v != []):
                p.add_row(v)
        
        self.values.text = p.get_string()

        
    def update_dataset(self):
        values = [ (getattr(self,"selects"+str(i))).value
                   for i in range(len(db.choice_keys)) ]
        self.cdsid.text = str(find_dsid(db, db.choice_keys, values))

            
    def update_data(self, cdsid, iline = None):
        """Updates data (x,y) and plots it
        """
        if iline is None:
            iline = int(self.lineselect.active)
        if(cdsid < 0):
            return
        x_axis = self.xaxisselect.labels[self.xaxisselect.active]
        d = self.source.data

        p = {}
        for i,k in enumerate(self._parameter_keys):
            p[k] = d['p' + str(iline)][i]
        neff, x, y = get_data(db, cdsid, x_axis,
                              p['sigma_i'],
                              p['tau_d'], p['EW_i'],
                              p['nbins'], p['nsmooth'])
        title = "tlac_web " + version
        self.plot.title = title
        
        self.source.data = dict()
        d['x' + str(iline)] = list(x)
        d['y' + str(iline)] = list(y)
        self.source.data = d 


    def update_xaxis(self, new):
        """Updates x values of existing curves to "new".
        """
        ## Also does not work...
        new_str = self.xaxisselect.labels[new]
        print new, new_str
        if new_str == 'v':
            newlbl = "Frequency (" + new_str + " in km/s)" 
            self.nsmooth.end = 100
        elif new_str == 'x':
            newlbl = "Frequency (" + new_str + ")" 
            self.nsmooth.end = 20


        ## Axis label
        ax =  self.plot.select(dict(name="x_axis"))
        ax[0].axis_label = newlbl

        
        #print "Current DS ids:"
        d = self.source.data
        for k in d.keys():
            if k[0] == 'p' and d[k] != []:
                cdsid = d[k][1]
                iline = d[k][0]
                #print d[k]
                #print cdsid, iline
                self.update_data(cdsid, iline)
        


    def create_data_file(self, val):
        """
        Creates downloadable data file with filename "`dl_prefix`_N.dat"
        """
        prefix = self.source.data['data_prefix']
        assert prefix is not None

        maxfiles = 20
        self.lbl_dl_button.text = "Saving..."
        
        for nfile in range(maxfiles):
            cname =  prefix + "_" + str(nfile) + ".dat"
            if not os.path.isfile(cname):
                break
        if nfile == maxfiles:
            # All files created
            logging.info("All files are created with prefix %s (i=%d)" %(prefix,
                                                                         nfile))
            self.lbl_dl_button.text = "Maximum number of files reached!"
            return

        ## Generate data to output
        try:
            f = open(cname, "w")
        except:
            msg = "Error creating file %s" %cname
            logging.error(error)
            return
        logging.info("Write data output to %s" %(cname))
        f.write("##########################################################\n")
        f.write("# Data created " + time.ctime() + "\n")
        f.write("# with tlac_web " + version + "\n")
        f.write("# Author: Max Gronke (maxbg~~astro.uio.no)\n")
        f.write("# Any comments/feedback is appreciated.\n")
        f.write("##########################################################\n")
        f.write("####### Database description #############################\n")
        f.write("# " + self.dbdesc.text.replace("\n","\n# "))
        f.write("##########################################################\n")
        f.write("####### Spectra parameters ###############################\n")
        f.write("# " + self.values.text.replace("\n","\n# ") + "\n")
        f.write("##########################################################\n")
        f.write("######## Data ############################################\n")
        d = self.source.data
        for iline in range(nlines):
            if d['x' + str(iline)] != []:
                x = d['x' + str(iline)]
                y = d['y' + str(iline)]
                f.write("# x" + str(iline) + "\n")
                f.write(' '.join(map(str, x)) + "\n")
                f.write("# y" + str(iline) + "\n")
                f.write(' '.join(map(str, y)) + "\n")
        f.write("##########################################################\n")
        f.close()
        logging.info("Done writing data output to %s" %(cname))

        self.lbl_dl_button.text = "%d files created" %(nfile + 1)
        
        

    def overplot_uploaded(self, obj, attrname, old, new):
        """When user selects a new file on the uploaded file list, we show it.
        """
        if old == new:
            return
        
        fn =  self.source.data['upload_prefix'] + new
        if not os.path.isfile(fn):
            new_x = []
            new_y = []
        else:
            try:
                dat = np.loadtxt(fn)
                s = np.shape(dat)
                if s[1] == 2: # columnwise
                    new_x = list(dat[:,0])
                    new_y = list(dat[:,1])
                else: #rowwise
                    new_x = list(dat[0,:])
                    new_y = list(dat[1,:])
            except:
                new_x = []
                new_y = []
                
        d = self.source.data
        self.source.data = dict()
        d['x_ul'] = new_x
        d['y_ul'] = new_y
        self.source.data = d
        

    def update_uploaded_data_list(self):
        """Updates the file list shown to the user to select from his uploaded
        data.
        """
        prefix = self.source.data['upload_prefix']
        
        filelst = [ i[len(prefix):] for i in glob.glob(prefix + "*") ]

        newlst = ["(none)"] + filelst
        
        if self.ul_select.options != newlst:
            self.ul_select.options = newlst

        

        
