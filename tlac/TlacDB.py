"""TlacDB main file
"""

import os, sys, glob, psutil
import numpy as np

import logging
from prettytable import PrettyTable

#import pickle
import cPickle as pickle   ## speed up

from TlacDatSet import TlacDatSet

class TlacDB:
    def __init__(self, path, verbose = True,
                 dedicated = False, add_all = False,
                 memory_management = None):
        """Initializes TlacDB in `path`.

        Opens database if existing, creates if not.

        If `dedicated` is true, uses dedicated directories for data,
        i.e., moves data files around!

        If `add_all` is True and DB does not exist yet,
        tries to add all files in current directory.

        If `memory_management` is set to float in [0-100], then memory
        will be checked each time __getitem__ and emtied if its bigger
        than the percentage given.
        """
        if(not os.path.isdir(path)):
            msg = "\"%s\" is not a valid path!" % path
            raise Exception(msg)
        
        ## Set logging levels
        if(verbose):
            logging.basicConfig(stream=sys.stdout, level=logging.INFO,
                                format='%(message)s')
        
        ## Assign variables
        if(path[-1] == "/"):
            self.path = path
        else:
            self.path = path + "/"

        self.db = [] # TlacDatSet Objects
        self.dirs = [] # List of directories
        
        self.db_file = self.path + "tlac_db.p"

        self.dedicated = dedicated

        self.lookup_table = {}

        self.memory_management = memory_management
        
        ## Check if already db
        if(os.path.isfile(self.db_file)):
            self.open(verbose)
        else:
            self.create()
            if(add_all):
                self.add(glob.glob(self.path + "*.cfg"))
                self.save()
                if(verbose):
                    self.summary()
            
    def __getitem__(self, i):
        """Returns corresponding TlacDatSet to index `i`
        """
        if self.memory_management is not None:
            if psutil.virtual_memory().percent > self.memory_management:
                logging.info("Managing memory. Unloading database.")
                self.unload_all()
            
        return self.db[i]
        

    def create(self):
        """Creates new database.
        """
        logging.info("Create database "+ self.db_file)
            
        try:
            nothing = []
            pickle.dump( nothing, open( self.db_file, "wb" ) )
        except:
            raise Exception("Could not create file " + self.db_file)
        

            
    def open(self, verbose = True):
        """Opens database.
        """
        logging.info("Open database %s.", self.db_file)
        
        fn = self.db_file
        try:
            p = pickle.load(open(fn, "rb" ) )
        except:
            raise Exception("Error opening (pickle) file: " + fn)

        # backwards compatibility
        if not isinstance(p, dict):
            if isinstance(p, list):
                logging.info("List object found in pickle file. "
                             "Old DB file format?")
                self._open_from_filenames(p)
                return

            raise Exception("Unknown content in pickle file!")

        if len(p['db']) > 0:
            if not isinstance(p['db'][0], TlacDatSet):
                logging.warning("No TlacDatSet found in pickle file. "
                                "Try anyway to proceed.")
            
        if len(p['db']) != len(p['filenames']):
            msg = "Stored number of filenames and DatSets are not the same "\
                  "length. (%d vs %d). Load from "\
                  "filenames." %(len(p['filenames']),
                                 len(p['db']))
            logging.warning(msg)
            self._open_from_filenames(p)
            return
            
        self.db = p['db']
        logging.info("%d data sets loaded.", len(self.db))

        if(verbose):
            self.summary()

    def add(self, filenames):
        """Adds `filenames` to datbase where `filenames` are the config (*.cfg)
        files of the tlac output.

        Adds to given datset if already in db, otherwise creates a new entry.
        Moves if self.dedicated is true
        """
        if(isinstance(filenames, basestring)):
            filenames = [ filenames ]

        before = len(self.db)
        logging.info("Adding %d files to the db...", len(filenames))

        existing_files = [item for sublist in self.get_filenames()
                          for item in sublist] 

        for i, cf in enumerate(filenames):
            if cf in existing_files:
                continue
                
            tmp = TlacDatSet(filelst = cf)

            found = self.find(tmp.header)
            if(found):
                logging.info("\tParameters already in db. Adding to set.")
                if(not found.add_files(cf)):
                    logging.info("\tNot added. Same random seed?")
            else:
                logging.info("\tParameters not found. Making new entry.")
                self.db.append(tmp)

        logging.info("...done! %d files added.", len(self.db) - before)




    def save(self):
        """Saves db to pickle file.
        """
        logging.info("Save database.")
        save_dict = {
            'filenames' : self.get_filenames(),
            'db' : self.db
        }
        pickle.dump( save_dict, open( self.db_file, "wb" ), -1 )
        logging.info("\t%d dat sets written", len(self.get_filenames()))
        
            
    def close(self):
        """Saves  & closes db
        """
        self.save()
    
    def summary(self, print_uv = 'auto', extra_columns = [],
                print_nphot = True, print_ndat = True):
        """Prints summary/statistics of db.

        print_uv decides wheter or not uv stuff gets printed.
        'auto' does not print if there is no UV photons
         
        extra_columns can be a list of [section, key] pairs
        """

        # Check parameters
        if(print_uv == 'auto'):
            n = self.get_ndat("UV")
            if(n == 0):
                print_uv = False
            else:
                print_uv = True
        elif(isinstance(print_uv, str)):
            raise Exception("Unrecognized string " + print_uv)


        # Make table
        tbl = PrettyTable()

        tbl.add_column("#", range(len(self.db)))
        if(print_nphot):
            tbl.add_column("nphot_lya", [i.get_nphot("Lya") for i in self.db])
        if(print_ndat):
            tbl.add_column("ndat_lya", [i.get_ndat("Lya") for i in self.db])

        if(print_uv):
            if(print_nphot):
                tbl.add_column("nphot_uv", [i.get_nphot("UV") for i in self.db])
            if(print_ndat):
                tbl.add_column("ndat_uv", [i.get_ndat("UV") for i in self.db])

        for csec, ckey in extra_columns:
            tbl.add_column(ckey, [i.header.get(csec, ckey) for i in self.db])

        print tbl
        


    def find(self, header, verbose = False, force = False,
             return_id = False, **kwargs):
        """Searches for TlacHeader `header` in database and returns entry to
        TlacDatSet if found.

        Else returns False

        If `force` is true, throws exception if not found!

        If `return_id` is true, not the TlacDatSet but the id of it will be
        returned.

        """
        for n, cds in enumerate(self.db):
            if(verbose):
                print "===>%d<===" %n
            if(cds.header.equal_parameters(header, verbose = verbose,
                                           **kwargs)):
                if(return_id):
                    return n
                else:
                    return cds

        if(force):
            raise Exception("No entry in database found!")
        
        return False


    def select(self, locate_dict, verbose = False):
        """Selects DataSets that are conformal with `locate dict`.

        Keywords:
          * locate_dict    -- Dictionary of header keys & corresponding
                              searching term. Regexp can be used for the
                              keys as well as the values!
          * verbose        -- Prints info if True (default: False)

        Returns:
          List of TlacDatSets which fit to the search pattern.

        Example:
          # Selects all the datasets with exacly 5e-5 hydrogen number density
          s = db.select({'grid:hydrogen' : 'n 5e-5'})
          # This probably won't work since the search string has to match the
          # one given in the header (which is something like "n 5.000e-05)

          # A wider selection criterion helps
          s = db.select('*hydrogen' : 'n 5\.?0+?e-0+?5')
        """

        if not isinstance(locate_dict, dict):
            raise ValueError("`locate_dict` has to be a dictionary!")

        r = []
        for cds in self.db:
            if cds.header.match_regexp(locate_dict):
                r.append(cds)


        return r
            


    def remove(self, ids, remove_files = False):
        """Removes DatSet from DB
        
        Keyword Arguments:
        ids            -- ID of TlacDatSet
        remove_files -- If True, removes also corresponding files (default not)
        """
        if remove_files:
            fns = [ glob.glob(os.path.splitext(i)[0] + ".*")
                    for i in  self.db[ids].get_filenames() ]
            fns = [item for sublist in fns for item in sublist] # flatten
            logging.info("Removing %d files." %(len(fns)))
            [ os.remove(i) for i in fns ]

        self.db.pop(ids)

        
    def get_filenames(self):
        """Returns list of lists of tlac dat filenames.
        (one entry per DatSet)
        """
        r = []
        for cds in self.db:
            r.append([os.path.basename(i) for i in cds.get_filenames()])

        return r


    def get_header_values(self, sec, key):
        """Returns list of header values in section `sec` and key `key`
        """
        return [ i.header.get(sec, key) for i in self.db ]
        

    def get_nphot(self, phot_type = "Lya"):
        return np.sum([ i.get_nphot(phot_type) for i in self.db ] )


    def get_ndat(self, phot_type = "Lya"):
        return np.sum([ i.get_ndat(phot_type) for i in self.db ] )

    def load_all(self, phot_type = 'both', verbose = False):
        """Loads all photons of `phot_type`.
        phot_type can be 'both', 'Lya' or 'UV'.
        """
        l = len(self.db)
        for i,cds in enumerate(self.db):
            if(verbose):
                print("Loading %d/%d..." %(i, l))
            cds.load_all(phot_type)
            if(verbose):
                print("done loading (%d/%d)!" %(i, l))


    def unload_all(self):
        for cds in self.db:
            cds.unload_all()

    def load(self, **kwargs):
        print "Deprecated function db.load. Use load_all or e.g. something "\
            "like db[0]['Lya','x']."

        
    def create_lookup_table(self, header_keys, callf = None, store_ids = False,
                            **kwargs):
        """Populates self.lookup_table with dataset header values by keys.

        Arguments:
        header_keys  -- List of lists contaning TlacHeader [section, key]-pairs
        callf        -- Function to be called with header values of header &
                        kwargs. Has to return a hashable object. If None is
                        given (default), a tuple of values is used.
        store_ids    -- If True not a pointer to a DatSet but the ID of it is
                        stored (default: False)
        
        Can be used to speedup search.
        """ 
        self.lookup_table = {}

        for i,cds in enumerate(self.db):
            vals = [ cds.header[sec,key] for sec,key in header_keys ]
            if(store_ids):
                l = i
            else:
                l = cds
            if callf is None:
                self.lookup_table[tuple(vals)] = l
            else:
                self.lookup_table[callf(vals, **kwargs)] = l

            
    def lookup(self, values):
        """Finds DatSet in the lookup table created with `create_lookup_table`.

        `values` must be in the same order as `header_keys` before.
        """
        if(self.lookup_table == {}):
            raise Exception("Call `create_lookup_table` first!")

        #if not isinstance(values, tuple):
        #    values = tuple(values)

        return self.lookup_table[values]


    def _open_from_filenames(self, filenames):
        """Opens database given filenames.
        """
        logging.info("Try opening from filenames")

        for i, cfns in enumerate(filenames):
            self.db.append(TlacDatSet())
            cfns_mod = [ self.path + cc  for cc in cfns ]
            self.db[-1].add_files(cfns_mod, check=False)

        logging.info("%d data sets loaded.", len(self.db))
