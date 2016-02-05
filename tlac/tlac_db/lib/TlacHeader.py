import os, re
import numpy as np

from numbers import Number

class TlacHeader:
    
    def __init__(self, filename, replace_list = []):
        """Class to read & handle tlac header data (stored in cfg files for
        each run).

        Initialize:
        Initializes TlacHeader class. Loads "*cfg" file created
        by tlac from `filename`.

        `replace_list` is to replace the "<X>" templates that can be used in the
        tlac config files.
        """
        self.header = {}

        self.cfg_load(filename)
        self.add_default_values()
        self.cfg_replace_templates(replace_list)

    def __getitem__(self, (section, key)):
        return self.get(section, key)

    def keys(self):
        r = {}
        for sec in self.header:
            r[sec] = []
            for key in self.header[sec]:
                r[sec].append(key)

        return r

        
    def add_entry(self, section, key, value):
        if( (not section) or (not key)):
            print section, key
            raise Exception("Key or section not given!")

        if(not section in self.header):
            self.header[section.lower()] = {}

        self.header[section.lower()][key.lower()] = value

    
    def get(self, section, key):
        try:
            r = self.header[section.lower()][key.lower()]
        except:
            r = ""
            
        return r
        


    def cfg_load(self, filename):
        """Loads cfg file and sets header data.
        """
        if(not os.path.isfile(filename)):
            raise Exception("File " + filename + " does not exist.")
        
        f = open(filename)

        csec = ""
        
        for line in f:
            line = line.strip()
            if(self.cfg_valid_line(line)):
                if(line[0] == '['):
                    csec = line[1:].split("]")[0].lower()
                else:
                    key, value = self.cfg_parse_line(line)
                    self.add_entry(csec, key, value)


    def add_default_values(self):
        """Adds default values if they do not exist yet.
        """
        defaults = {'tlac' : {'turbulence' : 0.0}}

        for k, v in defaults.iteritems():
            for k1, v1 in v.iteritems():
                if k1 not in self.header[k]:
                    #print "adding default value",k,k1,v1
                    self.add_entry(k, k1, v1)
        

                    

    def cfg_replace_templates(self, replace_list):

        nrepl = 0
        for csec in self.header:
            for ckey in self.header[csec]:
                cval = self.header[csec][ckey]
                if(isinstance(cval, Number)):
                    continue
                    
                if(cval[0] == '<'):
                    n = int(cval[1:].split(">")[0])
                    if(n >= len(replace_list)):
                        msg = "[Error] '%s' found in header but only %d"\
                              " replacements provided" %(cval,
                                                         len(replace_list))
                        raise Exception(msg)

                    crepl = replace_list[n]
                    #print("%s:%s --> %s") %(csec, ckey, str(crepl))
                    try:
                        crepl = float(crepl)
                    except:
                        crepl = crepl
                    self.header[csec][ckey] = crepl
                    nrepl += 1

        if(nrepl != len(replace_list)):
            print("nrepl = %d, len=%d") %(nrepl, len(replace_list))
            raise Exception("Mismatch found.")
                


    def cfg_parse_line(self, line):
        """Parses cfg line
        Returns [key, value]
        """        
        key = line.split("=")[0]
        value = line[(len(key) + 1):].strip()
        value = value.split("#")[0].strip().lower()
        if( (value[0] == "\"") or (value[0] == "'")):
            value = value[1:-1]
        
        try:
            value = float(value)
        except:
            value = value
        key = key.strip().lower()
        
        return [key,value]
        
    def cfg_valid_line(self, line):
        """Checks if cfg line contains information
        """
        if( (not line) or (len(line) == 0) ):
            return False

        fc = line[0]
        if( (fc == '#') or (fc == ' ') or (fc == '\n')):
            return False

        return True
        

    

    def equal_parameters(self, to_compare, reverse_check = False,\
                         verbose = False, return_unequal = False,
                         ignore_secs = ('result','output','run')):
        """Checks if tlac runtime parameters of self and TlacHeader instance
        `to_compare` are the same.
        Returns True if they are, otherwise False

        If `reverse_check` is True checks also "the other way round".

        If `return_unequal` is True, instead of True/False, the list
        of unequal parameters is returned.

        `ignore_secs` are the sections to be ignored when comparing.
        """
        ul = []

        for csec in self.header:
            csec = csec.lower()
            if csec in ignore_secs:
                continue

            if(not csec in to_compare.header):
                if(verbose):
                    print("Section %s not found.") %csec
                if(return_unequal):
                    ul.append([csec])
                else:
                    return False

                
            for ckey in self.header[csec]:
                ckey = ckey.lower()
                myval = self.header[csec][ckey]

                if(not ckey in to_compare.header[csec]):
                    if(verbose):
                        print("Key %s not found.") %ckey
                    if(return_unequal):
                        ul.append([csec, ckey])
                    else:
                        return False

                if(self.header[csec][ckey] != to_compare.header[csec][ckey]):
                    # for numerical values
                    try:
                        f1 = float(self.header[csec][ckey])
                        f2 = float(to_compare.header[csec][ckey])
                        
                        if(f1 == f2):
                            continue
                        if((f2 != 0) and (np.abs(f1 / f2 - 1.0) < 1e-4)):
                            continue
                    except:
                        f1 = 0 #nothing

                    # for strings
                    if(isinstance(self.header[csec][ckey], str)):
                        s = []
                        s.append(self.header[csec][ckey])
                        s.append(to_compare.header[csec][ckey])
                        quotations = [ "\"", "\'"]
                        for i in range(2):
                            if(s[i][0:1] in quotations):
                                s[i] = s[i][1:-1]
                        if(s[0] == s[1]):
                            continue
                            
                        # for density identifier (e.g. "c 2e+21 == c 2e21")
                        try:
                            f1 = float(s[0][2:])
                            f2 = float(s[1][2:])
                                
                            if(f1 == f2):
                                continue
                            if((f2 != 0) and (np.abs(f1 / f2 - 1.0) < 1e-4)):
                                continue
                        except:
                            f1 = 0 #nothign

                        # for lists ...
                        if (s[0][0]==s[1][0]=='(' and s[0][-1]==s[1][-1]==')')\
                           and (s[0].count(",") == s[1].count(","))\
                           and (s[0].count(")") == s[1].count(")"))\
                           and (s[0].count(")") == s[1].count(")")):
                            l1 = 1
                            l2 = 2
                            # ... of lists
                            if s[0][1]==s[1][1]=='(' and \
                                    s[0][-2]==s[1][-2]==')':
                                l1 = [ [ float(i) for i in j[1:-1].split(",")]
                                       for j in s[0][1:-1].split(",") ]
                                l2 = [ [ float(i) for i in j[1:-1].split(",")]
                                       for j in s[1][1:-1].split(",") ]
                            # ... of numbers (e.g. (20.0,20.0) == (20, 20))
                            else:
                                try:
                                    l1 = [ float(i) for i in s[0][1:-1].split(",")]
                                    l2 = [ float(i) for i in s[1][1:-1].split(",")]
                                except:
                                    l1 = s[0][1:-1].split(",")
                                    l2 = s[1][1:-1].split(",")

                            if(l1 == l2):
                                continue
                            
                    if(verbose):
                        print("Values for key " + ckey + " dont agree. ("  +
                              str(self.header[csec][ckey]) + " vs " +\
                              str(to_compare.header[csec][ckey]) + ")")
                    if(return_unequal):
                        ul.append([csec, ckey])
                    else:
                        return False


        if(reverse_check):
            return to_compare.equal_parameters(self, reverse_check = False,
                                               verbose = verbose)
        else:
            if(return_unequal):
                return ul
            else:
                return True
                    
        
    def match_regexp(self, match_dict):
        """Checks if header matches key & vals in `match_dict`.
        Regexp can be used.
        Uses category:key notation.

        Example:
        hdr.match_regexp({"grid:dust" : "n 5.0+e-06"})

        Returns True or False.
        """
        matched_once = False
        regexp_flags = re.IGNORECASE
        
        for k, v in match_dict.iteritems():
            for sk_meta, sv_meta in self.header.iteritems():
                for sk, sv in sv_meta.iteritems():
                    sk = sk_meta + ":" + sk
                    if re.match(k, sk, regexp_flags) is not None:
                        matched_once = True
                        if re.match(v, sv, regexp_flags) is None:
                            return False

        if not matched_once:
            raise ValueError("The `match_dict` provided has not matched "
                             "once. Please check keys and don't forget the "
                             "category:key notation. E.g. hydrogen --> "
                             "grid:hydrogen")

        return True


    def get_v_thermal(self):
        """Returns thermal velocity in m/s.
        """
        cc = 8248.99 # k_b / m_hydrogen in m^2/s^2 / K
        T = self.header['grid']['temperature']
        
        return np.sqrt( 2. * cc * T)
        
