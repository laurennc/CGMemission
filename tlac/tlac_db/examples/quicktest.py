#! /usr/bin/env python

import os, sys, inspect, glob
import numpy as np

lib_path = os.path.abspath('../lib')
sys.path.append(lib_path)

from TlacHeader import TlacHeader
from TlacDatFile import TlacDatFile
from TlacDatSet import TlacDatSet
from TlacDB import TlacDB

print "Make some files..."
os.system("mkdir -p dat/")
os.system("../../../tlac_spherical example.cfg 100 124456 666 \"t 1e5\"")
os.system("../../../tlac_spherical example.cfg 100 12256 666 \"t 1e5\"")
os.system("../../../tlac_spherical example.cfg 100 12256 666 \"t 1e4\"")

print "========= Header tests ============="
cfgfiles = glob.glob("dat/*.cfg")
print cfgfiles

fn1 = cfgfiles[-1]
fn2 = cfgfiles[-2]

h1 = TlacHeader(fn1)

h2 = TlacHeader(fn2)

if(h1.equal_parameters(h2, verbose=True)):
    print "equal!"


print "========= dat tests ============="
datfiles = glob.glob("dat/*.dat")
print datfiles
d1 = TlacDatFile(datfiles[-1])
d2 = TlacDatFile(datfiles[-2])

if(d1.header.equal_parameters(d2.header)):
    print "equal parameters!"

d1.load_data(verbose = True)
print np.shape(d1.data)

print "========= dat set tests ============="
ds1 = TlacDatSet(verbose=True)
ds1.add_files(datfiles)

ds1.load_data()


print "========= db tests ============="
db1 = TlacDB(".")
db1.add(datfiles[-1])
db1.add(datfiles)
db1.summary()
db1.close()

db2 = TlacDB(".")
db2.close()

os.system("rm tlac_db.p")

