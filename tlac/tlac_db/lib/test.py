import TlacASCII
import TlacHDF5
import TlacDat
import TlacDatSet

import TlacDB

import glob

ta = TlacASCII.TlacASCII("test.dat")

ta.load("pos0")
ta.load("x")

ta['x']
ta['pos1']


print "HDF5"
th = TlacHDF5.TlacHDF5("test.h5", "Lya")
print th['x']
print th['pos1']



print "TlacDat test"
fdir = "/Users/maxbg/Dropbox/uni_aktuell/1313_phd/tlac/source/tests/dat/"
fnlst = glob.glob(fdir + "output*.cfg")

for fn in fnlst:
    print fn
    td = TlacDat.TlacDat(fn)
    print td['Lya']
    print td['Lya']['x']
    print td['Lya','x']


print "TlacDatSet test"
tds = TlacDatSet.TlacDatSet(fdir + "output", verbose = True)


print "TlacDB test"
tdb = TlacDB.TlacDB(fdir, verbose = True)
tdb.add(fnlst)
tdb.summary()
