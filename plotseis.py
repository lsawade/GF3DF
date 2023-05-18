# %%
from obspy import read
from gf3d.source import CMTSOLUTION
from gf3d.seismograms import GFManager
from gf3d.plot.seismogram import plotseismogram
import matplotlib.pyplot as plt

#%% Get CMT source
# Actual values for half duration and timeshift
# time shift:     49.9800
# half duration:  33.4000

cmtfile = """
 PDEW2015  9 16 22 54 32.90 -31.5700  -71.6700  22.4 0.0 8.3 NEAR COAST OF CENTRAL CH
event name:     201509162254A
time shift:     0.00000
half duration:  0.00000
latitude:      -31.1300
longitude:     -72.0900
depth:          17.3500
Mrr:       1.950000e+28
Mtt:      -4.360000e+26
Mpp:      -1.910000e+28
Mrt:       7.420000e+27
Mrp:      -2.480000e+28
Mtp:       9.420000e+26
"""

cmt = CMTSOLUTION.read(cmtfile)

#%% Initialize the Green Function Manager

gfm = GFManager('single_element_not_fortran.h5')
gfm.load()

#%% Get Python seismograms

st = gfm.get_seismograms(cmt)

bfopy = st.select(station='BFO')

#%% Read seismograms from fortran
factor = 5
bfof = read('OUTPUT/II.BFO.*.sac')

for tr in bfof:
    tr.data *= factor

plotseismogram(bfopy, bfof, cmt)
plt.savefig('testplot.png', dpi=300)
