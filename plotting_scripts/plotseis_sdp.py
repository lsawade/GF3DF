# %%
import os
import sys
from obspy import read
from gf3d.source import CMTSOLUTION
from gf3d.seismograms import GFManager
from gf3d.plot.seismogram import plotseismogram
from gf3d.plot.frechet import plotfrechet
from gf3d.process import process_stream
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# %% Get CMT source
# Actual values for half duration and timeshift
# time shift:     49.9800
# half duration:  33.4000

cmt = CMTSOLUTION.read('CMTSOLUTION')

# %% Initialize the Green Function Manager

gfm = GFManager('single_element_not_fortran.h5')
gfm.load()

# %% Get Python seismograms

rp = gfm.get_seismograms(cmt)
dp = gfm.get_frechet(cmt, rtype=3)

pdppy = dict()
for key, _st in dp.items():
    pdppy[key] = _st.select(station='BFO')

# bfopy.differentiate()

# pbfopy = process_stream(bfopy, cmt=cmt)

# %% Read seismograms from fortran
pypars = ["Mrr", "Mtt", "Mpp", "Mrt", "Mrp",
          "Mtp", "latitude", "longitude", "depth", "time_shift", "hdur"]


pfp = read('OUTPUT_SDP/II.BFO.MX?.sem.sac')
pdfp = dict()

for _i, par in enumerate(
    ["mrr", "mtt", "mpp", "mrt", "mrp", "mtp",
     "lat", "lon", "dep", "cmt", "hdr"]):
    pdfof = read(f'OUTPUT_SDP/II.BFO.*.{par}.sem.sac')
    pdfp[pypars[_i]] = pdfof.select(station='BFO')

starttime = pfp[0].stats.starttime + 200
endtime = starttime + 10000
limits = (starttime.datetime, endtime.datetime)

plotfrechet(cmt, rp, dp, 'II', 'BFO', comp='z',
            limits=limits, outdir=".", lw=0.25,
            rp2=pfp, drp2=pdfp)

# factor = 1
# for tr in bfof:

#     tr.data *= factor

# if all([os.path.exists(f'OUTPUT_WIN/II.BFO.MX{C}.sem.sac') for C in ['N', 'E', 'Z']]):
#     pbfof_win = read('OUTPUT_WIN/II.BFO.*.sac')
#     # pbfof_win = process_stream(bfof_win, cmt=cmt)
# else:
#     pbfof_win = None

# pbfof = process_stream(bfof, cmt=cmt)
# plotseismogram(pbfopy, pbfof, cmt, newsyn=pbfof_win, nooffset=True, lw=0.25)
# plt.savefig('testplot_nooffset.pdf', dpi=300)
# plotseismogram(pbfopy, pbfof, cmt, newsyn=pbfof_win, nooffset=False, lw=0.25)
# plt.savefig('testplot.pdf', dpi=300)

# plt.close('all')
