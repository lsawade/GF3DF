# %%
import os
import sys
from obspy import read
from gf3d.source import CMTSOLUTION
from gf3d.seismograms import GFManager
from gf3d.plot.seismogram import plotseismogram
from gf3d.process import process_stream
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# %% Get CMT source
# Actual values for half duration and timeshift
# time shift:     49.9800
# half duration:  33.4000

cmtfile = """
 PDEW2015  9 16 22 54 32.90 -31.5700  -71.6700  22.4 0.0 8.3 NEAR COAST OF CENTRAL CH
event name:     201509162254A
time shift:      50.00000
half duration:  33.40000
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

# %% Initialize the Green Function Manager

# gfm = GFManager('single_element_not_fortran.h5')
gfm = GFManager('sample.h5')
gfm.load()

# %% Get Python seismograms

st = gfm.get_seismograms(cmt)

bfopy = st.select(station='BFO')

# bfopy.differentiate()

pbfopy = process_stream(bfopy, cmt=cmt)

# %% Read seismograms from fortran
bfof = read('OUTPUT/II.BFO.*.sac')

# factor = 1
# for tr in bfof:

#     tr.data *= factor

if all([os.path.exists(f'OUTPUT_WIN/II.BFO.MX{C}.sem.sac') for C in ['N', 'E', 'Z']]):
    bfof_win = read('OUTPUT_WIN/II.BFO.*.sac')
    pbfof_win = process_stream(bfof_win, cmt=cmt)
else:
    pbfof_win = None

pbfof = process_stream(bfof, cmt=cmt)
plotseismogram(pbfopy, pbfof, cmt, newsyn=pbfof_win, nooffset=True, lw=0.25)
plt.savefig('testplot_nooffset.pdf', dpi=300)
plotseismogram(pbfopy, pbfof, cmt, newsyn=pbfof_win, nooffset=False, lw=0.25)
plt.savefig('testplot.pdf', dpi=300)

plt.close('all')


# %%
stf = read('OUTPUT/STF.ERF.TIM.sem.sac')
plt.figure()
plt.close('all')
fig = plt.figure(figsize=(7, 7))

ax = plt.subplot(3, 1, 1)
plt.plot(stf[0].times("matplotlib"), stf[0].data,
         '-', c='k', lw=0.5, label='STF')
plt.xlabel('Time')
plt.ylabel('A')

ax.xaxis_date()
ax.xaxis.set_major_formatter(
    mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
ax.tick_params(labelleft=True, left=False)
ax.spines.right.set_visible(False)
# ax.spines.left.set_visible(False)
ax.spines.top.set_visible(False)

fstf = read('OUTPUT/STF.ERF.FRE.sem.sac')
ax = plt.subplot(3, 1, 2)
plt.plot(fstf[0].times("matplotlib"), fstf[0].data,
         '-', c='k', lw=0.5, label='STF')
plt.xlabel('Time')
plt.ylabel('A')

ax.xaxis_date()
ax.xaxis.set_major_formatter(
    mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
ax.tick_params(labelleft=True, left=False)
ax.spines.right.set_visible(False)
# ax.spines.left.set_visible(False)
ax.spines.top.set_visible(False)

istf = read('OUTPUT/STF.ERF.IFR.sem.sac')
ax = plt.subplot(3, 1, 3)
plt.plot(istf[0].times("matplotlib"), istf[0].data,
         '-', c='k', lw=0.5, label='STF')
# plt.plot(istf[0].times("matplotlib"), istf[0].data - stf[0].data, '-', c='b', lw=1.0, label='STF')
plt.xlabel('Time')
plt.ylabel('A')

ax.xaxis_date()
ax.xaxis.set_major_formatter(
    mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
ax.tick_params(labelleft=True, left=False)
ax.spines.right.set_visible(False)
# ax.spines.left.set_visible(False)
ax.spines.top.set_visible(False)


plt.savefig('STF.png', dpi=300)
