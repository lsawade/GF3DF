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

cmt = CMTSOLUTION.read('CMTSOLUTION')

# %% Initialize the Green Function Manager

gfm = GFManager('single_element_not_fortran.h5')
# gfm = GFManager('/scratch/gpfs/lsawade/large_subset.h5')

gfm.load()

# %% Get Python seismograms

st = gfm.get_seismograms(cmt)

pbfopy = st.select(station='BFO')

# bfopy.differentiate()

# pbfopy = process_stream(bfopy, cmt=cmt)

# %% Read seismograms from fortran
pbfof = read('OUTPUT/II.BFO.*.sem.sac')

# factor = 1
# for tr in bfof:

#     tr.data *= factor

if all([os.path.exists(f'OUTPUT_WIN/II.BFO.MX{C}.sem.sac') for C in ['N', 'E', 'Z']]):
    pbfof_win = read('OUTPUT_WIN/II.BFO.*.sac')
    # pbfof_win = process_stream(bfof_win, cmt=cmt)
else:
    pbfof_win = None

# pbfof = process_stream(bfof, cmt=cmt)
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
