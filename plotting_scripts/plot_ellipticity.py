# %%
import numpy as np
from numpy import *
from matplotlib.pyplot import *
from scipy.integrate import cumtrapz, simpson
from subprocess import check_call

check_call(f'./build/bin/make-ellipticity > test_ellipticity.dat', shell=True)

ellipticity_dat = loadtxt('test_ellipticity.dat')

# %%
ONE_CRUST = True
EARTH_RHOAV = 5514.3
EARTH_R = 6371000.0

# Number of sample points
NR_DENSITY = 640

PREM_RSURFACE = 6371000.0
PREM_ROCEAN = 6368000.0
PREM_RMIDDLE_CRUST = 6356000.0  # 15 km depth
PREM_RMOHO = PREM_RSURFACE - 24400.0  # 6346600.  at 24.4km depth
PREM_R80 = 6291000.0  # 80 km
PREM_R220 = 6151000.0
PREM_R400 = 5971000.0
PREM_R600 = 5771000.0
PREM_R670 = 5701000.0                 # at 670 km depth
PREM_R771 = 5600000.0
PREM_RTOPDDOUBLEPRIME = 3630000.0
PREM_RCMB = 3480000.0  # 2891 km depth
# note: SPECFEM versions up to 8.0 (Aug, 2021) used an inner core radius of 1221 km
#       based on Table 1 in Dziewonski & Anderson's PREM paper, the inner core radius is at 1221.5 km
#       double precision, parameter: : PREM_RICB = 1221000.            # old versions
PREM_RICB = 1221500.0

# PREM2 additional radii for modifications
PREM2_RDDOUBLEPRIME_UPPER = 3840000.0  # upper D'' region at 3840km radius
PREM2_ROC_LOWER = 1621500.0            # lower outer core at 1621.5km radius
PREM2_RIC_UPPER = 1010000.0            # upper inner core at 1010 km radius

# %%
R_EARTH_ELLIPTICITY = 6371000.0
ROCEAN_ELLIPTICITY = PREM_ROCEAN

RSURFACE = R_EARTH_ELLIPTICITY       # physical surface(Earth: 6371000, ..)
ROCEAN = ROCEAN_ELLIPTICITY
RMIDDLE_CRUST = PREM_RMIDDLE_CRUST
RMOHO = PREM_RMOHO
R80 = PREM_R80
R220 = PREM_R220
R400 = PREM_R400
R600 = PREM_R600
R670 = PREM_R670
R771 = PREM_R771
RTOPDDOUBLEPRIME = PREM_RTOPDDOUBLEPRIME
RCMB = PREM_RCMB
RICB = PREM_RICB

# %%

# non-dimensionalize
r_icb = RICB / RSURFACE
r_cmb = RCMB / RSURFACE
r_topddoubleprime = RTOPDDOUBLEPRIME / RSURFACE
r_771 = R771 / RSURFACE
r_670 = R670 / RSURFACE
r_600 = R600 / RSURFACE
r_400 = R400 / RSURFACE
r_220 = R220 / RSURFACE
r_80 = R80 / RSURFACE
r_moho = RMOHO / RSURFACE
r_middle_crust = RMIDDLE_CRUST / RSURFACE
r_ocean = ROCEAN / RSURFACE
r_0 = 1.0

# %%


def xrange(start, stop, step):
    return np.arange(start, stop+step, step)


r = zeros(NR_DENSITY)
rho = zeros(NR_DENSITY)
radau = zeros(NR_DENSITY)


for i in range(1, 164):
    r[i-1] = r_icb*(i-1)/(162)

# outer core
for i in range(164, 324):
    r[i-1] = r_icb+(r_cmb-r_icb)*(i-164)/(159)

# D''
for i in range(324, 337):
    r[i-1] = r_cmb+(r_topddoubleprime-r_cmb)*(i-324)/(12)

# D'' to 771
for i in range(337, 518):
    r[i-1] = r_topddoubleprime+(r_771-r_topddoubleprime)*(i-337)/(180)

# 771 to 670
for i in range(518, 531):
    r[i-1] = r_771+(r_670-r_771)*(i-518)/(12)

# 670 to 600
for i in range(531, 541):
    r[i-1] = r_670+(r_600-r_670)*(i-531)/(9)

# 600 to 400
for i in range(541, 566):
    r[i-1] = r_600+(r_400-r_600)*(i-541)/(24)

# 400 to 220
for i in range(566, 591):
    r[i-1] = r_400+(r_220-r_400)*(i-566)/(24)

  # 220 to 80
for i in range(591, 610):
    r[i-1] = r_220+(r_80-r_220)*(i-591)/(18)

# 80 to Moho
for i in range(610, 620):
    r[i-1] = r_80+(r_moho-r_80)*(i-610)/(9)

# Moho to middle crust
for i in range(620, 627):
    r[i-1] = r_moho+(r_middle_crust-r_moho)*(i-620)/(6)

# middle crust to ocean
for i in range(627, 634):
    r[i-1] = r_middle_crust+(r_ocean-r_middle_crust)*(i-627)/(6)

# ocean
for i in range(634, NR_DENSITY+1):   # NR_DENSITY = 640):
    r[i-1] = r_ocean+(r_0-r_ocean)*(i-634)/(6)

# %%


def prem_density(x):

    ONE_CRUST = True
    EARTH_RHOAV = 5514.3
    EARTH_R = 6371000.0

    R_PLANET = EARTH_R
    RHOAV = EARTH_RHOAV

    # compute real physical radius in meters
    r = x * R_PLANET

    # calculates density according to radius
    if (r <= PREM_RICB):
        rho = 13.0885 - 8.8381*x*x
    elif (r > PREM_RICB and r <= PREM_RCMB):
        rho = 12.5815 - 1.2638*x - 3.6426*x*x - 5.5281*x*x*x
    elif (r > PREM_RCMB and r <= PREM_RTOPDDOUBLEPRIME):
        rho = 7.9565 - 6.4761*x + 5.5283*x*x - 3.0807*x*x*x
    elif (r > PREM_RTOPDDOUBLEPRIME and r <= PREM_R771):
        rho = 7.9565 - 6.4761*x + 5.5283*x*x - 3.0807*x*x*x
    elif (r > PREM_R771 and r <= PREM_R670):
        rho = 7.9565 - 6.4761*x + 5.5283*x*x - 3.0807*x*x*x
    elif (r > PREM_R670 and r <= PREM_R600):
        rho = 5.3197 - 1.4836*x
    elif (r > PREM_R600 and r <= PREM_R400):
        rho = 11.2494 - 8.0298*x
    elif (r > PREM_R400 and r <= PREM_R220):
        rho = 7.1089 - 3.8045*x
    elif (r > PREM_R220 and r <= PREM_R80):
        rho = 2.6910 + 0.6924*x
    else:
        if (r > PREM_R80 and r <= PREM_RMOHO):
            rho = 2.6910 + 0.6924*x
        elif (r > PREM_RMOHO and r <= PREM_RMIDDLE_CRUST):
            if (ONE_CRUST):
                rho = 2.6  # takes upper crust value
            else:
                rho = 2.9
        elif (r > PREM_RMIDDLE_CRUST and r <= PREM_ROCEAN):
            rho = 2.6
        elif (r > PREM_ROCEAN):
            rho = 2.6  # extends upper crust

    # non-dimensionalizes
    rho = rho * 1000.0 / RHOAV

    return rho


# %%
for i in range(0, NR_DENSITY):
    rho[i] = prem_density(r[i])

rhorr = rho*r*r
rhorrrr = rhorr*r*r

eta = zeros(NR_DENSITY)
k = zeros(NR_DENSITY)

for i in range(1, NR_DENSITY):

    # Integral
    int_rho = simpson(rhorr[:i+1], r[:i+1])  # cumtrapz(rhorr,r)
    int_radau = simpson(rhorrrr[:i+1], r[:i+1])

    print(f"{int_radau:g}  {int_rho:g}")
    # Intermediate val
    z = (2.0/3.0) * int_radau / (int_rho * r[i] * r[i])
    # print(r[i], z)
    # This comes from equation (14.19) in Dahlen and Tromp (1998)

    # Compute eta
    eta[i] = (25.0/4.0) * ((1.0 - (3.0/2.0) * z)**2.0)-1.0
    # print(r[i], eta[i])
    k[i] = eta[i]/(r[i])

# %%
# day rotation
# bom = TWO_PI/(HOURS_PER_DAY*SECONDS_PER_HOUR)
bom = 2*pi/(24*3600)

EARTH_RHOAV = 5514.3
RHOAV = EARTH_RHOAV
GRAV = 6.67384e-11

# non-dimensionalized value
bom = bom/sqrt(pi*GRAV*RHOAV)

g_a = 4.0 * int_rho

R_UNIT_SPHERE = 1.0

# Blub
epsilonval = zeros(NR_DENSITY)

# this is the equation right above(14.21) in Dahlen and Tromp(1998)
epsilonval[NR_DENSITY-1] = (5.0/2.0)*(bom**2.0) * \
    R_UNIT_SPHERE / (g_a * (eta[NR_DENSITY-1]+2.0))

exponentval = zeros(NR_DENSITY)

krr = k*r*r

for i in range(0, NR_DENSITY):
    expval = trapz(k[i+1:NR_DENSITY+1], r[i+1:NR_DENSITY+1])
    print(expval)
    epsilonval[i] = epsilonval[NR_DENSITY-1] * \
        exp(-expval)  # * exp(-exponentval[-(i+1)])

# %%


def plot_boundaries(ax, yoffset=0.79, label=False, bottom=False, factor=1.75):
    ax.axvline(6371-R670/1000.0, label='670', c='k', lw=0.5)
    ax.axvline(6371-RCMB/1000.0, label='CMB', c='k', lw=0.5)
    ax.axvline(6371-RICB/1000.0, label='ICB', c='k', lw=0.5)

    if label:
        ax.text(6371-R670/1000.0+10, yoffset, '670',
                rotation=90, ha='left', va='top')
        ax.text(6371-RCMB/1000.0+10, yoffset, 'CMB',
                rotation=90, ha='left', va='top')
        ax.text(6371-RICB/1000.0+10, yoffset, 'ICB',
                rotation=90, ha='left', va='top')

    if bottom:
        ax.axvline(6371-R670/1000.0, ymax=yoffset *
                   factor, clip_on=False, c='k', lw=0.5)
        ax.axvline(6371-RCMB/1000.0, ymax=yoffset *
                   factor, clip_on=False, c='k', lw=0.5)
        ax.axvline(6371-RICB/1000.0, ymax=yoffset *
                   factor, clip_on=False, c='k', lw=0.5)


figure(figsize=(8, 4))

ax1 = subplot(2, 1, 1)
# plot(6371-r*6371, epsilonval)
plot(6371-ellipticity_dat[:, 0], ellipticity_dat[:, 1], 'k')
plot_boundaries(ax1, yoffset=0.00395, label=True)
xlim(0, 6371)
ylim(0.002, 0.004)
ylabel('$\epsilon$')
ax1.tick_params(axis='x', labelbottom=False)

ax2 = subplot(2, 1, 2)
# plot(6371-r*6371, eta)
plot(6371-ellipticity_dat[:, 0], ellipticity_dat[:, 2], 'k')
plot_boundaries(ax2, yoffset=0.79, label=False, bottom=True)

ylim(0, 0.8)
xlim(0, 6371)
xlabel('Depth [km]')
ylabel('$\eta$')
subplots_adjust(left=0.1, right=0.9, hspace=0.2)
savefig('test_ellipticity.pdf')


# %%
ion()
figure(figsize=(5, 4))

plot(r, arange(NR_DENSITY)/(NR_DENSITY-1), label='Norm. Index')
plot(r, rho, label='Non-dim. $\\rho(i)$')
plot(r, eta, label='Radau$(i)$')

xlabel('Radius $r[i-1]$/$R_{E}$')
ylabel('Normalized values')

legend(frameon=False, ncol=1, loc='center left')
