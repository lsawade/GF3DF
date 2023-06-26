# %%
from subprocess import check_call
from matplotlib.pyplot import *
from numpy import *
ion()

# %% Generate data
Nd = 15

check_call(f'./build/bin/test-interp {Nd} 1 > interp1d.dat', shell=True)
check_call(f'./build/bin/test-interp {Nd} 2 > spline1d.dat', shell=True)

interp1d_dat = loadtxt('interp1d.dat')
spline1d_dat = loadtxt('spline1d.dat')


# %%

xtrue = arange(0.0, 10.01, 0.01)
x = arange(0, Nd) * 10.0 / Nd

figure()
plot(xtrue, sin(xtrue)**2, c='grey', label='True')
plot(x, sin(x)**2, label='Data')
plot(interp1d_dat[:, 0], interp1d_dat[:, 1], 'o', label='LinearInterp')
plot(spline1d_dat[:, 0], spline1d_dat[:, 1], 'x', label='SplineInterp')
xlabel('x')
xlabel('y')
title('Interpolation test')
legend(loc='upper left', frameon=False, ncol=4)
ylim(-0.1, 1.2)
xlim(0, max(x))
# %%
