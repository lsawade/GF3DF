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

figure(figsize=(8, 2.5))
ax = gca()
plot(xtrue, sin(xtrue)**2, c='grey', label='True')
plot(x, sin(x)**2, '-', label='Data')

# Lines
iline, = plot(interp1d_dat[:, 0], interp1d_dat[:, 1], 'o',
              label='LinearInterp', fillstyle='none')
sline, = plot(spline1d_dat[:, 0], spline1d_dat[:, 1],
              'x', label='SplineInterp')

# Errors
plot(interp1d_dat[:, 0], abs(sin(interp1d_dat[:, 0])**2-interp1d_dat[:, 1]),
     c=iline.get_color(), label='|LinearErr|', fillstyle='none')
plot(spline1d_dat[:, 0], abs(sin(spline1d_dat[:, 0])**2-spline1d_dat[:, 1]),
     c=sline.get_color(), label='|SplineErr|', fillstyle='none')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

xlabel('x')
ylabel('y')
# title('Interpolation test')
legend(loc='upper center', frameon=False, ncol=3)
ylim(-0.1, 1.5)
xlim(0, 10)
subplots_adjust(left=0.075, right=0.925, bottom=0.175, top=0.95)
savefig('test_interp.pdf')

# %%
