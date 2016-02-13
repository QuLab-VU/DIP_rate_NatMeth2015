from two_state import model
from pysb.integrate import odesolve
from pysb.bng import run_ssa
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as ss
import ratio_norm

#####################################

def dirac(x,y,a):
    return 1./(a*np.sqrt(np.pi))*np.exp(-(x-y)**2/a**2)

# Z = X*Y
def prod_pdf(xvals, px, yvals, py, zvals, alpha):
    pz = []
    for z in zvals:
        print z
        pz.append(0)
        for i in range(len(xvals)):
            for j in range(len(yvals)):
                pz[-1] += px[i]*(xvals[1]-xvals[0]) * py[j]*(yvals[1]-yvals[0]) * dirac(z, xvals[i]*yvals[j], alpha)
    return pz

# Z = X/Y
def ratio_pdf(xvals, px, yvals, py, zvals, alpha):
    pz = []
    for z in zvals:
        print z
        pz.append(0)
        for i in range(len(xvals)):
            for j in range(len(yvals)):
                pz[-1] += px[i]*(xvals[1]-xvals[0]) * py[j]*(yvals[1]-yvals[0]) * dirac(z, xvals[i]/yvals[j], alpha)
    return pz

#####################################

fig, axs = plt.subplots(2,3,figsize=(12,8))

tspan = np.arange(72.1)
nsamples = 10000

# undrugged prolif rate
k0 = model.parameters['kdiv0'].value - model.parameters['kdth0'].value

# drugged prolif rate
drug = 0.003 # uM
model.parameters['Drug0'].value = drug # uM
x = odesolve(model, tspan, verbose=False)
slope, intercept, r_value, p_value, std_err = ss.linregress(tspan[10:], np.log(x['Cell_total'][10:]))
kdrug = slope

# CONTROL (DRUG=0)
# NORMALLY DISTRIBUTED INITIAL CELL NUMBERS
# NO SAMPLING TIME VARIATION
ax = axs[0][0]
model.parameters['Drug0'].value = 0.
sample_size = 10
ncells = np.random.normal(loc=1, scale=0.1, size=nsamples*sample_size)
final_control = []
for i,n in enumerate(ncells):
    print i
    model.parameters['Cell0_init'].value = n
    if i < 100:
        plt.figure('tc')
        x = odesolve(model, tspan, verbose=False)
        plt.plot(tspan, x['Cell_total'], color='0.5')
    else:
        x = odesolve(model, [0, tspan[-1]], verbose=False)
    final_control.append(x['Cell_total'][-1])
plt.xlabel('time (h)', fontsize=18)
plt.ylabel('N (count)', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
handles = [plt.Line2D((0,1),(0,0), color='0.5', lw=5),
           plt.Line2D((0,1),(0,0), color='k', lw=5)]
labels = ['Untreated', '%d nM' % int(drug*1000.)]
plt.legend(handles, labels, loc='upper left', fontsize=18)
#####
min = 0
max = sample_size
final_mean = []
for i in range(len(ncells)/sample_size):
    final_mean.append(np.mean(final_control[min:max]))
    min = max
    max = min + sample_size
ax.hist(final_mean, bins=30, range=(14,30), normed=True)
ax.annotate("mean: %6.3f" % np.mean(final_mean), (0.55, 0.90), xycoords='axes fraction', fontsize=14)
ax.annotate("sdev: %6.3f" % np.std(final_mean), (0.55, 0.83), xycoords='axes fraction', fontsize=14)
ax.annotate("skew: %6.3f" % ss.skew(final_mean), (0.55, 0.76), xycoords='axes fraction', fontsize=14)
#####
xvals = np.arange(16,24,0.1)
ax.plot(xvals, ss.norm.pdf(xvals, np.exp(k0*72.), 0.1*np.exp(k0*72.)/np.sqrt(sample_size)), color='r', lw=3)
#####
ax.set_xlim(14,30)
# ax.set_xlabel(r"$N_\mathrm{untreated}$ (72 h)")
ax.set_ylabel('density', fontsize=14)
ax.tick_params(labelsize=14)

# DRUG TREATED
# NORMALLY DISTRIBUTED INITIAL CELL NUMBERS
# NO SAMPLING TIME VARIATION
ax = axs[0][1]
model.parameters['Drug0'].value = drug # uM
ncells = np.random.normal(loc=1, scale=0.1, size=nsamples)
final_treated = []
for i,n in enumerate(ncells):
    print i
    model.parameters['Cell0_init'].value = n
    if i < 1000:
        plt.figure('tc')
        x = odesolve(model, tspan, verbose=False)
        plt.plot(tspan, x['Cell_total'], color='k')
    else:
        x = odesolve(model, [0, tspan[-1]], verbose=False)
    final_treated.append(x['Cell_total'][-1])
plt.savefig('timecourses.pdf', format='pdf', dpi=1000)
#####
ax.hist(final_treated, bins=30, normed=True, range=(4,12))
ax.annotate("mean: %6.3f" % np.mean(final_treated), (0.55, 0.90), xycoords='axes fraction', fontsize=14)
ax.annotate("sdev: %6.3f" % np.std(final_treated), (0.55, 0.83), xycoords='axes fraction', fontsize=14)
ax.annotate("skew: %6.3f" % ss.skew(final_treated), (0.55, 0.76), xycoords='axes fraction', fontsize=14)
#####
xvals = np.arange(4,12,0.1)
ax.plot(xvals, ss.norm.pdf(xvals, np.exp(kdrug*72.), 0.1*np.exp(kdrug*72.)), color='r', lw=3)
#####
ax.set_xlim(4,12)
# ax.set_xlabel(r"$N_\mathrm{12\ nM}$ (72 h)")
# ax.set_ylabel('density')
ax.tick_params(labelsize=14)

# RESPONSE RATIO (DRUG TREATED / <CONTROL>)
# NORMALLY DISTRIBUTED INITIAL CELL NUMBERS
# NO SAMPLING TIME VARIATION
ax = axs[0][2]
final_ratio = [x/y for x in final_treated for y in final_mean]
# final_ratio = [y/x for x in final_treated for y in final_mean]
ax.hist(final_ratio, bins=30, normed=True, range=(0.1,0.8))
ax.annotate("mean: %6.3f" % np.mean(final_ratio), (0.55, 0.90), xycoords='axes fraction', fontsize=14)
ax.annotate("sdev: %6.3f" % np.std(final_ratio), (0.55, 0.83), xycoords='axes fraction', fontsize=14)
ax.annotate("skew: %6.3f" % ss.skew(final_ratio), (0.55, 0.76), xycoords='axes fraction', fontsize=14)
#####
xvals = np.arange(0.1,0.8,0.001)
mean_drug = np.exp(kdrug*72.)
sdev_drug = 0.1*np.exp(kdrug*72.)
mean_ctrl = np.exp(k0*72.)
sdev_ctrl = 0.1*np.exp(k0*72.)/np.sqrt(sample_size)
ax.plot(xvals, ratio_norm.pdf(xvals, (mean_drug, sdev_drug), (mean_ctrl, sdev_ctrl)), color='r', lw=3)
# ax.plot(xvals, ratio_pdf(np.arange(1.,3.01,0.01), ss.norm.pdf(np.arange(1.,3.01,0.01), np.exp(kdrug*72.), 0.1*np.exp(kdrug*72.)), 
#                          np.arange(16, 24+0.1, 0.1), ss.norm.pdf(np.arange(16, 24+0.1, 0.1), np.exp(k0*72.), 0.1*np.exp(k0*72.)/np.sqrt(sample_size)), 
#                          xvals, 0.001), color='g', lw=6, ls='--')
#####
ax.set_xlim(0.1,0.8)
# ax.set_xlabel(r'$N_\mathrm{12\ nM} / N_\mathrm{untreated}$ (72 h)')
# ax.set_ylabel('density')
ax.tick_params(labelsize=14)

# plt.show()
# quit()

#####################################

# CONTROL (DRUG=0)
# NORMALLY DISTRIBUTED INITIAL CELL NUMBERS
# NORMALLY DISTRIBUTED SAMPLING TIMES
ax = axs[1][0]
model.parameters['Drug0'].value = 0.
ncells = np.random.normal(loc=1, scale=0.1, size=nsamples*sample_size)
times = np.random.normal(loc=72, scale=5, size=len(ncells))
final_control = []
for i,n in enumerate(ncells):
    print i
    model.parameters['Cell0_init'].value = n
    x = odesolve(model, [0, times[i]], verbose=False)
    final_control.append(x['Cell_total'][-1])
#####
min = 0
max = sample_size
final_mean = []
for i in range(len(ncells)/sample_size):
    final_mean.append(np.mean(final_control[min:max]))
    min = max
    max = min + sample_size
ax.hist(final_mean, bins=30, range=(14,30), normed=True)
ax.annotate("mean: %6.3f" % np.mean(final_mean), (0.55, 0.90), xycoords='axes fraction', fontsize=14)
ax.annotate("sdev: %6.3f" % np.std(final_mean), (0.55, 0.83), xycoords='axes fraction', fontsize=14)
ax.annotate("skew: %6.3f" % ss.skew(final_mean), (0.55, 0.76), xycoords='axes fraction', fontsize=14)
#####
# xvals = np.arange(0.1,2.,0.01)
# px = ss.norm.pdf(xvals, 1, 0.1)
# k = model.parameters['kdiv0'].value - model.parameters['kdth0'].value
# yvals = np.arange(1.,100.,0.1)
# py = ss.lognorm.pdf(yvals, k0*5., scale=np.exp(k0*72.))
# zvals_ctrl = np.arange(1.,50.1,1.)
# pz_ctrl = prod_pdf(xvals, px, yvals, py, zvals_ctrl, 0.25)
# ax.plot(zvals_ctrl, pz_ctrl, color='r', lw=3)
#####
ax.set_xlim(14,30)
ax.set_xlabel(r"$<N_\mathrm{untreated}>$ (72 h)", fontsize=14)
ax.set_ylabel('density', fontsize=14)
ax.tick_params(labelsize=14)

# DRUG TREATED
# NORMALLY DISTRIBUTED INITIAL CELL NUMBERS
# NORMALLY DISTRIBUTED SAMPLING TIMES
ax = axs[1][1]
model.parameters['Drug0'].value = drug # uM
ncells = np.random.normal(loc=1, scale=0.1, size=nsamples)
times = np.random.normal(loc=72, scale=5, size=len(ncells))
final_treated = []
for i,n in enumerate(ncells):
    print i
    model.parameters['Cell0_init'].value = n
    x = odesolve(model, [0, times[i]], verbose=False)
    final_treated.append(x['Cell_total'][-1])
#####
ax.hist(final_treated, bins=30, normed=True, range=(4,12))
ax.annotate("mean: %6.3f" % np.mean(final_treated), (0.55, 0.90), xycoords='axes fraction', fontsize=14)
ax.annotate("sdev: %6.3f" % np.std(final_treated), (0.55, 0.83), xycoords='axes fraction', fontsize=14)
ax.annotate("skew: %6.3f" % ss.skew(final_treated), (0.55, 0.76), xycoords='axes fraction', fontsize=14)
#####
# xvals = np.arange(0.1,2.,0.01)
# px = ss.norm.pdf(xvals, 1, 0.1)
# # koff_over_kon = model.parameters['koff'].value / model.parameters['kon'].value
# # k0 = model.parameters['kdiv0'].value - model.parameters['kdth0'].value
# # kmax = model.parameters['kdiv1'].value - model.parameters['kdth1'].value
# # k = (koff_over_kon*k0 + model.parameters['Drug0'].value*kmax) / (koff_over_kon + model.parameters['Drug0'].value)
# yvals = np.arange(1.,3.,0.01)
# py = ss.lognorm.pdf(yvals, kdrug*5., scale=np.exp(kdrug*72.))
# zvals_drug = np.arange(1.,3.,0.05)
# pz_drug = prod_pdf(xvals, px, yvals, py, zvals_drug, 0.025)
# ax.plot(zvals_drug, pz_drug, color='r', lw=3)
#####
ax.set_xlim(4,12)
ax.set_xlabel(r"$N_\mathrm{12\ nM}$ (72 h)", fontsize=14)
# ax.set_ylabel('density')
ax.tick_params(labelsize=14)

# RESPONSE RATIO (DRUG TREATED / <CONTROL>)
# NORMALLY DISTRIBUTED INITIAL CELL NUMBERS
# NORMALLY DISTRIBUTED SAMPLING TIMES
ax = axs[1][2]
# plt.figure()
final_ratio2 = [x/y for x in final_treated for y in final_mean]
ax.hist(final_ratio2, bins=30, normed=True, range=(0.01,0.8))
ax.annotate("mean: %6.3f" % np.mean(final_ratio2), (0.55, 0.90), xycoords='axes fraction', fontsize=14)
ax.annotate("sdev: %6.3f" % np.std(final_ratio2), (0.55, 0.83), xycoords='axes fraction', fontsize=14)
ax.annotate("skew: %6.3f" % ss.skew(final_ratio2), (0.55, 0.76), xycoords='axes fraction', fontsize=14)
#####
# zvals_ratio = np.arange(0.01,0.3,0.005)
# pz_ratio = ratio_pdf(zvals_drug, pz_drug, zvals_ctrl, pz_ctrl, zvals_ratio, 0.01)
# ax.plot(zvals_ratio, pz_ratio, color='r', lw=3)
#####
ax.set_xlabel(r'$N_\mathrm{12\ nM} / <N_\mathrm{untreated}>$ (72 h)', fontsize=14)
# ax.set_ylabel('density')
ax.set_xlim(0.1,0.8)
ax.tick_params(labelsize=14)

fig.tight_layout()
fig.savefig('distributions.pdf', format='pdf', dpi=1000)

#####
plt.figure()
plt.boxplot([final_ratio, final_ratio2], showfliers=False, showmeans=True, boxprops={'lw':2}, medianprops={'lw':2}, 
            whiskerprops={'lw':2}, capprops={'lw':2}, meanprops={'ms':8, 'markerfacecolor':'k', 'marker':'o'})
#x, notch, sym, vert, whis, positions, widths, patch_artist, bootstrap, usermedians, conf_intervals, meanline, showmeans, showcaps, 
#showbox, showfliers, boxprops, labels, flierprops, medianprops, meanprops, capprops, whiskerprops, manage_xticks, hold)
plt.ylabel('static response ratio (72 h)', fontsize=18)
plt.xticks([1, 2], ['seeding density', 'seeding density\n+ measurement time'], fontsize=18)
plt.yticks(fontsize=18)
plt.savefig('boxplots.pdf', format='pdf', dpi=1000)
#####

plt.show()
