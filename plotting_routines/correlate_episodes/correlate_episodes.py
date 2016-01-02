import numpy as np
import matplotlib.pyplot as plt
from datetime import date
import sys
today = str(date.today())
print today

stars_formed_thresh =0.01


if (len(sys.argv) < 2):
	print 'syntax  blah.py outflow_ep_file  vesc_ep_file'
	sys.exit()

finname1 = str(sys.argv[1])
finname2 = str(sys.argv[2])


f = open(finname1)
dars = np.loadtxt(f)
f.close()

f2 = open(finname2)
dars2  = np.loadtxt(f2)
f2.close()



Mstar = dars[:,8]
eta = dars[:,17]
formed_fraction = dars[:,19]
cut = formed_fraction>stars_formed_thresh
escape_v = dars2[:,3]
surf = dars2[:,7]
Mgained = dars[:,22]
tInflow_start = dars[:,23]
tSF_end =  dars[:,15]
duration_of_inflow =  tSF_end - tInflow_start
Inflow_rate = Mgained*-1e0 / duration_of_inflow



figure0 = plt.figure(figsize=(10, 8))
plt.xlabel('Mstar ', fontsize=20)
plt.ylabel('eta ', fontsize=20)
plt.xscale('log')
plt.yscale('log')
ax = plt.gca()
ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(18)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(18)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=2.5, length=8)
ax.yaxis.set_tick_params(width=2.5, length=8)

plt.plot(Mstar[cut], eta[cut], 'ok')
foutname=today+'Mstar-cor.pdf'
plt.savefig(foutname)
plt.clf()

figure1 = plt.figure(figsize=(10, 8))
plt.xlabel('v_esc ', fontsize=20)
plt.ylabel('eta ', fontsize=20)
plt.xscale('log')
plt.yscale('log')
ax = plt.gca()
ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(18)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(18)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=2.5, length=8)
ax.yaxis.set_tick_params(width=2.5, length=8)

plt.plot(escape_v[cut], eta[cut], 'ok')
foutname=today+'vesc-cor.pdf'
plt.savefig(foutname)
plt.clf()

figure2 = plt.figure(figsize=(10, 8))
plt.xlabel('surface density ', fontsize=20)
plt.ylabel(r'$\eta$', fontsize=20)
plt.xscale('log')
plt.yscale('log')
ax = plt.gca()
ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(18)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(18)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=2.5, length=8)
ax.yaxis.set_tick_params(width=2.5, length=8)

plt.plot(surf[cut], eta[cut], 'ok')
foutname=today+'surf.pdf'
plt.savefig(foutname)

weirdrat = surf/(escape_v**1.0)

figure3 = plt.figure(figsize=(10, 8))
plt.xlabel('ratio of surface to escape velocity ', fontsize=20)
plt.ylabel(r'$\eta$', fontsize=20)
plt.xscale('log')
plt.yscale('log')
ax = plt.gca()
ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(18)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(18)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=2.5, length=8)
ax.yaxis.set_tick_params(width=2.5, length=8)

plt.plot(weirdrat[cut], eta[cut], 'ok')
foutname=today+'weirdrat.pdf'
plt.savefig(foutname)

figure4 = plt.figure(figsize=(10, 8))
plt.xlabel('Inflow rate (Msun/yr)', fontsize=20)
plt.ylabel(r'$\eta$', fontsize=20)
plt.xscale('log')
plt.yscale('log')
ax = plt.gca()
ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(18)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(18)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=2.5, length=8)
ax.yaxis.set_tick_params(width=2.5, length=8)

plt.plot(Inflow_rate[cut], eta[cut], 'ok')
foutname=today+'inflow.pdf'
plt.savefig(foutname)
