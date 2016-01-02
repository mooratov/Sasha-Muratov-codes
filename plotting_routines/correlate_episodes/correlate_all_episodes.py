import numpy as np
import matplotlib.pyplot as plt
from datetime import date
import sys
today = str(date.today())
print today
import scipy.stats as stats
import scipy.optimize as optimize
import sasha_stats as SaSta

stars_formed_thresh =0.05
avgSFRthresh = 2
themarker = 10
acceptablefrac = 0.0
etathresh = 0.01
dobadplots = True

if (len(sys.argv) < 2):
	print 'syntax  blah.py list_outflow_ep_file  list_vesc_ep_file list_outflow_ep-inflowsig'
	sys.exit()

finname1 = str(sys.argv[1])
finname2 = str(sys.argv[2])
finname3 = str(sys.argv[3])


finname_ar = []
f = open(finname1, 'r')
dars = f.readlines()
for line in dars:
	xsd = line.split()
	finname_ar.append(str(xsd[0]))
f.close()	

finnameB_ar = []
f = open(finname2, 'r')
dars = f.readlines()
for line in dars:
	xsd = line.split()
	finnameB_ar.append(str(xsd[0]))
f.close()	

finnameC_ar = []
f = open(finname3, 'r')
dars = f.readlines()
for line in dars:
	xsd = line.split()
	finnameC_ar.append(str(xsd[0]))
f.close()	

Mstar_ar = []
eta_ar = []
formed_frac_ar = []
cut_ar = []
escape_v_ar = []
surf_ar = []
Mgained_ar = []
tInflow_start_ar = []
tSF_end_ar = []
duration_of_inflow_ar = []
Inflow_rate_ar  = []
SFR_ar = []

Mstar_ar_bad = []
eta_ar_bad = []
escape_v_ar_bad = []
surf_ar_bad = []
Inflow_rate_ar_bad = []
SFR_ar_bad = []

vesc_ext_ar = []
eta_ext_ar = []
surf_ext_ar = []
Inflow_ext_ar = []
SFR_ext_ar = []

thar = ['B1 med-z', 'B1 lo-z', 'm12i med-z', 'm12i lo-z', 'm12q med-z', 'm12q lo-z' ,  'CAFG 506', 'CAFG 600', 'Feld 37', 'Feld 206']

count = 0
while (count < len(finname_ar)):
	f = open(finname_ar[count])
	dars = np.loadtxt(f, ndmin=2)
	f.close()
	
	f2 = open(finnameB_ar[count])
	dars2  = np.loadtxt(f2, ndmin=2)
	f2.close()
	
	f3 = open(finnameC_ar[count])
	dars3  = np.loadtxt(f3, ndmin=2)
	f3.close()

	N = dars[:,0]
	Z = dars[:,1]
	Mstar = dars[:,8]
	eta = dars[:,17]
	etacut = eta > etathresh
	formed_fraction = dars[:,19]
	cut = formed_fraction>stars_formed_thresh
	#dars2:  3 - mean escape v. 4- 25th percentile. 5 - 50th percentile. 6 - 75th percentile
	escape_v = dars2[:,3]
	#dras2: 7 - mean weighted surface density, 8 - mean unweighted suface density
	surf = dars2[:,7]
	Mgained = dars[:,22]
	tInflow_start = dars[:,23]
	tSF_start = dars[:,13]
	Mformed = dars[:,3]
	tSF_end =  dars[:,15]
	
	maxinflow = -1*dars3[:,24]
	maxoutflow = dars3[:,12]
	
	outcut = maxoutflow > acceptablefrac*maxinflow
	
	duration_of_inflow =  tSF_end - tInflow_start
	duration_of_SF = tSF_end - tSF_start 
	avgSFR = Mformed / duration_of_SF
	#print avgSFR
	
	Inflow_rate = (Mgained*-1e0 / duration_of_inflow)
	SFRcut = avgSFR > avgSFRthresh
	cut*=SFRcut * outcut * etacut

	Mstar_ar.append(Mstar[cut])
	eta_ar.append(eta[cut])
	escape_v_ar.append(escape_v[cut])
	surf_ar.append(surf[cut])
	Inflow_rate_ar.append(Inflow_rate[cut])
	SFR_ar.append(avgSFR[cut])
	
	anticut = np.invert(cut)
	
	Mstar_ar_bad.append(Mstar[anticut])
	eta_ar_bad.append(eta[anticut])
	escape_v_ar_bad.append(escape_v[anticut])
	surf_ar_bad.append(surf[anticut])
	Inflow_rate_ar_bad.append(Inflow_rate[anticut])
	SFR_ar_bad.append(avgSFR[anticut])
	
	vesc_ext_ar.extend(escape_v[cut])
	eta_ext_ar.extend(eta[cut])
	surf_ext_ar.extend(surf[cut])
	Inflow_ext_ar.extend(Inflow_rate[cut])
	SFR_ext_ar.extend(avgSFR[cut])
	print finname_ar[count], N[cut], Z[cut]
	count +=1

Mstar_ar = np.array(Mstar_ar)
eta_ar = np.array(eta_ar)
escape_v_ar = np.array(escape_v_ar)
surf_ar = np.array(surf_ar)
Inflow_rate_ar = np.array(Inflow_rate_ar)
SFR_ar = np.array(SFR_ar)

Mstar_ar_bad = np.array(Mstar_ar_bad)
eta_ar_bad = np.array(eta_ar_bad)
escape_v_ar_bad = np.array(escape_v_ar_bad)
surf_ar_bad = np.array(surf_ar_bad)
Inflow_rate_ar_bad = np.array(Inflow_rate_ar_bad)
SFR_ar_bad = np.array(SFR_ar_bad)



vesc_ext_ar = np.array(vesc_ext_ar)
eta_ext_ar = np.array(eta_ext_ar)
surf_ext_ar = np.array(surf_ext_ar)
Inflow_ext_ar = np.array(Inflow_ext_ar)
SFR_ext_ar = np.array(SFR_ext_ar)



vesc_ext_ar_log = np.log10(vesc_ext_ar)
eta_ext_ar_log = np.log10(eta_ext_ar)
surf_ext_ar_log = np.log10(surf_ext_ar)
Inflow_ext_ar_log = np.log10(Inflow_ext_ar)
SFR_ext_ar_log = np.log10(SFR_ext_ar)
p3 = [-1, 1, 1, 1]
p1 = [-1, 1]

def residuals_skinny(b , y, myvesc):
    err = y - (b[0] + b[1]*myvesc)
    return err
    
print len(eta_ext_ar_log)
print len(vesc_ext_ar_log)
print eta_ext_ar_log
print vesc_ext_ar_log
vescfit = optimize.leastsq(residuals_skinny,p1[:], args=(eta_ext_ar_log, vesc_ext_ar_log), full_output=True)
print 'the vesc fit ',vescfit[0]

vescfitcov = vescfit[1]
perrS = np.sqrt(np.diag(vescfitcov))
print 'full perrS', perrS

surffit = optimize.leastsq(residuals_skinny,p1[:], args=(eta_ext_ar_log, surf_ext_ar_log), full_output=True)
print 'the surf fit ',surffit[0]

surffitcov = surffit[1]
perrS = np.sqrt(np.diag(surffitcov))
print 'full perrS', perrS

Inflowfit = optimize.leastsq(residuals_skinny,p1[:], args=(eta_ext_ar_log, Inflow_ext_ar_log), full_output=True)
print 'the Inflow fit ',Inflowfit[0]

inflowfitcov = Inflowfit[1]
perrS = np.sqrt(np.diag(inflowfitcov))
print 'full perrS', perrS

SFRfit = optimize.leastsq(residuals_skinny,p1[:], args=(eta_ext_ar_log, SFR_ext_ar_log), full_output=True)
print 'the SFR fit ',SFRfit[0]

SFRfitcov = SFRfit[1]
perrS = np.sqrt(np.diag(SFRfitcov))
print 'full perrS', perrS


def residuals(b , y, myvesc, mysurf, myflow):
    err = y - (b[0] + b[1]*myvesc +  b[2]*mysurf + b[3]*myflow)
    return err


megafit = optimize.leastsq(residuals,p3[:], args=(eta_ext_ar_log, vesc_ext_ar_log, surf_ext_ar_log, Inflow_ext_ar_log), full_output=True)

print 'the all triple fit ',megafit[0]

megafitcov = megafit[1]
#print 'megafit cov ',megafitcov
perrS = np.sqrt(np.diag(megafitcov))
print 'full perrS', perrS
#print 'perrS ',perrS[0]

megafit2 = optimize.leastsq(residuals,p3[:], args=(eta_ext_ar_log, vesc_ext_ar_log, surf_ext_ar_log, SFR_ext_ar_log), full_output=True)

print 'the triple fit 2(SFR)',megafit2[0]

megafit2cov = megafit2[1]
#print 'megafit cov 2',megafitcov
perrS2 = np.sqrt(np.diag(megafit2cov))
print 'full perrS 2', perrS2
#print 'perrS ',perrS[0]

p2 = [-1, 1, 1]
def residuals_just2(b , y, myvesc, mysurf):
    err = y - (b[0] + b[1]*myvesc +  b[2]*mysurf)
    return err
hugefit = optimize.leastsq(residuals_just2,p2[:], args=(eta_ext_ar_log, vesc_ext_ar_log, surf_ext_ar_log), full_output=True)

print 'huge fit vesc and sruf',hugefit[0]

hugefitcov = hugefit[1]
perrS3 = np.sqrt(np.diag(hugefitcov))
print 'full perrS 3', perrS3

hugefit2 = optimize.leastsq(residuals_just2,p2[:], args=(eta_ext_ar_log, vesc_ext_ar_log, Inflow_ext_ar_log), full_output=True)

print 'huge fit 2 vesc and inflow',hugefit2[0]

hugefitcov2 = hugefit2[1]
perrS4 = np.sqrt(np.diag(hugefitcov2))
print 'full perrS 4', perrS4


hugefit3 = optimize.leastsq(residuals_just2,p2[:], args=(eta_ext_ar_log, surf_ext_ar_log, Inflow_ext_ar_log), full_output=True)

print 'huge fit 3 surf and inflow',hugefit3[0]

hugefitcov3 = hugefit3[1]
perrS5 = np.sqrt(np.diag(hugefitcov3))
print 'full perrS 5', perrS5


color_ar = ['k', 'r', 'b', 'g', 'm', 'y', 'DarkSlateBlue', 'LightGreen', 'DarkGoldenRod', 'DeepPink']
color_ar.extend(color_ar)
color_ar.extend(color_ar)

marker_ar = ['^', 'v', 'o', 's', '>', '<', '*', 'd', '8', 'p',  ]
marker_ar.extend(marker_ar)
marker_ar.extend(marker_ar)


figure0 = plt.figure(figsize=(10, 8))
plt.xlabel('Mstar (Msun/h)', fontsize=20)
plt.ylabel(r'$\eta$', fontsize=20)
plt.xscale('log')
plt.yscale('log')
#plt.xlim([7e9, 2.5e10])
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

count = 0
fart_ar = []
while (count < len(Mstar_ar)):
	fart, = plt.plot(Mstar_ar[count], eta_ar[count], color=color_ar[count], marker=marker_ar[count], ms=themarker,  ls='none')
	fart_ar.append(fart)
	if (dobadplots):
		plt.plot(Mstar_ar_bad[count], eta_ar_bad[count], color=color_ar[count], marker='_', ms=themarker,  ls='none')
	count += 1
foutname=today+'Mstar-cor.pdf'
plt.legend(fart_ar,thar, loc = 'best', prop={'size':10})
plt.savefig(foutname)
plt.clf()


figure1 = plt.figure(figsize=(10, 8))
plt.xlabel('escape velocity (km/s)', fontsize=20)
plt.ylabel(r'$\eta$', fontsize=20)
plt.title(r'$\alpha=$'+str(vescfit[0][1]))
plt.xscale('log')
plt.yscale('log')
#plt.xlim([250, 1000])
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

count = 0
fart_ar = []
while (count < len(Mstar_ar)):
	fart, =plt.plot(escape_v_ar[count], eta_ar[count], color=color_ar[count], marker=marker_ar[count], ms=themarker,  ls='none')
	fart_ar.append(fart)
	if (dobadplots):
		plt.plot(escape_v_ar_bad[count], eta_ar_bad[count], color=color_ar[count], marker='_', ms=themarker,  ls='none')
	count += 1
foutname=today+'vesc-cor.pdf'
plt.legend(fart_ar,thar, loc = 'best', prop={'size':10})
plt.savefig(foutname)
plt.clf()

weirdquant = 10**vescfit[0][0]*vesc_ext_ar**vescfit[0][1]
weirdquant_ext_log = np.log10(weirdquant)
weirdquant_ext_log 
[chisquared_fit1, RMSE_fit] = SaSta.sasha_chisquared(weirdquant_ext_log, eta_ext_ar_log, degree_of_freedom=1)
print 'chi squared, RMSE of vesc fit:',chisquared_fit1, RMSE_fit



figure2 = plt.figure(figsize=(10, 8))
plt.xlabel(r'$\Sigma_{SFR}$ (Msun/kpc^2/yr)', fontsize=20)
plt.ylabel(r'$\eta$', fontsize=20)
plt.title(r'$\alpha=$'+str(surffit[0][1]))
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

count = 0
fart_ar = []
while (count < len(Mstar_ar)):
	fart, =plt.plot(surf_ar[count], eta_ar[count], color=color_ar[count], marker=marker_ar[count], ms=themarker, ls='none')
	fart_ar.append(fart)
	if (dobadplots):
		plt.plot(surf_ar_bad[count], eta_ar_bad[count], color=color_ar[count], marker='_', ms=themarker,  ls='none')
	count += 1
foutname=today+'surf.pdf'
plt.legend(fart_ar,thar, loc = 'best', prop={'size':10})
plt.savefig(foutname)
weirdquant = 10**surffit[0][0]*surf_ext_ar**surffit[0][1]
weirdquant_ext_log = np.log10(weirdquant)
weirdquant_ext_log 
[chisquared_fit1, RMSE_fit] = SaSta.sasha_chisquared(weirdquant_ext_log, eta_ext_ar_log, degree_of_freedom=1)
print 'chi squared, RMSE of surf fit:',chisquared_fit1, RMSE_fit




figure4 = plt.figure(figsize=(10, 8))
plt.xlabel(r'$\dot{M}_{in}$ (Msun/yr)', fontsize=20)
plt.ylabel(r'$\eta$', fontsize=20)
plt.title(r'$\alpha=$'+str(Inflowfit[0][1]))
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

count = 0
fart_ar = []
while (count < len(Mstar_ar)):
	fart, =plt.plot(Inflow_rate_ar[count], eta_ar[count], color=color_ar[count], marker=marker_ar[count], ms=themarker,  ls='none')
	fart_ar.append(fart)
	if (dobadplots):
		plt.plot(Inflow_rate_ar_bad[count], eta_ar_bad[count], color=color_ar[count], marker='_', ms=themarker,  ls='none')
	count += 1
foutname=today+'inflow.pdf'
plt.legend(fart_ar,thar, loc = 'best', prop={'size':10})
plt.savefig(foutname)

weirdquant = 10**Inflowfit[0][0]*Inflow_ext_ar**Inflowfit[0][1]
weirdquant_ext_log = np.log10(weirdquant)
weirdquant_ext_log 
[chisquared_fit1, RMSE_fit] = SaSta.sasha_chisquared(weirdquant_ext_log, eta_ext_ar_log, degree_of_freedom=1)
print 'chi squared, RMSE of Inflow fit:',chisquared_fit1, RMSE_fit



figure5 = plt.figure(figsize=(10, 8))
plt.xlabel(r'$\dot{M}_{*}$ (Msun/yr)', fontsize=20)
plt.ylabel(r'$\eta$', fontsize=20)
plt.title(r'$\alpha=$'+str(SFRfit[0][1]))

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

count = 0
fart_ar = []
while (count < len(Mstar_ar)):
	fart, =plt.plot(SFR_ar[count], eta_ar[count], color=color_ar[count], marker=marker_ar[count], ms=themarker,  ls='none')
	fart_ar.append(fart)
	if (dobadplots):
		plt.plot(SFR_ar_bad[count], eta_ar_bad[count], color=color_ar[count], marker='_', ms=themarker,  ls='none')
	count += 1
foutname=today+'SFR.pdf'
plt.legend(fart_ar,thar, loc = 'best', prop={'size':10})
plt.savefig(foutname)
weirdquant = 10**SFRfit[0][0]*SFR_ext_ar**SFRfit[0][1]
weirdquant_ext_log = np.log10(weirdquant)
weirdquant_ext_log 
[chisquared_fit1, RMSE_fit] = SaSta.sasha_chisquared(weirdquant_ext_log, eta_ext_ar_log, degree_of_freedom=1)
print 'chi squared, RMSE of SFR fit:',chisquared_fit1, RMSE_fit



xspace = np.linspace(1,1000,50)
yspace = xspace
figure6 = plt.figure(figsize=(10, 8))
plt.ylabel(r'$\dot{M}_{*}$ (Msun/yr)', fontsize=20)
plt.xlabel(r'$\dot{M}_{in}$ (Msun/yr)', fontsize=20)
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

count = 0
fart_ar = []
while (count < len(Mstar_ar)):
	fart, = plt.plot(Inflow_rate_ar[count], SFR_ar[count], color=color_ar[count], marker=marker_ar[count], ms=themarker,  ls='none')
	fart_ar.append(fart)
	count += 1
plt.plot(xspace, yspace, ':r')
foutname=today+'SFR_vs_Inflow.pdf'
plt.legend(fart_ar,thar, loc = 'best', prop={'size':10})
plt.savefig(foutname)



###start of weird quants

weirdquant = 10**hugefit[0][0] * escape_v_ar**hugefit[0][1] * surf_ar**(hugefit[0][2]) 
weirdquant_bad = 10**hugefit[0][0] * escape_v_ar_bad**hugefit[0][1] * surf_ar_bad**(hugefit[0][2]) 

weirdquant_ext = 10**hugefit[0][0] * vesc_ext_ar**hugefit[0][1] * surf_ext_ar**(hugefit[0][2]) 
weirdquant_ext_log = np.log10(weirdquant_ext)

figure8 = plt.figure(figsize=(10, 8))
plt.ylabel(r'$\eta$', fontsize=20)
plt.xlabel(r'C$v_{esc}^\alpha \Sigma_{SFR}^\beta$',fontsize=20)
plt.title(r'$\alpha=$'+str(hugefit[0][1])+r' $\beta=$'+str(hugefit[0][2]))

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

count = 0
fart_ar = []
while (count < len(Mstar_ar)):
	fart, = plt.plot(weirdquant[count], eta_ar[count], color=color_ar[count], marker=marker_ar[count], ms=themarker,  ls='none')
	fart_ar.append(fart)
	if (dobadplots):
		plt.plot(weirdquant_bad[count], eta_ar_bad[count], color=color_ar[count], marker='_', ms=themarker,  ls='none')
	count += 1
foutname=today+'fit2_vesc-surf.pdf'
plt.legend(fart_ar,thar, loc = 'best', prop={'size':10})
plt.savefig(foutname)

[chisquared_fit1, RMSE_fit] = SaSta.sasha_chisquared(weirdquant_ext_log, eta_ext_ar_log)
print 'chi squared, RMSE of vesc surf fit:',chisquared_fit1, RMSE_fit



weirdquant = 10**hugefit2[0][0] * escape_v_ar**(hugefit2[0][1]) * Inflow_rate_ar**(hugefit2[0][2]) 
weirdquant_bad = 10**hugefit2[0][0] * escape_v_ar_bad**(hugefit2[0][1]) * Inflow_rate_ar_bad**(hugefit2[0][2]) 

weirdquant_ext = 10**hugefit2[0][0] * vesc_ext_ar**(hugefit2[0][1]) * Inflow_ext_ar**(hugefit2[0][2]) 
weirdquant_ext_log = np.log10(weirdquant_ext)

figure9 = plt.figure(figsize=(10, 8))
plt.ylabel(r'$\eta$', fontsize=20)
plt.xlabel(r'C$v_{esc}^\alpha \dot{M}_{in}^\beta$',fontsize=20)
plt.title(r'$\alpha=$'+str(hugefit2[0][1])+r' $\beta=$'+str(hugefit2[0][2]))

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

count = 0
fart_ar = []
while (count < len(Mstar_ar)):
	fart, = plt.plot(weirdquant[count], eta_ar[count], color=color_ar[count], marker=marker_ar[count], ms=themarker,  ls='none')
	fart_ar.append(fart)
	if (dobadplots):
		plt.plot(weirdquant_bad[count], eta_ar_bad[count], color=color_ar[count], marker='_', ms=themarker,  ls='none')
	count += 1
foutname=today+'fit3_vesc-flow.pdf'
plt.legend(fart_ar,thar, loc = 'best', prop={'size':10})
plt.savefig(foutname)

[chisquared_fit1, RMSE_fit] = SaSta.sasha_chisquared(weirdquant_ext_log, eta_ext_ar_log)
print 'chi squared, RMSE of vesc inflow fit:',chisquared_fit1, RMSE_fit



weirdquant = 10**hugefit3[0][0]* surf_ar**(hugefit3[0][1])  * Inflow_rate_ar**(hugefit3[0][2]) 
weirdquant_bad = 10**hugefit3[0][0]* surf_ar_bad**(hugefit3[0][1])  * Inflow_rate_ar_bad**(hugefit3[0][2]) 

weirdquant_ext = 10**hugefit3[0][0]* surf_ext_ar**(hugefit3[0][1])  * Inflow_ext_ar**(hugefit3[0][2]) 
weirdquant_ext_log = np.log10(weirdquant_ext)

figure10 = plt.figure(figsize=(10, 8))
plt.ylabel(r'$\eta$', fontsize=20)
plt.xlabel(r'C$\Sigma_{SFR}^\alpha \dot{M}_{in}^\beta$',fontsize=20)
plt.title(r'$\alpha=$'+str(hugefit3[0][1])+r' $\beta=$'+str(hugefit3[0][2]))

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

count = 0
fart_ar = []
while (count < len(Mstar_ar)):
	fart, = plt.plot(weirdquant[count], eta_ar[count], color=color_ar[count], marker=marker_ar[count], ms=themarker,  ls='none')
	fart_ar.append(fart)
	if (dobadplots):
		plt.plot(weirdquant_bad[count], eta_ar_bad[count], color=color_ar[count], marker='_', ms=themarker,  ls='none')
	count += 1
foutname=today+'fit4_surf-flow.pdf'
plt.legend(fart_ar,thar, loc = 'best', prop={'size':10})
plt.savefig(foutname)
[chisquared_fit1, RMSE_fit] = SaSta.sasha_chisquared(weirdquant_ext_log, eta_ext_ar_log)
print 'chi squared, RMSE of surf inflow fit:',chisquared_fit1, RMSE_fit

weirdquant = 10**megafit[0][0]*escape_v_ar**(megafit[0][1]) * surf_ar**(megafit[0][2]) * Inflow_rate_ar**(megafit[0][3] )
weirdquant_bad = 10**megafit[0][0]*escape_v_ar_bad**(megafit[0][1]) * surf_ar_bad**(megafit[0][2]) * Inflow_rate_ar_bad**(megafit[0][3] )

weirdquant_ext = 10**megafit[0][0]*vesc_ext_ar**(megafit[0][1]) * surf_ext_ar**(megafit[0][2]) * Inflow_ext_ar**(megafit[0][3] )
weirdquant_ext_log = np.log10(weirdquant_ext)

figure11 = plt.figure(figsize=(10, 8))
plt.ylabel(r'$\eta$', fontsize=20)
plt.xlabel(r'C$v_{esc}^\alpha \Sigma_{SFR}^\beta \dot{M}_{in}^\gamma$',fontsize=20)
plt.title(r'$\alpha=$'+str(megafit[0][1])+r' $\beta=$'+str(megafit[0][2])+r' $\gamma=$'+str(megafit[0][3]))
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

count = 0
fart_ar = []
while (count < len(Mstar_ar)):
	fart, = plt.plot(weirdquant[count], eta_ar[count], color=color_ar[count], marker=marker_ar[count], ms=themarker,  ls='none')
	fart_ar.append(fart)
	if (dobadplots):
		plt.plot(weirdquant_bad[count], eta_ar_bad[count], color=color_ar[count], marker='_', ms=themarker,  ls='none')
	count += 1
foutname=today+'fit5_allthree.pdf'
plt.legend(fart_ar,thar, loc = 'best', prop={'size':10})
plt.savefig(foutname)

[chisquared_fit1, RMSE_fit] = SaSta.sasha_chisquared(weirdquant_ext_log, eta_ext_ar_log,degree_of_freedom=3)

print 'chi squared, RMSE of mega fit :',chisquared_fit1, RMSE_fit



weirdquant = 10**megafit2[0][0]*escape_v_ar**(megafit2[0][1]) * surf_ar**(megafit2[0][2]) * SFR_ar**(megafit2[0][3] )
weirdquant_bad = 10**megafit2[0][0]*escape_v_ar_bad**(megafit2[0][1]) * surf_ar_bad**(megafit2[0][2]) * SFR_ar_bad**(megafit2[0][3] )

weirdquant_ext = 10**megafit2[0][0]*vesc_ext_ar**(megafit2[0][1]) * surf_ext_ar**(megafit2[0][2]) * SFR_ext_ar**(megafit2[0][3] )
weirdquant_ext_log = np.log10(weirdquant_ext)

figure12 = plt.figure(figsize=(10, 8))
plt.ylabel(r'$\eta$', fontsize=20)
plt.xlabel(r'C$v_{esc}^\alpha \Sigma_{SFR}^\beta \dot{M}_{*}^\gamma$',fontsize=20)
plt.title(r'$\alpha=$'+str(megafit2[0][1])+r' $\beta=$'+str(megafit2[0][2])+r' $\gamma=$'+str(megafit2[0][3]))

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

count = 0
fart_ar = []
while (count < len(Mstar_ar)):
	fart, = plt.plot(weirdquant[count], eta_ar[count], color=color_ar[count], marker=marker_ar[count], ms=themarker,  ls='none')
	fart_ar.append(fart)
	if (dobadplots):
		plt.plot(weirdquant_bad[count], eta_ar_bad[count], color=color_ar[count], marker='_', ms=themarker,  ls='none')
	count += 1
foutname=today+'fit6_withSFR.pdf'
plt.legend(fart_ar,thar, loc = 'best', prop={'size':10})
plt.savefig(foutname)

[chisquared_fit1, RMSE_fit] = SaSta.sasha_chisquared(weirdquant_ext_log, eta_ext_ar_log,degree_of_freedom=3)

print 'chi squared, RMSE of mega fit with SFR :',chisquared_fit1, RMSE_fit


figure12b = plt.figure(figsize=(10, 8))
outflow_ar = eta_ar * SFR_ar
outflow_ar_bad = eta_ar_bad *SFR_ar_bad
plt.ylabel('Outflow rate (Msun/yr)', fontsize=20)
plt.xlabel(r'C$v_{esc}^\alpha \Sigma_{SFR}^\beta \dot{M}_{*}^\gamma$',fontsize=20)
plt.title(r'$\alpha=$'+str(megafit2[0][1])+r' $\beta=$'+str(megafit2[0][2])+r' $\gamma=$'+str(megafit2[0][3]))

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

count = 0
fart_ar = []
while (count < len(Mstar_ar)):
	fart, = plt.plot(weirdquant[count], outflow_ar[count], color=color_ar[count], marker=marker_ar[count], ms=themarker,  ls='none')
	fart_ar.append(fart)
	if (dobadplots):
		plt.plot(weirdquant_bad[count], outflow_ar_bad[count], color=color_ar[count], marker='_', ms=themarker,  ls='none')
	count += 1
foutname=today+'fit6_withSFRoutflow.pdf'
plt.legend(fart_ar,thar, loc = 'best', prop={'size':10})
plt.savefig(foutname)

[chisquared_fit1, RMSE_fit] = SaSta.sasha_chisquared(weirdquant_ext_log, eta_ext_ar_log,degree_of_freedom=3)

print 'chi squared, RMSE of mega fit with SFR :',chisquared_fit1, RMSE_fit