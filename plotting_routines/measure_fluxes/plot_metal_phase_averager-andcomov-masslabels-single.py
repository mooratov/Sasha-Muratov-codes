import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.interpolate import interp1d
from datetime import date

today = date.today()
print today
today = str(today)

plt.rcParams['ps.fonttype'] = 42 
plt.rcParams['axes.linewidth'] = 2.0 #set the value globally
plt.rcParams['legend.numpoints'] = 1
Mstar_wanted = [1e7, 1e8, 1e9, 1e10, 3e10]
Mstar_wanted_str =[r'$10^7$', r'$10^8$', r'$10^9$', r'$10^{10}$', r'$3\times10^{10}$']


legsize=17
themarkersize = 20
empirical_yield = True
no_zach = False
z2mode = True
if (z2mode):
	today+='z2'
	Mstar_wanted = [1e6, 1e7, 1e8, 1e9, 1e10]
	Mstar_wanted_str =[r'$10^6$', r'$10^7$', r'$10^8$', r'$10^9$', r'$10^{10}$']

#finname = 'MetalTableAlpha.txt'

zstart = 2.0
zend = 0.0

if (len(sys.argv) < 3):
	print 'syntax  blah.py filename filename2 [zstart zend]'
	print 'ideally filename is ', 'guide_of_files.txt'
	print 'filneame2 is guide_of_files-comov.txt'
	sys.exit()

finname = str(sys.argv[1])
finname2 = str(sys.argv[2])

if  (len(sys.argv) > 3):
	zstart = float(sys.argv[3])
	zend = float(sys.argv[4])

astart = 1.0 / (1.0 + zstart)
aend = 1.0 / (1.0 + zend)

#finname  = 'guide_of_files.txt'


f = open(finname)
shorp = f.readlines()

Mstar = []
Mexpected = []
Minstars = []
MinISM = []
MinHot =[]
MinO6 = []
MinLowIon = []
MinCold = []
thea = []
theyield = []


for name in shorp:
	dars = np.loadtxt(name.rstrip(), ndmin=2)
	saag = dars.shape[1]
	athea = dars[:,-2]
	cuthigh = (athea <= aend)
	cutlow = (athea >= astart)
	cuts  = cuthigh*cutlow
	colcount = 0
	sline = []
	while (colcount < saag):
		fordocut = (dars[:,colcount][cuts] > 0 ) * np.isreal(dars[:,colcount][cuts])
		if (len(dars[:,colcount][cuts][fordocut] > 0)):
			newavgval = np.mean(dars[:,colcount][cuts][fordocut])
		else:
			newavgval = 0
		#newavgval = 10**np.mean(np.log10(dars[:,colcount][cuts][fordocut]))
		sline.append(newavgval)
		colcount+=1
	Mstar.append(sline[0])
	Mexpected.append(sline[1])
	Minstars.append(sline[2])
	MinISM.append(sline[3])
	MinHot.append(sline[4])
	MinO6.append(sline[5])
	MinLowIon.append(sline[6])
	MinCold.append(sline[7])
	thea.append(sline[-2])	
	theyield.append(sline[-3])


Mstar = np.array(Mstar)
Mexpected = np.array(Mexpected)
Minstars = np.array(Minstars)
MinISM = np.array(MinISM)
MinHot = np.array(MinHot)
MinO6 = np.array(MinO6)
MinLowIon = np.array(MinLowIon)
MinCold = np.array(MinCold)
thea = np.array(thea)
theyield  = np.array(theyield)

Mstar_label_ar = []
for ara in Mstar:
	Mstar_label_ar.append(str(float('{:0.3g}'.format(np.log10(ara)))))
Mstar_label_ar = np.array(Mstar_label_ar)

Mexp_label_ar = []
for ara in Mexpected:
	Mexp_label_ar.append(str(float('{:0.3g}'.format(np.log10(ara)))))
Mexp_label_ar = np.array(Mexp_label_ar)


MinHalo = MinHot + MinO6 + MinLowIon + MinCold
TotalM = MinHalo + MinISM + Minstars
f.close()

#now get comoving stuff
f = open(finname2)
shorp = f.readlines()

CMstar = []
CMexpected = []
CMinstars = []
CMinISM = []
CMinHot =[]
CMinO6 = []
CMinLowIon = []
CMinCold = []
Cthea = []
Ctheyield = []


for name in shorp:
	dars = np.loadtxt(name.rstrip(), ndmin=2)
	saag = dars.shape[1]
	athea = dars[:,-2]
	cuthigh = (athea <= aend)
	cutlow = (athea >= astart)
	cuts  = cuthigh*cutlow
	colcount = 0
	sline = []
	while (colcount < saag):
		fordocut = (dars[:,colcount][cuts] > 0 ) * np.isreal(dars[:,colcount][cuts])
		if (len(dars[:,colcount][cuts][fordocut] > 0)): 
			newavgval = np.mean(dars[:,colcount][cuts][fordocut])
		else:
			newavgval = 0
		sline.append(newavgval)
		colcount+=1
	CMstar.append(sline[0])
	CMexpected.append(sline[1])
	CMinstars.append(sline[2])
	CMinISM.append(sline[3])
	CMinHot.append(sline[4])
	CMinO6.append(sline[5])
	CMinLowIon.append(sline[6])
	CMinCold.append(sline[7])
	Cthea.append(sline[-2])	
	Ctheyield.append(sline[-3])


CMstar = np.array(CMstar)
CMexpected = np.array(CMexpected)
CMinstars = np.array(CMinstars)
CMinISM = np.array(CMinISM)
CMinHot = np.array(CMinHot)
CMinO6 = np.array(CMinO6)
CMinLowIon = np.array(CMinLowIon)
CMinCold = np.array(CMinCold)
Cthea = np.array(Cthea)
Ctheyield  = np.array(Ctheyield)
	
	
CMinHalo = CMinHot + CMinO6 + CMinLowIon + CMinCold
CTotalM = CMinHalo + CMinISM + CMinstars

labelfont = 25
ticklength = 12
tickwidth = 3.2
titlefont = 26



nums = np.arange(len(Mstar))


if (no_zach):
	labels = [ r'm10', r'm11', r'm12v', r'm12i', r'm12q', r'${\rm m11_{h383}}$', r'${\rm m10_{h573}}$', r'${\rm m10_{h1146}}$', r'${\rm m10_{h1297}}$']
else:
	labels = [ r'${\rm m10}$', r'${\rm m11}$', r'${\rm m12v}$', r'${\rm m12i}$', r'${\rm m12q}$', r'${\rm m11_{h383}}$', r'${\rm m10_{h573}}$', r'${\rm m10_{h1146}}$', r'${\rm m10_{h1297}}$', r'${\rm m11.4a}$', r'${\rm m11.9a}$']

if (z2mode):
	labels = labels + [r'$z2h_{506}$', r'$z2h_{830}$']


labels = np.array(labels)
Mstar/=0.7
CMstar/=0.7
if (empirical_yield):
	Mexpected = theyield*Mstar
	CMexpected = Ctheyield*CMstar

fartoon = np.argsort(Mstar)

totalMeta = Minstars + MinHalo + MinISM

fexpected = totalMeta/Mexpected

fartoon2 = np.argsort(fexpected)
fexp_wanted = [0.4, 0.6, 0.8, 1.0]
Mstarexp_wanted_str =['40%', '60%', '80%', '100%']


fig1 = plt.figure(figsize=(10,8))
plt.plot(nums, Mexpected, 'sk', fillstyle='none', ms=themarkersize)
plt.plot(nums, Minstars, 'og', ms=themarkersize)
plt.plot(nums, MinISM, '*k', ms =themarkersize)
plt.plot(nums, MinHot, '^r', ms=themarkersize)
plt.plot(nums, MinO6, '^k', ms=themarkersize)
plt.plot(nums, MinLowIon, '^b', ms=themarkersize)
plt.plot(nums, MinCold, '^m', ms=themarkersize)
plt.xticks(nums, labels, rotation='vertical')
plt.margins(0.2)
plt.yscale('log')
plt.ylim([10**1, 10**9.5])
plt.ylabel(r'Mass of Metals ($M_\odot$)', fontsize=titlefont)
plt.xlabel(r'${\rm log_10} (M_* /\odot)')
plt.yscale('log')

ax = plt.gca()

ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax.yaxis.set_tick_params(width=tickwidth, length=ticklength)


thar22b =['Expected','Locked in Stars','ISM (r<0.1 Rvir)','Hot Phase(logT>5.3)','OVI-like ions (4.7<logT<5.3)', 'Low Ions (4.0<logT<4.7)', 'Cold (logT<4.0)']

#plt.legend(thar22b, loc = 'best',  prop={'size':12})
plt.legend(thar22b, loc = 'best',  prop={'size':12})
ax = plt.gca()

ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax.yaxis.set_tick_params(width=tickwidth, length=ticklength)
plt.savefig(today + 'metaldiagram'+'.pdf')
plt.clf()


fig1 = plt.figure(figsize=(10,8))
plt.plot(Mstar, Mexpected, 'sk', fillstyle='none', ms=themarkersize)
plt.plot(Mstar, Minstars, 'og', ms=themarkersize)
plt.plot(Mstar, MinISM, '*k', ms =themarkersize)
plt.plot(Mstar, MinHot, '^r', ms=themarkersize)
plt.plot(Mstar, MinO6, '^k', ms=themarkersize)
plt.plot(Mstar, MinLowIon, '^b', ms=themarkersize)
plt.plot(Mstar, MinCold, '^m', ms=themarkersize)
#plt.xticks(nums, labels, rotation='vertical')
plt.margins(0.2)
plt.yscale('log')
plt.xscale('log')
plt.ylim([10**1, 10**9.5])
plt.xlim([10**6, 10**11])
plt.ylabel(r'Mass of Metals ($M_\odot$)', fontsize=titlefont)
plt.xlabel(r'Stellar Mass at z=0 ($M_\odot$)', fontsize=titlefont)
plt.yscale('log')

ax = plt.gca()

ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax.yaxis.set_tick_params(width=tickwidth, length=ticklength)


thar22b =['Expected','Locked in Stars','ISM (r<0.1 Rvir)','Hot Phase(logT>5.3)','OVI-like Ions (4.7<logT<5.3)', 'Low Ions (4.0<logT<4.7)', 'Cold (logT<4.0)']

#plt.legend(thar22b, loc = 'best',  prop={'size':12})
plt.legend(thar22b, loc = 'best',  prop={'size':12})
ax = plt.gca()

ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax.yaxis.set_tick_params(width=tickwidth, length=ticklength)

plt.savefig(today + 'Mstar_metaldiagram'+'.pdf')

fig1 = plt.figure(figsize=(10,9))

totalMeta = Minstars + MinHalo + MinISM

count = 0
print 'name Mstar Mexpected, Minstars, Minhalo, MinISM, totalMeta'
while (count < len(totalMeta)):
	print labels[count], "{0:.6g}".format(Mstar[count]), "{0:.6g}".format(Mexpected[count]), "{0:.6g}".format(Minstars[count]), "{0:.6g}".format(MinHalo[count]), "{0:.6g}".format(MinISM[count]), "{0:.6g}".format(totalMeta[count])
	count += 1 

plt.plot(Mstar, Mexpected, 'sk', fillstyle='none', ms=themarkersize,mew=2)
plt.plot(Mstar, Minstars, 'og', ms=themarkersize, alpha=0.75)
plt.plot(Mstar, MinHalo, 'pm', ms=themarkersize, alpha=0.75)
plt.plot(Mstar, MinISM, '*b', ms =themarkersize, alpha=0.75)
plt.plot(Mstar, totalMeta, '+r', ms=themarkersize, mew=2)


#plt.xticks(nums, labels, rotation='vertical')
plt.margins(0.2)
plt.yscale('log')
plt.xscale('log')
plt.ylim([10**2.5, 10**9.5])
plt.xlim([10**6, 10**11])
plt.ylabel(r'Mass of Metals ($M_\odot$)', fontsize=titlefont)
plt.xlabel(r'Stellar Mass at z=0 ($M_\odot$)', fontsize=titlefont)
plt.yscale('log')

ax = plt.gca()

ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax.yaxis.set_tick_params(width=tickwidth, length=ticklength)


thar22b =['Expected','Locked in Stars','Halo CGM','ISM (r<0.1 Rvir)', 'total in halo']

#plt.legend(thar22b, loc = 'best',  prop={'size':12})
plt.legend(thar22b, loc = 'best',  prop={'size':legsize})
ax = plt.gca()

ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax.yaxis.set_tick_params(width=tickwidth, length=ticklength)

plt.savefig(today + 'Mstar_metalbasics'+'.pdf')
plt.clf()



fig1 = plt.figure(figsize=(10,9))

totalMeta = Minstars + MinHalo + MinISM

count = 0
print 'name Mstar Mexpected, Minstars, Minhalo, MinISM, totalMeta'
while (count < len(totalMeta)):
	print labels[count], "{0:.6g}".format(Mstar[count]), "{0:.6g}".format(Mexpected[count]), "{0:.6g}".format(Minstars[count]), "{0:.6g}".format(MinHalo[count]), "{0:.6g}".format(MinISM[count]), "{0:.6g}".format(totalMeta[count])
	count += 1 

plt.plot(Mexpected, Minstars, 'og', ms=themarkersize, alpha=0.75)
plt.plot(Mexpected, MinHalo, 'pm', ms=themarkersize, alpha=0.75)
plt.plot(Mexpected, MinISM, '*b', ms =themarkersize, alpha=0.75)
plt.plot(Mexpected, totalMeta, '+r', ms=themarkersize, mew=2)

xspace = np.arange(3, 12, 0.1)
xspace = 10**xspace
yspace = xspace
plt.plot(xspace, yspace, '--r')
yspace = yspace/2.0
plt.plot(xspace, yspace, ':b')

#plt.xticks(nums, labels, rotation='vertical')
plt.margins(0.2)
plt.yscale('log')
plt.xscale('log')
plt.ylim([10**2.5, 10**9.5])
plt.xlim([10**4, 10**11])
plt.ylabel(r'Mass of Metals ($M_\odot$)', fontsize=titlefont)
plt.xlabel(r'Expected Metals ($M_\odot$)', fontsize=titlefont)
plt.yscale('log')

ax = plt.gca()

ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax.yaxis.set_tick_params(width=tickwidth, length=ticklength)


thar22b =['Locked in Stars','Halo CGM','ISM (r<0.1 Rvir)', 'total in halo']

#plt.legend(thar22b, loc = 'best',  prop={'size':12})
plt.legend(thar22b, loc = 'best',  prop={'size':legsize})
ax = plt.gca()

ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax.yaxis.set_tick_params(width=tickwidth, length=ticklength)

plt.savefig(today + 'Mexp_metalbasics'+'.pdf')
plt.clf()





fig1 = plt.figure(figsize=(10,8))
plt.plot(Mstar, MinHot/MinHalo, '^r', ms=themarkersize)
plt.plot(Mstar, MinO6/MinHalo, '^k', ms=themarkersize)
plt.plot(Mstar, MinLowIon/MinHalo, '^b', ms=themarkersize)
plt.plot(Mstar, MinCold/MinHalo, '^m', ms=themarkersize)
#plt.xticks(nums, labels, rotation='vertical')
plt.margins(0.2)
plt.yscale('log')
plt.xscale('log')
#plt.ylim([10**1, 10**9.5])
plt.ylim([1e-3,1.0])
plt.xlim([10**6, 10**11])
plt.xlabel(r'Stellar Mass at z=0 ($M_\odot$)', fontsize=titlefont)
plt.ylabel('fraction of CGM metals per phase', fontsize=titlefont)
plt.yscale('log')

ax = plt.gca()

ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax.yaxis.set_tick_params(width=tickwidth, length=ticklength)


thar22b =['Hot Phase(logT>5.3)','OVI-like Ions (4.7<logT<5.3)', 'Low Ions (4.0<logT<4.7)', 'Cold (logT<4.0)']

#plt.legend(thar22b, loc = 'best',  prop={'size':12})
plt.legend(thar22b, loc = 'best',  prop={'size':legsize})
ax = plt.gca()

ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax.yaxis.set_tick_params(width=tickwidth, length=ticklength)

plt.savefig(today + 'Mstar_halodiagram'+'.pdf')


plt.clf()

CtotalMeta = CMinstars + CMinHalo + CMinISM


fig1b = plt.figure(figsize=(11.5,12.5))
mex1, = plt.plot(nums, Mexpected[fartoon], 'sk', fillstyle='none',  ms=themarkersize, mew=2)
mst1, = plt.plot(nums, Minstars[fartoon], 'ob', ms=themarkersize)
mcgm1, = plt.plot(nums, MinHalo[fartoon], '^b', ms=themarkersize)
mism1, = plt.plot(nums, MinISM[fartoon], 'vb', ms =themarkersize)
mtot1, = plt.plot(nums, totalMeta[fartoon], '+r', ms=themarkersize, mew=2)

#mex2, = plt.plot(nums+0.5, CMexpected[fartoon], 'sk', fillstyle='none', ms=themarkersize)
#mst2, = plt.plot(nums+0.5, CMinstars[fartoon], 'ob',fillstyle='none', ms=themarkersize)
#mcgm2, = plt.plot(nums+0.5, CMinHalo[fartoon], '^b',fillstyle='none', ms=themarkersize)
#mism2, = plt.plot(nums+0.5, CMinISM[fartoon], 'vb',fillstyle='none', ms =14)
#mtot2, = plt.plot(nums+0.5, CtotalMeta[fartoon], '+r',fillstyle='none', ms=themarkersize)

#plt.xticks(nums, labels, rotation='vertical')
plt.yscale('log')
plt.ylim([10**2.5, 10**9.5])
plt.ylabel(r'Mass of Metals $(M_\odot)$', fontsize=titlefont)
#plt.xlabel('Stellar Mass at z=0 (Msun)', fontsize=20)
plt.yscale('log')
plt.xticks(nums, labels[fartoon], rotation='vertical', fontsize=titlefont)
#plt.xlabel(r'${\rm log_{10}} (M_* /M_{\odot})$', fontsize=titlefont)

plt.margins(0.05)

ax = plt.gca()

ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax.yaxis.set_tick_params(width=tickwidth, length=ticklength)


thar22b =['Expected','Locked in Stars','CGM','ISM', 'Total in Halo']

#plt.legend(thar22b, loc = 'best',  prop={'size':12})

plt.legend( [(mex1), (mst1), (mcgm1), (mism1), (mtot1)],thar22b, loc = 'best',  prop={'size':legsize})
#plt.legend(thar22b, loc = 'best',  prop={'size':legsize})

ax = plt.gca()
ax2 = ax.twiny()

num_for_Mstar = interp1d(Mstar[fartoon], nums, bounds_error=False)
num_wanted = num_for_Mstar(Mstar_wanted)

ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax.yaxis.set_tick_params(width=tickwidth, length=ticklength)

ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(num_wanted)
ax2.set_xticklabels(Mstar_wanted_str)
ax2.set_xlabel(r'$M_* (M_{\odot})$', fontsize=titlefont)
ticklabels = ax2.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax2.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ax2.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax2.yaxis.set_tick_params(width=tickwidth, length=ticklength)   

plt.savefig(today + 'label_metalbasics'+'.pdf')
plt.clf()



fig1b = plt.figure(figsize=(11.5,12))
#mex1, = plt.plot(nums, Mexpected[fartoon], 'sk', fillstyle='none',  ms=themarkersize, mew=2)
mst1, = plt.plot(nums, Minstars[fartoon]/Mexpected[fartoon], 'og', ms=themarkersize*1.2)
mcgm1, = plt.plot(nums, MinHalo[fartoon]/Mexpected[fartoon], 'pm', ms=themarkersize*1.2)
mism1, = plt.plot(nums, MinISM[fartoon]/Mexpected[fartoon], '*b', ms =themarkersize*1.2)
mtot1, = plt.plot(nums, totalMeta[fartoon]/Mexpected[fartoon], '+r', ms=themarkersize*1.2, mew=3)

#mex2, = plt.plot(nums+0.5, CMexpected[fartoon], 'sk', fillstyle='none', ms=themarkersize)
#mst2, = plt.plot(nums+0.5, CMinstars[fartoon], 'ob',fillstyle='none', ms=themarkersize)
#mcgm2, = plt.plot(nums+0.5, CMinHalo[fartoon], '^b',fillstyle='none', ms=themarkersize)
#mism2, = plt.plot(nums+0.5, CMinISM[fartoon], 'vb',fillstyle='none', ms =14)
#mtot2, = plt.plot(nums+0.5, CtotalMeta[fartoon], '+r',fillstyle='none', ms=themarkersize)

#plt.xticks(nums, labels, rotation='vertical')
plt.yscale('log')
plt.ylim([10**-2.0, 10**0.1])
#plt.ylabel(r'Mass of Metals $(M_\odot)$', fontsize=titlefont)
plt.ylabel(r'fraction of expected metals', fontsize=titlefont)

#plt.xlabel('Stellar Mass at z=0 (Msun)', fontsize=20)
plt.yscale('log')
plt.xticks(nums, Mstar_label_ar[fartoon], rotation='vertical', fontsize=titlefont)
plt.xlabel(r'${\rm log_{10}} (M_* /M_{\odot})$', fontsize=titlefont)

plt.margins(0.05)

ax = plt.gca()

ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax.yaxis.set_tick_params(width=tickwidth, length=ticklength)


thar22b =['Locked in Stars','CGM','ISM', 'Total in Halo']

#plt.legend(thar22b, loc = 'best',  prop={'size':12})

plt.legend( [(mst1), (mcgm1), (mism1), (mtot1)],thar22b, loc = 'best',  prop={'size':legsize})
#plt.legend(thar22b, loc = 'best',  prop={'size':legsize})

ax = plt.gca()
ax2 = ax.twiny()

num_for_Mstar = interp1d(Mstar[fartoon], nums, bounds_error=False)
num_wanted = num_for_Mstar(Mstar_wanted)

ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax.yaxis.set_tick_params(width=tickwidth, length=ticklength)
ax.yaxis.set_tick_params(which='minor', length=ticklength/2.0)

ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(nums)
ax2.set_xticklabels(Mexp_label_ar[fartoon], rotation='vertical')
ax2.set_xlabel(r'${\rm log_{10}} (M_{Z,exp} /M_{\odot})$', fontsize=titlefont)
ax2.xaxis.set_label_coords(0.5, 1.12)
ticklabels = ax2.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax2.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ax2.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax2.yaxis.set_tick_params(width=tickwidth*2.5, length=ticklength*2.5)   

plt.tight_layout()
#plt.subplots_adjust(top=0.2)
plt.savefig(today + 'double_label_metalfracbasics'+'.pdf')
plt.clf()





fig1b = plt.figure(figsize=(11.5,12))
#mex1, = plt.plot(nums, Mexpected[fartoon], 'sk', fillstyle='none',  ms=themarkersize, mew=2)
mst1, = plt.plot(nums, Minstars[fartoon]/Mexpected[fartoon], 'og', ms=themarkersize*1.2)
mcgm1, = plt.plot(nums, MinHalo[fartoon]/Mexpected[fartoon], 'pm', ms=themarkersize*1.2)
mism1, = plt.plot(nums, MinISM[fartoon]/Mexpected[fartoon], '*b', ms =themarkersize*1.2)
mtot1, = plt.plot(nums, totalMeta[fartoon]/Mexpected[fartoon], '+r', ms=themarkersize*1.2, mew=3)

#mex2, = plt.plot(nums+0.5, CMexpected[fartoon], 'sk', fillstyle='none', ms=themarkersize)
#mst2, = plt.plot(nums+0.5, CMinstars[fartoon], 'ob',fillstyle='none', ms=themarkersize)
#mcgm2, = plt.plot(nums+0.5, CMinHalo[fartoon], '^b',fillstyle='none', ms=themarkersize)
#mism2, = plt.plot(nums+0.5, CMinISM[fartoon], 'vb',fillstyle='none', ms =14)
#mtot2, = plt.plot(nums+0.5, CtotalMeta[fartoon], '+r',fillstyle='none', ms=themarkersize)

#plt.xticks(nums, labels, rotation='vertical')
plt.yscale('linear')
plt.ylim([0, 1.1])
#plt.ylabel(r'Mass of Metals $(M_\odot)$', fontsize=titlefont)
plt.ylabel(r'fraction of expected metals', fontsize=titlefont)

#plt.xlabel('Stellar Mass at z=0 (Msun)', fontsize=20)
plt.yscale('linear')
plt.xticks(nums, Mstar_label_ar[fartoon], rotation='vertical', fontsize=titlefont)
plt.xlabel(r'${\rm log_{10}} (M_* /M_{\odot})$', fontsize=titlefont)

plt.margins(0.05)

ax = plt.gca()

ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax.yaxis.set_tick_params(width=tickwidth, length=ticklength)


thar22b =['Locked in Stars','CGM','ISM', 'Total in Halo']

#plt.legend(thar22b, loc = 'best',  prop={'size':12})

plt.legend( [(mst1), (mcgm1), (mism1), (mtot1)],thar22b, loc = 'best',  prop={'size':legsize})
#plt.legend(thar22b, loc = 'best',  prop={'size':legsize})

ax = plt.gca()
ax2 = ax.twiny()

num_for_Mstar = interp1d(Mstar[fartoon], nums, bounds_error=False)
num_wanted = num_for_Mstar(Mstar_wanted)

ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax.yaxis.set_tick_params(width=tickwidth, length=ticklength)
ax.yaxis.set_tick_params(which='minor', length=ticklength/2.0)

ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(nums)
ax2.set_xticklabels(Mexp_label_ar[fartoon], rotation='vertical')
ax2.set_xlabel(r'${\rm log_{10}} (M_{Z,exp} /M_{\odot})$', fontsize=titlefont)
ax2.xaxis.set_label_coords(0.5, 1.12)
ticklabels = ax2.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax2.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ax2.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax2.yaxis.set_tick_params(width=tickwidth*2.5, length=ticklength*2.5)   

plt.tight_layout()
#plt.subplots_adjust(top=0.2)
plt.savefig(today + 'double_label_metalfracbasics_linear'+'.pdf')
plt.clf()



fig1b = plt.figure(figsize=(11.5,12))
#mex1, = plt.plot(nums, Mexpected[fartoon], 'sk', fillstyle='none',  ms=themarkersize, mew=2)
#mst1, = plt.plot(nums, Minstars[fartoon]/Mexpected[fartoon], 'og', ms=themarkersize*1.2)
#mcgm1, = plt.plot(nums, MinHalo[fartoon]/Mexpected[fartoon], 'pm', ms=themarkersize*1.2)
#mism1, = plt.plot(nums, MinISM[fartoon]/Mexpected[fartoon], '*b', ms =themarkersize*1.2)
#mtot1, = plt.plot(nums, totalMeta[fartoon]/Mexpected[fartoon], '+r', ms=themarkersize*1.2, mew=3)

width = 0.6



p1 = plt.bar(nums, Minstars[fartoon]/Mexpected[fartoon], width, color='g', alpha=1)
p2 = plt.bar(nums,  MinISM[fartoon]/Mexpected[fartoon], width, color='b',
             bottom= Minstars[fartoon]/Mexpected[fartoon], alpha=1)
p3 = plt.bar(nums, MinHalo[fartoon]/Mexpected[fartoon], width, color='m',
             bottom=( Minstars[fartoon]/Mexpected[fartoon] + MinISM[fartoon]/Mexpected[fartoon]), alpha=1)
             
#p4 = plt.bar(nums, MinCold[fartoon]/MinHalo[fartoon], width, color='m',
#             bottom=(MinHot[fartoon]/MinHalo[fartoon] + MinO6[fartoon]/MinHalo[fartoon] + MinLowIon[fartoon]/MinHalo[fartoon]), alpha=1)




#plt.xticks(nums, labels, rotation='vertical')
plt.yscale('linear')
plt.ylim([0, 1.1])
#plt.ylabel(r'Mass of Metals $(M_\odot)$', fontsize=titlefont)
plt.ylabel(r'fraction of expected metals', fontsize=titlefont)

#plt.xlabel('Stellar Mass at z=0 (Msun)', fontsize=20)
plt.yscale('linear')
plt.xticks(nums+width/2.0, Mstar_label_ar[fartoon], rotation='vertical', fontsize=titlefont)
plt.xlabel(r'${\rm log_{10}} (M_* /M_{\odot})$', fontsize=titlefont)

plt.margins(0.05)

ax = plt.gca()

ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax.yaxis.set_tick_params(width=tickwidth, length=ticklength)


thar22b =['Locked in Stars','ISM', 'CGM']

#plt.legend(thar22b, loc = 'best',  prop={'size':12})

plt.legend( [(p1), (p2), (p3)],thar22b, loc = 'best',  prop={'size':legsize})
#plt.legend(thar22b, loc = 'best',  prop={'size':legsize})

ax = plt.gca()
ax2 = ax.twiny()

num_for_Mstar = interp1d(Mstar[fartoon], nums, bounds_error=False)
num_wanted = num_for_Mstar(Mstar_wanted)

ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax.yaxis.set_tick_params(width=tickwidth, length=ticklength)
ax.yaxis.set_tick_params(which='minor', length=ticklength/2.0)

ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(nums+width/2.0)
ax2.set_xticklabels(Mexp_label_ar[fartoon], rotation='vertical')
ax2.set_xlabel(r'${\rm log_{10}} (M_{Z,exp} /M_{\odot})$', fontsize=titlefont)
ax2.xaxis.set_label_coords(0.5, 1.12)
ticklabels = ax2.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax2.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ax2.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax2.yaxis.set_tick_params(width=tickwidth*2.5, length=ticklength*2.5)   

plt.tight_layout()
#plt.subplots_adjust(top=0.2)
plt.savefig(today + 'double_label_metalfracbars_linear'+'.pdf')
plt.clf()





plt.clf()

CtotalMeta = CMinstars + CMinHalo + CMinISM


fig1b = plt.figure(figsize=(11.5,12.5))
mex1, = plt.plot(nums, Mexpected[fartoon2], 'sk', fillstyle='none',  ms=themarkersize, mew=2)
mst1, = plt.plot(nums, Minstars[fartoon2], 'og', ms=themarkersize)
mcgm1, = plt.plot(nums, MinHalo[fartoon2], 'pm', ms=themarkersize)
mism1, = plt.plot(nums, MinISM[fartoon2], '*b', ms =themarkersize)
mtot1, = plt.plot(nums, totalMeta[fartoon2], '+r', ms=themarkersize, mew=2)

#mex2, = plt.plot(nums+0.5, CMexpected[fartoon], 'sk', fillstyle='none', ms=themarkersize)
#mst2, = plt.plot(nums+0.5, CMinstars[fartoon], 'ob',fillstyle='none', ms=themarkersize)
#mcgm2, = plt.plot(nums+0.5, CMinHalo[fartoon], '^b',fillstyle='none', ms=themarkersize)
#mism2, = plt.plot(nums+0.5, CMinISM[fartoon], 'vb',fillstyle='none', ms =14)
#mtot2, = plt.plot(nums+0.5, CtotalMeta[fartoon], '+r',fillstyle='none', ms=themarkersize)

#plt.xticks(nums, labels, rotation='vertical')
plt.yscale('log')
plt.ylim([10**2.5, 10**9.5])
plt.ylabel(r'Mass of Metals $(M_\odot)$', fontsize=titlefont)
#plt.xlabel('Stellar Mass at z=0 (Msun)', fontsize=20)
plt.yscale('log')
plt.xticks(nums, Mstar_label_ar[fartoon2], rotation='vertical', fontsize=titlefont)
plt.xlabel(r'${\rm log_{10}} (M_* /M_{\odot})$', fontsize=titlefont)

plt.margins(0.05)

ax = plt.gca()

ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax.yaxis.set_tick_params(width=tickwidth, length=ticklength)


thar22b =['Expected','Locked in Stars','CGM','ISM', 'Total in Halo']

#plt.legend(thar22b, loc = 'best',  prop={'size':12})

plt.legend( [(mex1), (mst1), (mcgm1), (mism1), (mtot1)],thar22b, loc = 'best',  prop={'size':legsize})
#plt.legend(thar22b, loc = 'best',  prop={'size':legsize})

ax = plt.gca()
ax2 = ax.twiny()

num_for_fexp = interp1d(fexpected[fartoon2], nums, bounds_error=False)
num_wanted = num_for_fexp(fexp_wanted)

ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax.yaxis.set_tick_params(width=tickwidth, length=ticklength)

ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(num_wanted)
ax2.set_xticklabels(Mstarexp_wanted_str)
ax2.set_xlabel(r'fraction of expected', fontsize=titlefont)
ticklabels = ax2.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax2.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ax2.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax2.yaxis.set_tick_params(width=tickwidth, length=ticklength)   

plt.savefig(today + 'label_metalbasics_fexp'+'.pdf')
plt.clf()







width = Mstar*0.5


fig5 = plt.figure(figsize=(10,8))
p1 = plt.bar(Mstar, MinHot/MinHalo, width, color='r', alpha=1)
p2 = plt.bar(Mstar, MinO6/MinHalo, width, color='k',
             bottom=MinHot/MinHalo, alpha=1)
p3 = plt.bar(Mstar, MinLowIon/MinHalo, width, color='b',
             bottom=(MinHot/MinHalo + MinO6/MinHalo), alpha=1)
p4 = plt.bar(Mstar, MinCold/MinHalo, width, color='m',
             bottom=(MinHot/MinHalo + MinO6/MinHalo + MinLowIon/MinHalo), alpha=1)
#Mstar[-1]-=1e10
 
            
#plt.plot(Mstar, MinHot/MinHalo, '^r', ms=themarkersize)
#plt.plot(Mstar, MinO6/MinHalo, '^k', ms=themarkersize)
#plt.plot(Mstar, MinLowIon/MinHalo, '^b', ms=themarkersize)
#plt.plot(Mstar, MinCold/MinHalo, '^m', ms=themarkersize)
#plt.xticks(nums, labels, rotation='vertical')
plt.margins(0.2)
plt.yscale('linear')
plt.xscale('log')
#plt.ylim([10**1, 10**9.5])
#plt.ylim([0,1.0])
plt.xlim([10**6, 10**11])
plt.ylim([0,1.14])
plt.xlabel(r'Stellar Mass ($M_\odot$)', fontsize=titlefont)
plt.ylabel('fraction of halo CGM metals', fontsize=titlefont)
#plt.yscale('log')

ax = plt.gca()
#ax2 = ax.twiny()

ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax.yaxis.set_tick_params(width=tickwidth, length=ticklength)


thar22b =['Hot Phase (logT>5.3)','Warm Ions (4.7<logT<5.3)', 'Low Ions (4.0<logT<4.7)', 'Cold (logT<4.0)']

#plt.legend(thar22b, loc = 'best',  prop={'size':12})
plt.legend(thar22b, loc = 'upper left',  prop={'size':legsize}, ncol=2)
#ax = plt.gca()

plt.savefig(today + 'Mstar_halobars'+'.pdf')



missing_metals = Mexpected - totalMeta 
print 'missing metals ',missing_metals

plt.clf()
width = 0.6
print nums[fartoon]
print labels[fartoon]
fig5b = plt.figure(figsize=(11.5,12.5))
p1 = plt.bar(nums, MinHot[fartoon]/MinHalo[fartoon], width, color='r', alpha=1)
p2 = plt.bar(nums, MinO6[fartoon]/MinHalo[fartoon], width, color='k',
             bottom=MinHot[fartoon]/MinHalo[fartoon], alpha=1)
p3 = plt.bar(nums, MinLowIon[fartoon]/MinHalo[fartoon], width, color='b',
             bottom=(MinHot[fartoon]/MinHalo[fartoon] + MinO6[fartoon]/MinHalo[fartoon]), alpha=1)
p4 = plt.bar(nums, MinCold[fartoon]/MinHalo[fartoon], width, color='m',
             bottom=(MinHot[fartoon]/MinHalo[fartoon] + MinO6[fartoon]/MinHalo[fartoon] + MinLowIon[fartoon]/MinHalo[fartoon]), alpha=1)
#Mstar[-1]-=1e10
 
#p1a = plt.bar(nums+width, CMinHot[fartoon]/CMinHalo[fartoon], width, color='r', alpha=1, fill = False, hatch="/", ecolor='r',edgecolor="r")
#p2a = plt.bar(nums+width, CMinO6[fartoon]/CMinHalo[fartoon], width, color='k', fill =False, hatch="/",
#             bottom=CMinHot[fartoon]/CMinHalo[fartoon], alpha=1, ecolor='k',edgecolor="k")
#p3a = plt.bar(nums+width, CMinLowIon[fartoon]/CMinHalo[fartoon], width, color='b', fill = False, hatch="/", 
#             bottom=(CMinHot[fartoon]/CMinHalo[fartoon] + CMinO6[fartoon]/CMinHalo[fartoon]), alpha=1, ecolor='b',edgecolor="b")
#p4a = plt.bar(nums+width, CMinCold[fartoon]/CMinHalo[fartoon], width, color='m', fill = False, hatch="/",
#             bottom=(CMinHot[fartoon]/CMinHalo[fartoon] + CMinO6[fartoon]/CMinHalo[fartoon] + CMinLowIon[fartoon]/CMinHalo[fartoon]), alpha=1, ecolor='m',edgecolor="m")

#p1a = plt.bar(nums+width, CMinHot[fartoon]/CMinHalo[fartoon], width, color='r', alpha=0.75, fill = True, hatch="/", ecolor='k',edgecolor="k")
#p2a = plt.bar(nums+width, CMinO6[fartoon]/CMinHalo[fartoon], width, color='g', fill =True, hatch="/",
#             bottom=CMinHot[fartoon]/CMinHalo[fartoon], alpha=0.75, ecolor='k',edgecolor="k")
#p3a = plt.bar(nums+width, CMinLowIon[fartoon]/CMinHalo[fartoon], width, color='b', fill = True, hatch="/", 
#             bottom=(CMinHot[fartoon]/CMinHalo[fartoon] + CMinO6[fartoon]/CMinHalo[fartoon]), alpha=0.75, ecolor='k',edgecolor="k")
#p4a = plt.bar(nums+width, CMinCold[fartoon]/CMinHalo[fartoon], width, color='m', fill = True, hatch="/",
#             bottom=(CMinHot[fartoon]/CMinHalo[fartoon] + CMinO6[fartoon]/CMinHalo[fartoon] + CMinLowIon[fartoon]/CMinHalo[fartoon]), alpha=0.75, ecolor='k',edgecolor="k")

            
#plt.plot(Mstar, MinHot/MinHalo, '^r', ms=themarkersize)
#plt.plot(Mstar, MinO6/MinHalo, '^k', ms=themarkersize)
#plt.plot(Mstar, MinLowIon/MinHalo, '^b', ms=themarkersize)
#plt.plot(Mstar, MinCold/MinHalo, '^m', ms=themarkersize)
plt.xticks(nums+(width/2.0), labels[fartoon], rotation='vertical')
plt.margins(0.05)
plt.yscale('linear')
#plt.xscale('linear')
#plt.xscale('log')
#plt.ylim([10**1, 10**9.5])
#plt.ylim([0,1.0])
#plt.xlim([10**6, 10**11])
plt.ylim([0,1.14])
#plt.xlabel('Stellar Mass at z=0 (Msun)', fontsize=20)
plt.ylabel('fraction of halo CGM metals', fontsize=titlefont)
#plt.yscale('log')

thar22b =['Hot Phase (logT>5.3)','Warm Ions (4.7<logT<5.3)', 'Low Ions (4.0<logT<4.7)', 'Cold (logT<4.0)']
thar22b =[r'Hot Phase(T>$10^{5.3}$K)',r'Warm Ions ($10^{4.7}$K<T<$10^{5.3}$K)', r'Low Ions ($10^{4.0}$K<T<$10^{4.7}$K)', r'Cold Phase (T<$10^{4.0}$K)']

#plt.legend(thar22b, loc = 'best',  prop={'size':12})
plt.legend(thar22b, loc = 'best',  prop={'size':legsize}, ncol=2, frameon=False)


ax = plt.gca()
ax2 = ax.twiny()

ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax.yaxis.set_tick_params(width=tickwidth, length=ticklength)

ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(num_wanted)
ax2.set_xticklabels(Mstar_wanted_str)
ax2.set_xlabel(r'$M_* (M_{\odot})$', fontsize=titlefont)
ticklabels = ax2.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax2.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ax2.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax2.yaxis.set_tick_params(width=tickwidth, length=ticklength)   


#ax = plt.gca()

plt.savefig(today + 'label_halobars'+'.pdf')








plt.clf()
width = 0.6
fig5b = plt.figure(figsize=(11.5,12.5))
p1 = plt.bar(nums, MinHot[fartoon]/MinHalo[fartoon], width, color='r', alpha=1)
p2 = plt.bar(nums, MinO6[fartoon]/MinHalo[fartoon], width, color='k',
             bottom=MinHot[fartoon]/MinHalo[fartoon], alpha=1)
p3 = plt.bar(nums, MinLowIon[fartoon]/MinHalo[fartoon], width, color='b',
             bottom=(MinHot[fartoon]/MinHalo[fartoon] + MinO6[fartoon]/MinHalo[fartoon]), alpha=1)
p4 = plt.bar(nums, MinCold[fartoon]/MinHalo[fartoon], width, color='m',
             bottom=(MinHot[fartoon]/MinHalo[fartoon] + MinO6[fartoon]/MinHalo[fartoon] + MinLowIon[fartoon]/MinHalo[fartoon]), alpha=1)
#Mstar[-1]-=1e10
 
#p1a = plt.bar(nums+width, CMinHot[fartoon]/CMinHalo[fartoon], width, color='r', alpha=1, fill = False, hatch="/", ecolor='r',edgecolor="r")
#p2a = plt.bar(nums+width, CMinO6[fartoon]/CMinHalo[fartoon], width, color='k', fill =False, hatch="/",
#             bottom=CMinHot[fartoon]/CMinHalo[fartoon], alpha=1, ecolor='k',edgecolor="k")
#p3a = plt.bar(nums+width, CMinLowIon[fartoon]/CMinHalo[fartoon], width, color='b', fill = False, hatch="/", 
#             bottom=(CMinHot[fartoon]/CMinHalo[fartoon] + CMinO6[fartoon]/CMinHalo[fartoon]), alpha=1, ecolor='b',edgecolor="b")
#p4a = plt.bar(nums+width, CMinCold[fartoon]/CMinHalo[fartoon], width, color='m', fill = False, hatch="/",
#             bottom=(CMinHot[fartoon]/CMinHalo[fartoon] + CMinO6[fartoon]/CMinHalo[fartoon] + CMinLowIon[fartoon]/CMinHalo[fartoon]), alpha=1, ecolor='m',edgecolor="m")

#p1a = plt.bar(nums+width, CMinHot[fartoon]/CMinHalo[fartoon], width, color='r', alpha=0.75, fill = True, hatch="/", ecolor='k',edgecolor="k")
#p2a = plt.bar(nums+width, CMinO6[fartoon]/CMinHalo[fartoon], width, color='g', fill =True, hatch="/",
#             bottom=CMinHot[fartoon]/CMinHalo[fartoon], alpha=0.75, ecolor='k',edgecolor="k")
#p3a = plt.bar(nums+width, CMinLowIon[fartoon]/CMinHalo[fartoon], width, color='b', fill = True, hatch="/", 
#             bottom=(CMinHot[fartoon]/CMinHalo[fartoon] + CMinO6[fartoon]/CMinHalo[fartoon]), alpha=0.75, ecolor='k',edgecolor="k")
#p4a = plt.bar(nums+width, CMinCold[fartoon]/CMinHalo[fartoon], width, color='m', fill = True, hatch="/",
#             bottom=(CMinHot[fartoon]/CMinHalo[fartoon] + CMinO6[fartoon]/CMinHalo[fartoon] + CMinLowIon[fartoon]/CMinHalo[fartoon]), alpha=0.75, ecolor='k',edgecolor="k")

            
#plt.plot(Mstar, MinHot/MinHalo, '^r', ms=themarkersize)
#plt.plot(Mstar, MinO6/MinHalo, '^k', ms=themarkersize)
#plt.plot(Mstar, MinLowIon/MinHalo, '^b', ms=themarkersize)
#plt.plot(Mstar, MinCold/MinHalo, '^m', ms=themarkersize)
plt.xticks(nums+(width/2.0), Mstar_label_ar[fartoon], rotation='vertical')
plt.margins(0.05)
plt.yscale('linear')
#plt.xscale('linear')
#plt.xscale('log')
#plt.ylim([10**1, 10**9.5])
#plt.ylim([0,1.0])
#plt.xlim([10**6, 10**11])
plt.ylim([0,1.14])
#plt.xlabel('Stellar Mass at z=0 (Msun)', fontsize=20)
plt.ylabel('fraction of CGM metals', fontsize=titlefont)
plt.xlabel(r'${\rm log_{10}} (M_* /M_{\odot})$', fontsize=titlefont)

#plt.yscale('log')

thar22b =['Hot Phase (logT>5.3)','Warm Ions (4.7<logT<5.3)', 'Low Ions (4.0<logT<4.7)', 'Cold (logT<4.0)']
thar22b =[r'Hot Phase(T>$10^{5.3}$K)',r'Warm Ions ($10^{4.7}$K<T<$10^{5.3}$K)', r'Low Ions ($10^{4.0}$K<T<$10^{4.7}$K)', r'Cold Phase (T<$10^{4.0}$K)']

#plt.legend(thar22b, loc = 'best',  prop={'size':12})
plt.legend(thar22b, loc = 'best',  prop={'size':legsize}, ncol=2, frameon=False)


ax = plt.gca()
ax2 = ax.twiny()

ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax.yaxis.set_tick_params(width=tickwidth, length=ticklength)

ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(nums+width/2.0)
ax2.set_xticklabels(labels[fartoon], rotation='vertical')
#ax2.set_xlabel(r'$M_* (M_{\odot})$', fontsize=titlefont)
ticklabels = ax2.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax2.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ax2.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax2.yaxis.set_tick_params(width=tickwidth, length=ticklength)   
plt.tight_layout()


#ax = plt.gca()

plt.savefig(today + 'label_massbars'+'.pdf')









fig12 = plt.figure(figsize=(10,8))
plt.plot(Mstar, missing_metals, '^r', ms=themarkersize)

#plt.xticks(nums, labels, rotation='vertical')
plt.margins(0.1)
plt.yscale('log')
plt.xscale('log')
#plt.ylim([10**1, 10**9.5])
#plt.ylim([1e-3,1.0])
plt.xlim([10**6, 10**11])
plt.xlabel(r'Stellar Mass at z=0 ($M_\odot$)', fontsize=titlefont)
plt.ylabel('mass of missing metals', fontsize=titlefont)
plt.yscale('log')

ax = plt.gca()

ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax.yaxis.set_tick_params(width=tickwidth, length=ticklength)


#plt.legend(thar22b, loc = 'best',  prop={'size':12})
ax = plt.gca()

plt.savefig(today + 'missingmetals'+'.pdf')

missing_metals_frac = missing_metals / Mexpected

fig13 = plt.figure(figsize=(10,8))
plt.plot(Mstar, missing_metals_frac, 'sb', ms=themarkersize)

#plt.xticks(nums, labels, rotation='vertical')
plt.margins(0.2)
plt.yscale('linear')
plt.xscale('log')
#plt.ylim([10**1, 10**9.5])
plt.ylim([0,1.0])
plt.xlim([10**6, 10**11])
plt.xlabel(r'Stellar Mass at z=0 ($M_\odot$)', fontsize=titlefont)
plt.ylabel('frac mass of missing metals', fontsize=titlefont)
#plt.yscale('log')

ax = plt.gca()

ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(labelfont)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=tickwidth, length=ticklength)
ax.yaxis.set_tick_params(width=tickwidth, length=ticklength)


#plt.legend(thar22b, loc = 'best',  prop={'size':12})
ax = plt.gca()

plt.savefig(today + 'fracmissingmetals'+'.pdf')

#MissingMet_Outfit = optimize.leastsq(residuals,p2[:], args=(metrat_ar_ext[stoopcut * fillcut_ar_ext * garthcut * hiz_cut], ISMmet_ar_ext[stoopcut * fillcut_ar_ext * garthcut * hiz_cut] ))

#ISMmet_int = ISM_Outfit[0][0]
#ISMmet_slope = ISM_Outfit[0][1]
#print 'heres ISM metrat fit ',ISMmet_int,ISMmet_slope


#p1 = plt.bar(ind, menMeans, width, color='r', yerr=menStd)
#p2 = plt.bar(ind, womenMeans, width, color='y',
#             bottom=menMeans, yerr=womenStd)
