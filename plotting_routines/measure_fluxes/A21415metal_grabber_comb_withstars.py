#only works for 91814 and beyond
#does not work right for 92414 and beyond! lets see why

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from datetime import date
import sys
import scipy.optimize as optimize

today = date.today()
print today

#these parameters should be set by command prompt
zstart =999
zend = 999
ShellUseInner = '2'
ShellUseOuter = '9'
InnerISM_only = True #use only Metals within 0.1 Rvir as ISM (otherwise 0.2 Rvir) - almost irrelevant now that i have this other mode for ISM anyway
ISMMass_Is_01Rvir = True #false if you want physical val
ISMMass_physical_border = 10.0 # physical kpc
fracwewant = 0.8
usev0cut = True
lozmedzcomb = False
naked = True
jacket = False
usesmoothe = True
specialshell = True #use specialval Extra Shell instead of whatever input was
doublespecial = True
specialval = 2 #for 25 kpc, specialval = 1. for 1.0 Rvir, specialval = 5
Dspecialval = 1  #0 for 10 kpc, 1 for 25 kpc, 2 for 50 kpc - confusing but this is inner region
speciallabel = '50 kpc'
doublespeciallabel = '25 kpc'
#speciallabel = '25 kpc'
#doublespeciallabel = '10 kpc'
barfolemew = True
#speciallabel = r'1.0 $R_{vir}$'

SFRplotboost = 1
Metalplotboost = 100

SFR_correct_factor = 1.15  #to correct for stars formed between snapshots 
little_h = 0.7 
legsize = 16
Smoothing_interval=400.0
smooth_half = Smoothing_interval/2.0


def smoothe_over400(iT, iquant):
	returnquant = []
	count = 0 
	while (count < len(iT)):
		iTmin = iT[count]-smooth_half
		iTmax = iT[count]+smooth_half
		iTcut = (iT > iTmin) * (iT < iTmax)
#		if (len(iquant[iTcut*np.isreal(iquant)*(np.fabs(iquant)>1e-10)]) > 0):
		if (len(iquant[iTcut*np.isreal(iquant)]) > 0):
			iTempQuant = np.mean(iquant[iTcut*np.isreal(iquant)])
			returnquant.append(iTempQuant)
		else:
			returnquant.append(0)
		count+=1
		
	return np.array(returnquant)

if (len(sys.argv) < 2):
	print 'syntax  blah.py Nhalo directory date zstart zend use_v0 name starhistory_finname'
	print 'dir: outflows_91814/  date: 2014-09-24'
	sys.exit()
haloN = str(sys.argv[1])
thedir = str(sys.argv[2])
thedate = str(sys.argv[3])
zstart = float(sys.argv[4])
zend = float(sys.argv[5])
use_v0_cut = str(sys.argv[6])
name = str(sys.argv[7])
starhistory = str(sys.argv[8])

if (use_v0_cut == 'y'):
	usev0cut = True
else:
	usev0cut = False
print 'got zstart, zend', zstart, zend

filename = name+'/'+thedir + thedate + 'OutflowN'+haloN+'.txt'
filename2 =name+'/'+ thedir + thedate + 'B_OutflowN'+haloN+'.txt' 
filename3 =name+'/'+ thedir + thedate + 'InflowN'+haloN+'.txt'
filename4 = name+'/'+thedir + thedate + 'B_InflowN'+haloN+'.txt'
if (usev0cut): 
	filename =name+'/'+ thedir + thedate + '0_OutflowN'+haloN+'.txt'
	filename2 =name+'/'+ thedir + thedate + '0B_OutflowN'+haloN+'.txt' 
	filename3 =name+'/'+ thedir + thedate + '0_InflowN'+haloN+'.txt'
	filename4 =name+'/'+ thedir + thedate + '0B_InflowN'+haloN+'.txt'

if (barfolemew):
	filename5 = name+'/'+thedir + thedate + 'OutflowN'+haloN+'.txt'
	filename6 =name+'/'+ thedir + thedate + 'B_OutflowN'+haloN+'.txt' 



starf = open(starhistory)
stardars = np.loadtxt(starf, ndmin=2)

amed = stardars[:,0] + stardars[:,1]
amed /= 2.0
starmetalbins = stardars[:,2]
starmassbins = stardars[:,3]
starmetalcum = np.cumsum(starmetalbins) * 1e10 / little_h


def line_to_string(line, newline = 'y'):
	linestring = ''
	for q in line:
		linestring += "{0:.6g}".format(q) + '  '
	if (newline == 'y'):
		linestring+=' \n'
	return linestring


if (zstart > 3.5):
	z_wanted = np.array([4.0, 3.5, 3.0, 2.5, 2.0])
	zmedpoint = 3.0
if (zstart > 1 and zstart < 3.5):
	z_wanted = np.array([2.0, 1.5, 1.25, 1.0, 0.75, 0.5])
	zmedpoint = 1.25
if (zstart < 1):
	z_wanted = np.array([0.5, 0.25, 0.1, 0.0])
	zmedpoint = 0.25

#z_wanted = np.array([4.0, 2.0, 1.0, 0.5, 0.25, 0.0])
#zmedpoint = 1.0

z_wanted_str = []
for q in z_wanted:
	z_wanted_str.append(str(q))
a_wanted = 1.0 / (1 + z_wanted)

finname = name+'/'+'tasz.txt'
f = open(finname)
dars = np.loadtxt(f, ndmin=2)
theas = dars[:,0]
thets = dars[:,1]
f.close()

finname = filename
f = open(finname)
dars = np.loadtxt(f, ndmin=2)
Zs = dars[:,1]
cut1 = Zs >= zend
cut2 = Zs <= zstart
cuts = cut1*cut2

Zs = Zs[cuts] 
zend = max(zend, np.min(Zs))

As = 1.0 / (1.0 + Zs)

tfora = interp1d(theas, thets)
afort = interp1d(thets, theas)
Ts = tfora(As)

starTs = tfora(amed) * 1000

Ts *= 1000.0

T_wanted = tfora(a_wanted)
T_wanted *= 1000

Rvirused = dars[:,2][cuts]
SFRInner = dars[:,3][cuts]
SFRouter = dars[:,4][cuts]
NewSFR_timescale = dars[:,11][cuts]

#SFR1 = SFRInner + SFRouter
#for now excluding SF outside 0.2 Rvir
SFR1 = SFRInner * SFR_correct_factor

NumShells = 20
Flow_Shells = []

#This 'first swipe' only gets outflow rates (no metals)
starter = -5 
count = 0
while (count < NumShells):
	anticount = NumShells - count - 1
	riteindex = -6*anticount+starter
	TheFlow = dars[:,riteindex][cuts]
	#print 'The Flow of Shell ', count,' is index ',riteindex
	Flow_Shells.append(TheFlow)
	count += 1

count = 0
NumExtraShells = 6
Extra_FlowShells = []
while (count < (NumShells + NumExtraShells)):
	anticount = NumShells + NumExtraShells - count - 1
	riteindex = -6*anticount+starter
	TheFlow= dars[:,riteindex][cuts]
	Extra_FlowShells.append(TheFlow)
	count += 1


NPartCross_ar = []
NPartFlux_ar = []
NpartShell_ar = []
starter1 = -2
starter2 = -3
starter3 = -4 
count = 0

while (count <= NumShells):
	anticount = NumShells - count - 1
	riteindex = -6*anticount+starter1
	riteindex2 = -6*anticount+starter2
	riteindex3 = -6*anticount+starter3
	NpartInOutflow_Cross = dars[:,riteindex][cuts]
	NpartInOutflow_Flux = dars[:,riteindex2][cuts]
	NpartInShell = dars[:,riteindex3][cuts]
	
	NPartCross_ar.append(NpartInOutflow_Cross)
	NPartFlux_ar.append(NpartInOutflow_Flux)
	NpartShell_ar.append(NpartInShell)
	#print 'The Flow of Shell ', count,' is index ',riteindex
	count += 1

#adding this part on 2-11-14
ExtraNPartCross_ar = []
ExtraNPartFlux_ar = []
ExtraNpartShell_ar = []
count = 0
NumExtraShells = 6
while (count < (NumShells + NumExtraShells)):
	anticount = NumShells + NumExtraShells - count - 1
	riteindex = -6*anticount+starter1
	riteindex2 = -6*anticount+starter2
	riteindex3 = -6*anticount+starter3
	NpartInOutflow_Cross = dars[:,riteindex][cuts]
	NpartInOutflow_Flux = dars[:,riteindex2][cuts]
	NpartInShell = dars[:,riteindex3][cuts]
	
	ExtraNPartCross_ar.append(NpartInOutflow_Cross)
	ExtraNPartFlux_ar.append(NpartInOutflow_Flux)
	ExtraNpartShell_ar.append(NpartInShell)
	count+=1
	#print 'The Flow of Shell ', count,' is index ',riteindex


NPartCross_ar=np.array(NPartCross_ar)
NPartFlux_ar=np.array(NPartFlux_ar)
NpartShell_ar=np.array(NpartShell_ar)
starter = -1 
count = 0
Cross_Shells = []
while (count < NumShells):
	anticount = NumShells - count - 1
	riteindex = -6*anticount+starter
	TheCross = (dars[:,riteindex][cuts]* 1e10 / little_h) / NewSFR_timescale
	#print 'The Cross of Shell ', count,' is index ',riteindex
	Cross_Shells.append(TheCross)
	count+=1

Extra_CrossShells = []
while (count < (NumShells + NumExtraShells)):
	anticount = NumShells + NumExtraShells - count - 1
	riteindex = -6*anticount+starter
	TheCross= (dars[:,riteindex][cuts]* 1e10 / little_h) / NewSFR_timescale
	Extra_CrossShells.append(TheCross)
	count += 1


#time to get the outflow metals - from file 2  B_outflow
finname = filename2
f = open(finname)
dars2 = np.loadtxt(f, ndmin=2)

count = 0
starter1 = -3
starter2 = -4
starter3 = -2
starter4 = -1

MetalFlowShells = []
MetalFlowMassShells = []
ISMShellMetal = []
ISMShellMass = []
while (count < NumShells):
	anticount = NumShells - count - 1
	riteindex = -5*anticount+starter1
	riteindex2 = -5*anticount+starter2
	riteindex3 = -5*anticount+starter3
	riteindex4 = -5*anticount+starter4

	MetalFlow = dars2[:,riteindex][cuts]
	MassMetalOutflow= dars2[:,riteindex2][cuts]
	TotalShellMetal = dars2[:,riteindex3][cuts]
	TotalShellMass = dars2[:,riteindex4][cuts]
	#print 'The Flow of Shell ', count,' is index ',riteindex
	MetalFlowShells.append(MetalFlow)	
	MetalFlowMassShells.append(MassMetalOutflow)
	ISMShellMetal.append(TotalShellMetal)
	ISMShellMass.append(TotalShellMass)
	count += 1
ExtraMetalFlowShells = []
ExtraMetalMassFlowShells = []
ExtraISMShellMetal = []
ExtraISMShellMass = []
count = 0
while (count < NumShells + NumExtraShells):
	anticount = NumShells + NumExtraShells - count - 1
	riteindex = -5*anticount+starter1
	riteindex2 = -5*anticount+starter2
	riteindex3 = -5*anticount+starter3
	riteindex4 = -5*anticount+starter4
	MetalFlow = dars2[:,riteindex][cuts]
	MassMetalOutflow= dars2[:,riteindex2][cuts]
	TotalShellMetal = dars2[:,riteindex3][cuts]
	TotalShellMass = dars2[:,riteindex4][cuts]
	#print 'The Flow of Shell ', count,' is index ',riteindex
	ExtraMetalFlowShells.append(MetalFlow)	
	ExtraMetalMassFlowShells.append(MassMetalOutflow)
	ExtraISMShellMetal.append(TotalShellMetal)
	ExtraISMShellMass.append(TotalShellMass)
	count += 1

if (ISMMass_Is_01Rvir):
	ISMMass = np.copy(ISMShellMetal[0])
	ISM_TMass = np.copy(ISMShellMass[0])
else:
	count = 0
	gatherISM = []
	gatherMass = []
	while (count < len(Rvirused)):
		physical_boundary  = ISMMass_physical_border / (Rvirused[count] * As[count] / little_h) #gives fraction of Rvir that should be used
		print 'physical boundary ',physical_boundary
		tempISMMass = 0 
		tempTMass = 0
		lastbin = int(np.floor(physical_boundary*10.0))
		if (lastbin > (NumShells-1)):
			print 'error too big '
			lastbin = NumShells-1
		countb = 0
		while (countb < lastbin):
			tempISMMass += ISMShellMetal[countb][count]
			tempTMass += ISMShellMass[countb][count]
			countb+=1
		remainderfrac = physical_boundary*10 - np.floor(physical_boundary*10)
		print 'remainder frac ',remainderfrac
		tempISMMass += ISMShellMetal[lastbin][count]*remainderfrac
		tempTMass += ISMShellMass[lastbin][count]*remainderfrac
		gatherISM.append(tempISMMass)
		gatherMass.append(tempTMass)
		count+=1
	ISMMass = np.array(gatherISM)
	ISM_TMass = np.array(gatherMass)
	
	
	

ISMMass *= 1e10
ISM_TMass *= 1e10

#now getting regular inflow stuff
finname = filename3
f = open(finname)
dars3 = np.loadtxt(f, ndmin=2)

InFlow_Shells = []

starter = -5 
count = 0
while (count < NumShells):
	anticount = NumShells - count  -1 
	riteindex = -6*anticount+starter
	TheFlow = dars3[:,riteindex][cuts]
	InFlow_Shells.append(TheFlow)
	count +=1

count = 0
Extra_InFlowShells = []
while (count < (NumShells + NumExtraShells)):
	anticount = NumShells + NumExtraShells - count - 1
	riteindex = -6*anticount+starter
	TheFlow= dars3[:,riteindex][cuts]
	Extra_InFlowShells.append(TheFlow)
	count += 1



starter = -1 
count = 0
InCross_Shells = []
while (count < NumShells):
	anticount = NumShells - count  -1
	riteindex = -6*anticount+starter
	TheCross = (dars3[:,riteindex][cuts]* 1e10 / little_h) / NewSFR_timescale * -1 
	#print 'The Cross of Shell ', count,' is index ',riteindex
	InCross_Shells.append(TheCross)
	count+=1
	
count = 0
Extra_InCrossShells = []
while (count < (NumShells + NumExtraShells)):
	anticount = NumShells + NumExtraShells - count - 1
	riteindex = -6*anticount+starter
	TheCross= (dars3[:,riteindex][cuts]* 1e10 / little_h) / NewSFR_timescale
	Extra_InCrossShells.append(TheCross)
	count += 1


finname = filename4
f = open(finname)
dars4 = np.loadtxt(f, ndmin=2)

count = 0
starter1 = -1
starter2 = -2

InMetalFlowShells = []
InMetalMassFlowShells = []
while (count < NumShells):
	anticount = NumShells - count - 1
	riteindex = -3*anticount+starter1
	riteindex2 = -3*anticount+starter2
	InMetalFlow = dars4[:,riteindex][cuts]
	InMassMetalOutflow= dars4[:,riteindex2][cuts]
	InMetalFlowShells.append(InMetalFlow)	
	InMetalMassFlowShells.append(InMassMetalOutflow)
	count += 1

count = 0
ExtraInMetalFlowShells = []
ExtraInMetalMassFlowShells = []
while (count < NumShells + NumExtraShells):
	anticount = NumShells + NumExtraShells - count - 1
	riteindex = -3*anticount+starter1
	riteindex2 = -3*anticount+starter2
	InMetalFlow = dars4[:,riteindex][cuts]
	InMassMetalOutflow= dars4[:,riteindex2][cuts]
	ExtraInMetalFlowShells.append(InMetalFlow)	
	ExtraInMetalMassFlowShells.append(InMassMetalOutflow)
	count += 1

	
MetalFlowShells = np.array(MetalFlowShells)
Flow_Shells = np.array(Flow_Shells)
InMetalFlowShells = np.array(InMetalFlowShells)
InFlow_Shells = np.array(InFlow_Shells)


if (barfolemew):
	finname = filename5
	f = open(finname)
	dars5 = np.loadtxt(f, ndmin=2)
	barfFlow_Shells = []

	#This 'first swipe' only gets outflow rates (no metals)
	starter = -5 
	count = 0
	while (count < NumShells):
		anticount = NumShells - count - 1
		riteindex = -6*anticount+starter
		TheFlow = dars5[:,riteindex][cuts]
		#print 'The Flow of Shell ', count,' is index ',riteindex
		barfFlow_Shells.append(TheFlow)
		count += 1

	count = 0
	barfExtra_FlowShells = []
	while (count < (NumShells + NumExtraShells)):
		anticount = NumShells + NumExtraShells - count - 1
		riteindex = -6*anticount+starter
		TheFlow= dars5[:,riteindex][cuts]
		barfExtra_FlowShells.append(TheFlow)
		count += 1



	#adding this part on 2-11-14

	#time to get the outflow metals - from file 2  B_outflow
	finname = filename6
	f = open(finname)
	dars6 = np.loadtxt(f, ndmin=2)

	count = 0
	starter1 = -3
	starter2 = -4
	starter3 = -2
	starter4 = -1

	barfMetalFlowShells = []
	barfMetalFlowMassShells = []
	barfISMShellMetal = []
	barfISMShellMass = []
	while (count < NumShells):
		anticount = NumShells - count - 1
		riteindex = -5*anticount+starter1
		riteindex2 = -5*anticount+starter2
		riteindex3 = -5*anticount+starter3
		riteindex4 = -5*anticount+starter4

		MetalFlow = dars6[:,riteindex][cuts]
		MassMetalOutflow= dars6[:,riteindex2][cuts]
		TotalShellMetal = dars6[:,riteindex3][cuts]
		TotalShellMass = dars6[:,riteindex4][cuts]
		#print 'The Flow of Shell ', count,' is index ',riteindex
		barfMetalFlowShells.append(MetalFlow)	
		barfMetalFlowMassShells.append(MassMetalOutflow)
		barfISMShellMetal.append(TotalShellMetal)
		barfISMShellMass.append(TotalShellMass)
		count += 1
	barfExtraMetalFlowShells = []
	barfExtraMetalMassFlowShells = []
	barfExtraISMShellMetal = []
	barfExtraISMShellMass = []
	count = 0
	while (count < NumShells + NumExtraShells):
		anticount = NumShells + NumExtraShells - count - 1
		riteindex = -5*anticount+starter1
		riteindex2 = -5*anticount+starter2
		riteindex3 = -5*anticount+starter3
		riteindex4 = -5*anticount+starter4
		MetalFlow = dars6[:,riteindex][cuts]
		MassMetalOutflow= dars6[:,riteindex2][cuts]
		TotalShellMetal = dars6[:,riteindex3][cuts]
		TotalShellMass = dars6[:,riteindex4][cuts]
		#print 'The Flow of Shell ', count,' is index ',riteindex
		barfExtraMetalFlowShells.append(MetalFlow)	
		barfExtraMetalMassFlowShells.append(MassMetalOutflow)
		barfExtraISMShellMetal.append(TotalShellMetal)
		barfExtraISMShellMass.append(TotalShellMass)
		count += 1



Flow_in = Flow_Shells[int(ShellUseInner)]
Flow_out = Flow_Shells[int(ShellUseOuter)]

Cross_in = Cross_Shells[int(ShellUseInner)]
Cross_out = Cross_Shells[int(ShellUseOuter)]

Flow_in2 = InFlow_Shells[int(ShellUseInner)]
Flow_out2 = InFlow_Shells[int(ShellUseOuter)]

Cross_in2 = InCross_Shells[int(ShellUseInner)] 
Cross_out2 = InCross_Shells[int(ShellUseOuter)]

MetalFlow_in = MetalFlowShells[int(ShellUseInner)]
MetalFlow_out = MetalFlowShells[int(ShellUseOuter)]

MetalFlow_in2 = InMetalFlowShells[int(ShellUseInner)]
MetalFlow_out2 = InMetalFlowShells[int(ShellUseOuter)]

if (barfolemew):
	barfFlow_in = barfFlow_Shells[int(ShellUseInner)]
	barfFlow_out = barfFlow_Shells[int(ShellUseOuter)]
	barfMetalFlow_in = barfMetalFlowShells[int(ShellUseInner)]
	barfMetalFlow_out = barfMetalFlowShells[int(ShellUseOuter)]	


if (specialshell):
	Flow_out = Extra_FlowShells[int(specialval)]
	Cross_out = Extra_CrossShells[int(specialval)]
	Flow_out2 = Extra_InFlowShells[int(specialval)]
	Cross_out2 = Extra_InCrossShells[int(specialval)]
	MetalFlow_out = ExtraMetalFlowShells[int(specialval)]
	MetalFlow_out2 = ExtraInMetalFlowShells[int(specialval)]
	ShellUseOuter = 'S'+str(specialval)
	if (barfolemew):
		barfFlow_out = barfExtra_FlowShells[int(specialval)]
		barfMetalFlow_out = barfExtraMetalFlowShells[int(specialval)]	

		
if (doublespecial):
	Flow_in = Extra_FlowShells[int(Dspecialval)]
	Cross_in = Extra_CrossShells[int(Dspecialval)]	
	Flow_in2 = Extra_InFlowShells[int(Dspecialval)]
	Cross_in2 = Extra_InCrossShells[int(Dspecialval)]
	MetalFlow_in = ExtraMetalFlowShells[int(Dspecialval)]
	MetalFlow_in2 = ExtraInMetalFlowShells[int(Dspecialval)]
	ShellUseInner = 'S'+str(Dspecialval)
	if (barfolemew):
		barfFlow_in = barfExtra_FlowShells[int(Dspecialval)]
		barfMetalFlow_in = barfExtraMetalFlowShells[int(Dspecialval)]	

OutflowMetalFrac2 = MetalFlow_in / Flow_in
OutflowMetalFrac9  = MetalFlow_out / Flow_out
	

Zsolar = 0.02

if (barfolemew):
	barfOutflowMetalFrac2 = (barfMetalFlow_in / barfFlow_in) / Zsolar
	barfOutflowMetalFrac9 = (barfMetalFlow_out / barfFlow_out) / Zsolar	


OutflowMetalFrac2 /= Zsolar
OutflowMetalFrac9 /= Zsolar

InflowMetalFrac2 = MetalFlow_in2 / Flow_in2
InflowMetalFrac9  = MetalFlow_out2 / Flow_out2

InflowMetalFrac2 /= Zsolar 
InflowMetalFrac9 /= Zsolar 

ISM_metal_fraction = (ISMShellMetal[0] + ISMShellMetal[1]) / (ISMShellMass[0] + ISMShellMass[1]) / Zsolar


if (InnerISM_only):
	#ISM_metal_fraction = (ISMShellMetal[0] / ISMShellMass[0]) / Zsolar
	#print 'now what the fuck ',ISM_metal_fraction
	ISM_metal_fraction = ISMMass / ISM_TMass / Zsolar
	print 'now what the fuck ',ISM_metal_fraction

Halo_metal_FS = ( ISMShellMetal[1] + ISMShellMetal[2] + ISMShellMetal[3] + ISMShellMetal[4] + ISMShellMetal[5] + ISMShellMetal[6] + ISMShellMetal[7] + ISMShellMetal[8] + ISMShellMetal[9]) * 1e10

Halo_mass_FS = ( ISMShellMass[1] + ISMShellMass[2] + ISMShellMass[3] + ISMShellMass[4] + ISMShellMass[5] + ISMShellMass[6] + ISMShellMass[7] + ISMShellMass[8] + ISMShellMass[9]) * 1e10

Halo_metal_fraction = (Halo_metal_FS/ Halo_mass_FS ) /Zsolar

benchmark_mass = np.median((ISMShellMass[0] + ISMShellMass[1]))

print 'benchmark !',benchmark_mass

metallicity_cut = (ISMShellMass[0] + ISMShellMass[1]) > 0.1*benchmark_mass

benchmark_mass = np.median((ISMShellMass[0]))

metallicity_cut = (ISMShellMass[0]) > 0.1*benchmark_mass


fig1 = plt.figure(figsize=(10,8))

OutFrac2cut = MetalFlow_in > 0.1*np.median(MetalFlow_in)
OutFrac9cut = MetalFlow_out > 0.1*np.median(MetalFlow_out)
InFrac2cut = MetalFlow_in2 < 0.1 * np.median(MetalFlow_in2)
InFrac9cut = MetalFlow_out2 < 0.1 * np.median(MetalFlow_out2)

BeefyOutflow2 = MetalFlow_in > 0.1*np.max(MetalFlow_in)
BeefyOutflow9 = MetalFlow_out > 0.1*np.max(MetalFlow_out)
BeefyInflow2 = MetalFlow_in2 < 0.1 * np.min(MetalFlow_in2)
BeefyInflow9 = MetalFlow_out2 < 0.1 * np.min(MetalFlow_out2)

realcut = np.isreal(ISM_metal_fraction)
poscut = ISM_metal_fraction > 1e-10 
realcut *= poscut

MetalMassFlow_in = NewSFR_timescale * MetalFlow_in
sortorder = np.argsort(MetalMassFlow_in)
TotalMetalMassFlow_in = np.sum(MetalMassFlow_in)
count = 0
cumTotal = 0
beeflist = []
print 'total we are trying to achieve is :',TotalMetalMassFlow_in
while (count < len(MetalMassFlow_in)):
	theind = -1 - count
	cumTotal += MetalMassFlow_in[sortorder][theind]
	print 'added z= ', Zs[sortorder][theind], MetalMassFlow_in[sortorder][theind]
	beeflist.append(Zs[sortorder][theind])
	print 'cum total so far is ',cumTotal
	if (cumTotal > fracwewant * TotalMetalMassFlow_in):
		print 'we have enough!'
		break
	count +=1
beeflist = np.array(beeflist)
Qualify_Flow2 = np.in1d(Zs, beeflist)
print 'f', Qualify_Flow2
print Zs[Qualify_Flow2]
#print 'meh ',Zs[BeefyOutflow2]

MetalMassFlow_out = NewSFR_timescale * MetalFlow_out
sortorder = np.argsort(MetalMassFlow_out)
TotalMetalMassFlow_out = np.sum(MetalMassFlow_out)
count = 0
cumTotal = 0
beeflist = []
print 'total we are trying to achieve is :',TotalMetalMassFlow_out
while (count < len(MetalMassFlow_out)):
	theind = -1 - count
	cumTotal += MetalMassFlow_out[sortorder][theind]
	print 'added z= ', Zs[sortorder][theind], MetalMassFlow_out[sortorder][theind]
	beeflist.append(Zs[sortorder][theind])
	print 'cum total so far is ',cumTotal
	if (cumTotal > fracwewant * TotalMetalMassFlow_out):
		print 'we have enough!'
		break
	count +=1
beeflist = np.array(beeflist)
Qualify_Flow9 = np.in1d(Zs, beeflist)
print 'f', Qualify_Flow9
print Zs[Qualify_Flow9]
#print 'meh ',Zs[BeefyOutflow2]


MetalMassFlow_in2 = NewSFR_timescale * MetalFlow_in2
sortorder = np.argsort(MetalMassFlow_in2)
TotalMetalMassFlow_in2 = np.sum(MetalMassFlow_in2)
count = 0
cumTotal = 0
beeflist = []
print 'total we are trying to achieve is :',TotalMetalMassFlow_in2
while (count < len(MetalMassFlow_in2)):
	theind = count
	cumTotal += MetalMassFlow_in2[sortorder][theind]
	print 'added z= ', Zs[sortorder][theind], MetalMassFlow_in2[sortorder][theind]
	beeflist.append(Zs[sortorder][theind])
	print 'cum total so far is ',cumTotal
	if (cumTotal < fracwewant * TotalMetalMassFlow_in2):
		print 'we have enough!'
		break
	count +=1
beeflist = np.array(beeflist)
Qualify_InFlow2 = np.in1d(Zs, beeflist)
print 'f', Qualify_InFlow2
print Zs[Qualify_InFlow2]

MetalMassFlow_out2 = NewSFR_timescale * MetalFlow_out2
sortorder = np.argsort(MetalMassFlow_out2)
TotalMetalMassFlow_out2 = np.sum(MetalMassFlow_out2)
count = 0
cumTotal = 0
beeflist = []
print 'total we are trying to achieve is :',TotalMetalMassFlow_out2
while (count < len(MetalMassFlow_out2)):
	theind = count
	cumTotal += MetalMassFlow_out2[sortorder][theind]
	print 'added z= ', Zs[sortorder][theind], MetalMassFlow_out2[sortorder][theind]
	beeflist.append(Zs[sortorder][theind])
	print 'cum total so far is ',cumTotal
	if (cumTotal < fracwewant * TotalMetalMassFlow_out2):
		print 'we have enough!'
		break
	count +=1
beeflist = np.array(beeflist)
Qualify_InFlow9 = np.in1d(Zs, beeflist)
print 'f', Qualify_InFlow9
print Zs[Qualify_InFlow9]

#in1, = plt.plot(Ts[OutFrac2cut], OutflowMetalFrac2[OutFrac2cut], color='r', ls='-', lw=2, alpha=0.2)
if (usesmoothe):
	OutflowMetalFrac2smooth = smoothe_over400(Ts[OutFrac2cut], OutflowMetalFrac2[OutFrac2cut])
	#OutflowMetalFrac2smooth = smoothe_over400(Ts, OutflowMetalFrac2)

	in1, = plt.plot(Ts[OutFrac2cut], OutflowMetalFrac2smooth, color='r', ls='-', lw=2, alpha=0.2)
	#in1, = plt.plot(Ts, OutflowMetalFrac2smooth, color='r', ls='-', lw=2, alpha=0.2)

else:
	in1, = plt.plot(Ts[OutFrac2cut], OutflowMetalFrac2[OutFrac2cut], color='r', ls='-', lw=2, alpha=0.2)

if (usesmoothe):
	OutflowMetalFrac9smooth = smoothe_over400(Ts[OutFrac9cut], OutflowMetalFrac9[OutFrac9cut])
	#OutflowMetalFrac9smooth = smoothe_over400(Ts, OutflowMetalFrac9)

	in2, = plt.plot(Ts[OutFrac9cut], OutflowMetalFrac9smooth, color='r',  ls=':', lw=2, alpha=0.2)
	#in2, = plt.plot(Ts, OutflowMetalFrac9smooth, color='r',  ls=':', lw=2, alpha=0.2)

else:
	in2, = plt.plot(Ts[OutFrac9cut], OutflowMetalFrac9[OutFrac9cut], color='r',  ls=':', lw=2, alpha=0.2)

#out1, = plt.plot(Ts[InFrac2cut], InflowMetalFrac2[InFrac2cut], color='b', ls='-', lw=2, alpha=0.2)
if (usesmoothe):
	InflowMetalFrac2smooth = smoothe_over400(Ts[InFrac2cut], InflowMetalFrac2[InFrac2cut])
	#InflowMetalFrac2smooth = smoothe_over400(Ts, InflowMetalFrac2)
	out1, = plt.plot(Ts[InFrac2cut], InflowMetalFrac2smooth, color='b', ls='-', lw=2, alpha=0.2)
	#out1, = plt.plot(Ts, InflowMetalFrac2smooth, color='b', ls='-', lw=2, alpha=0.2)

else:
	out1, = plt.plot(Ts[InFrac2cut], InflowMetalFrac2[InFrac2cut], color='b', ls='-', lw=2, alpha=0.2)


#out2, = plt.plot(Ts[InFrac9cut], InflowMetalFrac9[InFrac9cut], color='b', ls=':', lw=2, alpha=0.2)
if (usesmoothe):
	InflowMetalFrac9smooth = smoothe_over400(Ts[InFrac9cut], InflowMetalFrac9[InFrac9cut])
	#InflowMetalFrac9smooth = smoothe_over400(Ts, InflowMetalFrac9)
	#out2, = plt.plot(Ts, InflowMetalFrac9smooth, color='b', ls=':', lw=2, alpha=0.2)

	out2, = plt.plot(Ts[InFrac9cut], InflowMetalFrac9smooth, color='b', ls=':', lw=2, alpha=0.2)
else:
	out2, = plt.plot(Ts[InFrac9cut], InflowMetalFrac9[InFrac9cut], color='b', ls=':', lw=2, alpha=0.2)


benchmark = np.mean(ISMShellMass[0] + ISMShellMass[0])

#ism1, = plt.plot(Ts[metallicity_cut], ISM_metal_fraction[metallicity_cut], '.k', alpha=0.5, ms=3)



fakeT = np.copy(Ts)
fakeY = smoothe_over400(Ts,Halo_metal_fraction)
#badcut = np.invert(metallicity_cut*)
#fakeT[badcut] = np.nan
#fakeY[badcut] = np.nan
ism3, = plt.plot(fakeT, fakeY, '--k', alpha=1.0, lw=3, ms=3)


#ism3, = plt.plot(Ts[metallicity_cut], Halo_metal_fraction[metallicity_cut], '_k', alpha=0.5, ms=3)



#ism1, = plt.plot(Ts[realcut], ISM_metal_fraction[realcut], '.k', ms=6)
#ism3, = plt.plot(Ts[realcut], Halo_metal_fraction[realcut], 'xk', ms=3)
#ins1, = plt.plot(Ts[Qualify_Flow2], OutflowMetalFrac2[Qualify_Flow2], '^r', ms=5, fillstyle='full', alpha=0.5 )
#ins2, = plt.plot(Ts[Qualify_Flow9], OutflowMetalFrac9[Qualify_Flow9], 'or', ms=5, fillstyle='full', alpha=0.5 )


#ins1, = plt.plot(Ts[Qualify_Flow2], OutflowMetalFrac2[Qualify_Flow2], '-r', ms=5, fillstyle='full', alpha=0.5 )

if (usesmoothe):
	OutflowMetalFrac2smooth = smoothe_over400(Ts, OutflowMetalFrac2)
	
fakeT = np.copy(Ts[OutFrac2cut])
fakeY = np.copy(OutflowMetalFrac2)
if (usesmoothe):
	OutflowMetalFrac2smooth = smoothe_over400(Ts[OutFrac2cut], OutflowMetalFrac2[OutFrac2cut])
	fakeY = np.copy(OutflowMetalFrac2smooth)
fakeY2 = np.copy(OutflowMetalFrac2)
if (barfolemew):
	fakeY2 = np.copy(barfOutflowMetalFrac2[OutFrac2cut])
badcut = np.invert(Qualify_Flow2[OutFrac2cut])
fakeT[badcut] = np.nan
fakeY[badcut] = np.nan
fakeY2[badcut] = np.nan

ins1, = plt.plot(fakeT, fakeY, '-r', lw=3, alpha=1.0 )
ins1, = plt.plot(fakeT, fakeY, '.r', ms=0.3, lw=3, alpha=1.0 )
plt.fill_between(fakeT, fakeY, fakeY2, alpha=0.5, color='y')



#ins2, = plt.plot(Ts[Qualify_Flow9], OutflowMetalFrac9[Qualify_Flow9], ':r', ms=5, fillstyle='full', alpha=0.5 )
fakeT = np.copy(Ts[OutFrac9cut])
fakeY = np.copy(OutflowMetalFrac9)
if (usesmoothe):
	OutflowMetalFrac9smooth = smoothe_over400(Ts[OutFrac9cut], OutflowMetalFrac9[OutFrac9cut])
	fakeY = np.copy(OutflowMetalFrac9smooth)
badcut = np.invert(Qualify_Flow9[OutFrac9cut])
fakeT[badcut] = np.nan
fakeY[badcut] = np.nan
ins2, = plt.plot(fakeT, fakeY, ':r', lw=3, alpha=1.0 )
ins2, = plt.plot(fakeT, fakeY, '.r', ms=0.3, lw=3, alpha=1.0 )

#outs1, = plt.plot(Ts[Qualify_InFlow2], InflowMetalFrac2[Qualify_InFlow2], '^b', ms=5, fillstyle='full', alpha=0.5 )
#outs1, = plt.plot(Ts[Qualify_InFlow2], InflowMetalFrac2[Qualify_InFlow2], '-b', ms=5, fillstyle='full', alpha=0.5 )



#EXPERIMENT 1
#fakeT = np.copy(Ts)
#fakeY = np.copy(InflowMetalFrac2)
#fakeY2 = np.copy(InflowMetalFrac2)
#if (usesmoothe):
#	InflowMetalFrac2smooth = smoothe_over400(Ts, InflowMetalFrac2)
#	fakeY = np.copy(InflowMetalFrac2smooth)
#badcut = np.invert(Qualify_InFlow2)
#fakeT[badcut] = np.nan
#fakeY[badcut] = np.nan
#fakeY2[badcut] = np.nan
#outs1, = plt.plot(fakeT[InFrac2cut], fakeY[InFrac2cut], '-b', lw=3, alpha=1.0 )
#outs1, = plt.plot(fakeT[InFrac2cut], fakeY[InFrac2cut], '.b', ms=0.3, lw=3, alpha=1.0 )

fakeT = np.copy(Ts[InFrac2cut])
fakeY = np.copy(InflowMetalFrac2)
fakeY2 = np.copy(InflowMetalFrac2)
if (usesmoothe):
	InflowMetalFrac2smooth = smoothe_over400(Ts[InFrac2cut], InflowMetalFrac2[InFrac2cut])
	fakeY = np.copy(InflowMetalFrac2smooth)
badcut = np.invert(Qualify_InFlow2[InFrac2cut])
fakeT[badcut] = np.nan
fakeY[badcut] = np.nan
fakeY2[badcut] = np.nan
outs1, = plt.plot(fakeT, fakeY, '-b', lw=3, alpha=1.0 )
outs1, = plt.plot(fakeT, fakeY, '.b', ms=0.3, lw=3, alpha=1.0 )




#plt.fill_between(fakeT[InFrac2cut], fakeY2[InFrac2cut], fakeY[InFrac2cut], alpha=0.2, color='b')

#outs2, = plt.plot(Ts[Qualify_InFlow9], InflowMetalFrac9[Qualify_InFlow9], 'ob', ms=5, fillstyle='full', alpha=0.5 )

fakeT = np.copy(Ts[InFrac9cut])
fakeY = np.copy(InflowMetalFrac9)
if (usesmoothe):
	InflowMetalFrac9smooth = smoothe_over400(Ts[InFrac9cut], InflowMetalFrac9[InFrac9cut])
	fakeY = np.copy(InflowMetalFrac9smooth)
	print 'gurara ',InflowMetalFrac9smooth
	print 'bavara ',InflowMetalFrac9
badcut = np.invert(Qualify_InFlow9[InFrac9cut])
fakeT[badcut] = np.nan
fakeY[badcut] = np.nan
outs2, = plt.plot(fakeT, fakeY, ':b', ms=3, alpha=1.0, lw=3)
outs2, = plt.plot(fakeT, fakeY, '.b', ms=1.0, alpha=1.0)

#outs2, = plt.plot(Ts[Qualify_InFlow9], InflowMetalFrac9[Qualify_InFlow9], ':b', ms=5, fillstyle='full', alpha=0.5 )


#ism2, = plt.plot(Ts[metallicity_cut], ISM_metal_fraction[metallicity_cut], ':k', ms=3, alpha=0.25)
fakeT = np.copy(Ts)
fakeY = smoothe_over400(Ts,ISM_metal_fraction)
print len(metallicity_cut), len(realcut)
badcut = np.invert(metallicity_cut*realcut)
print badcut
fakeT[badcut] = np.nan
fakeY[badcut] = np.nan
#ism2, = plt.plot(Ts[metallicity_cut], ISM_metal_fraction[metallicity_cut], '.k', ms=3, alpha=0.25)
ism2, = plt.plot(fakeT, fakeY, '-k', lw=3, alpha=0.5)



plt.yscale('log')
plt.ylim(5e-4, 1.5e0)

plt.xlabel('cosmic time (Myr)',fontsize=20)
plt.ylabel(r'Z/${\rm Z}_\odot$', fontsize=20)
Shells = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95]

innerstring = '0.'+ShellUseInner+r'5 $R_{\rm vir}$'
outerstring = '0.'+ShellUseOuter+r'5 $R_{\rm vir}$ '
if (specialshell):
	outerstring = speciallabel
if (doublespecial):
	innerstring = doublespeciallabel


#thar22b =['0.'+ShellUseInner+r'5 $R_{\rm vir}$ outflow', outerstring+r' $R_{\rm vir}$ outflow', '0.'+ShellUseInner+r'5 $R_{\rm vir}$ inflow', outerstring+r' $R_{\rm vir}$ inflow', 'ISM metallicity', 'halo metallicity']
thar22b =[innerstring+r' outflow', outerstring+r' outflow', innerstring+r' inflow', outerstring+r' inflow', 'ISM metallicity', 'CGM metallicity']


#plt.legend(thar22b, loc = 'best',  prop={'size':12})
plt.legend( [(in1, ins1), (in2, ins2), (out1, outs1), (out2, outs2), ism2, ism3],thar22b, loc = 'upper left',  prop={'size':legsize}, ncol=2)
ax = plt.gca()

x_range = ax.get_xlim()
y_range = ax.get_ylim()

ax2 = ax.twiny()


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

ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(T_wanted)
ax2.set_xticklabels(z_wanted_str)
ax2.set_xlabel('z', fontsize=20)
ticklabels = ax2.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(18)
    label.set_family('serif')
ticklabels = ax2.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(18)
    label.set_family('serif')
ax2.xaxis.set_tick_params(width=2.5, length=8)
ax2.yaxis.set_tick_params(width=2.5, length=8)    

plt.savefig(name+'/'+ str(today) + 'inoutmetal'+haloN+'.pdf')
plt.clf()


fig2 = plt.figure(figsize=(10,8))
plt.plot(Ts, Flow_in, color='b', lw=2)
plt.plot(Ts, MetalFlow_in*1000, color='k',  ls='--', lw=2)

plt.plot(Ts, Flow_out, color='r', lw=2)
plt.plot(Ts, MetalFlow_out*1000, color='g',  ls='--', lw=2)

plt.plot(Ts, SFR1*10, color='k',  lw=2)


plt.plot(Ts, Flow_in2, color='b',  lw=2)
plt.plot(Ts, MetalFlow_in2*1000, color='k', ls='--', lw=2)

plt.plot(Ts, Flow_out2, color='r',  lw=2)
plt.plot(Ts, MetalFlow_out2*1000, color='g', ls='--', lw=2)


plt.xlabel('cosmic time (Myr)',fontsize=20)
plt.ylabel(r'Outflow, Inflow  ($M_\odot$${\rm yr}^{-1}$)',fontsize=20)
#plt.yscale('log')
#ShellUseInner = '2'
#ShellUseOuter = '9'

thar22b =['0.'+ShellUseInner+r'5 $R_{\rm vir}$ out/inflow', '0.'+ShellUseInner+r'5 $R_{\rm vir}$ metal out/inflow $\times$ 1000', outerstring+r' $R_{\rm vir}$ out/inflow', outerstring+r' $R_{\rm vir}$ metal out/inflow $\times$ 1000']
thar22b =[innerstring+' out/inflow', innerstring+r' metal out/inflow $\times$ 1000', outerstring+r' out/inflow', outerstring+r' metal out/inflow $\times$ 1000', r'SFR $\times$ 10']

plt.legend(thar22b, loc = 'best',  prop={'size':legsize})

ax = plt.gca()

x_range = ax.get_xlim()
y_range = ax.get_ylim()

ax2 = ax.twiny()


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

ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(T_wanted)
ax2.set_xticklabels(z_wanted_str)
ax2.set_xlabel('z', fontsize=20)
ticklabels = ax2.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(18)
    label.set_family('serif')
ticklabels = ax2.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(18)
    label.set_family('serif')
ax2.xaxis.set_tick_params(width=2.5, length=8)
ax2.yaxis.set_tick_params(width=2.5, length=8)    

plt.savefig(name+'/'+ str(today) + 'variousflows'+haloN+'.pdf')
plt.clf()



fig2b = plt.figure(figsize=(11,8.5))

Flow_insmooth = smoothe_over400(Ts, Flow_in)
Flow_outsmooth = smoothe_over400(Ts, Flow_out)
Flow_in2smooth = smoothe_over400(Ts, Flow_in2)
Flow_out2smooth = smoothe_over400(Ts, Flow_out2)
MetalFlow_insmooth = smoothe_over400(Ts, MetalFlow_in)
MetalFlow_outsmooth = smoothe_over400(Ts, MetalFlow_out)
MetalFlow_in2smooth = smoothe_over400(Ts, MetalFlow_in2)
MetalFlow_out2smooth = smoothe_over400(Ts, MetalFlow_out2)
SFRsmooth = smoothe_over400(Ts, SFR1)

#plt.plot(Ts, Flow_in, color='b', lw=2)
#plt.plot(Ts, MetalFlow_in*1000, color='k',  ls='--', lw=2)

#plt.plot(Ts, Flow_out, color='r', lw=2)
#plt.plot(Ts, MetalFlow_out*1000, color='g',  ls='--', lw=2)

#plt.plot(Ts, Flow_in2, color='b',  lw=2)
#plt.plot(Ts, MetalFlow_in2*1000, color='k', ls=':', lw=2)

#plt.plot(Ts, Flow_out2, color='r',  lw=2)
#plt.plot(Ts, MetalFlow_out2*1000, color='g', ls=':', lw=2)


plt.plot(Ts, Flow_insmooth, color='b', lw=2, alpha=0.5)
plt.plot(Ts, MetalFlow_insmooth*Metalplotboost, color='b',  ls='--', lw=2)

plt.plot(Ts, Flow_outsmooth, color='r', lw=2, alpha=0.5)
plt.plot(Ts, MetalFlow_outsmooth*Metalplotboost, color='r',  ls='--', lw=2)

plt.plot(Ts, SFRsmooth*SFRplotboost, color='k',  lw=2, alpha=0.5)
#plt.plot(Ts, SFR1*10, color='y', ls =':', lw=2)


plt.plot(Ts, Flow_in2smooth, color='b',  lw=2, alpha=0.5)
plt.plot(Ts, MetalFlow_in2smooth*Metalplotboost, color='b', ls='--', lw=2)

plt.plot(Ts, Flow_out2smooth, color='r',  lw=2, alpha=0.5)
plt.plot(Ts, MetalFlow_out2smooth*Metalplotboost, color='r', ls='--', lw=2)

plt.xlabel('cosmic time (Myr)',fontsize=20)
plt.ylabel(r'Outflow, Inflow  ($M_\odot$${\rm yr}^{-1}$)',fontsize=20)
#plt.yscale('log')
#ShellUseInner = '2'
#ShellUseOuter = '9'

thar22b =['0.'+ShellUseInner+r'5 $R_{\rm vir}$ out/inflow', '0.'+ShellUseInner+r'5 $R_{\rm vir}$ metal out/inflow $\times$ 1000', outerstring+r' $R_{\rm vir}$ out/inflow', outerstring+r' $R_{\rm vir}$ metal out/inflow $\times$ 1000']
thar22b =[innerstring+' out/inflow', innerstring+r' metal out/inflow$\times$'+str(Metalplotboost), outerstring+r' out/inflow', outerstring+r' metal out/inflow$\times$'+str(Metalplotboost), r'SFR$\times$'+str(SFRplotboost)]

plt.legend(thar22b, loc = 'best',  prop={'size':legsize})

ax = plt.gca()

x_range = ax.get_xlim()
y_range = ax.get_ylim()

ax2 = ax.twiny()


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

ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(T_wanted)
ax2.set_xticklabels(z_wanted_str)
ax2.set_xlabel('z', fontsize=20)
ticklabels = ax2.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(18)
    label.set_family('serif')
ticklabels = ax2.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(18)
    label.set_family('serif')
ax2.xaxis.set_tick_params(width=2.5, length=8)
ax2.yaxis.set_tick_params(width=2.5, length=8)    

plt.savefig(name+'/'+ str(today) + 'variousflowsB'+haloN+'.pdf')
plt.clf()


import matplotlib.ticker as ticker
def y_fmt(x,y):
    return "{0:.6g}".format(x)

fig2bb = plt.figure(figsize=(11.5,8.5))



plt.plot(Ts, SFRsmooth*SFRplotboost, color='k',  lw=2, alpha=1)


plt.plot(Ts, MetalFlow_insmooth*Metalplotboost, color='r',  ls='--', lw=3)

plt.plot(Ts, MetalFlow_outsmooth*Metalplotboost, color='r',  ls=':', lw=3)



plt.plot(Ts, MetalFlow_in2smooth*Metalplotboost, color='b', ls='--', lw=3)

plt.plot(Ts, MetalFlow_out2smooth*Metalplotboost, color='b', ls=':', lw=3)

plt.xlabel('cosmic time (Myr)',fontsize=20)
plt.ylabel(r'Outflow, Inflow  ($M_\odot$${\rm yr}^{-1}$)',fontsize=20)
#plt.yscale('log')
#ShellUseInner = '2'
#ShellUseOuter = '9'

thar22b =[r'SFR$\times$'+str(SFRplotboost), innerstring+r' metal outflow $\times$'+str(Metalplotboost), outerstring+r' metal outflow$\times$'+str(Metalplotboost), innerstring+r' metal inflow$\times$'+str(Metalplotboost), outerstring+r' metal inflow $\times$'+str(Metalplotboost), outerstring+r' outflow', outerstring+r' metal inflow$\times$'+str(Metalplotboost) ]

plt.legend(thar22b, loc = 'best', ncol=1,  prop={'size':legsize})

ax = plt.gca()

x_range = ax.get_xlim()
y_range = ax.get_ylim()

ax2 = ax.twiny()


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
ax.yaxis.set_major_formatter(ticker.FuncFormatter(y_fmt))
ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(T_wanted)
ax2.set_xticklabels(z_wanted_str)
ax2.set_xlabel('z', fontsize=20)
ticklabels = ax2.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(18)
    label.set_family('serif')
ticklabels = ax2.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(18)
    label.set_family('serif')
ax2.xaxis.set_tick_params(width=2.5, length=8)
ax2.yaxis.set_tick_params(width=2.5, length=8)    

plt.savefig(name+'/'+ str(today) + 'variousflowsD'+haloN+'.pdf')
plt.clf()







fig2c = plt.figure(figsize=(10,8))

Flow_insmooth = smoothe_over400(Ts, Flow_in)
Flow_outsmooth = smoothe_over400(Ts, Flow_out)
Flow_in2smooth = smoothe_over400(Ts, Flow_in2)
Flow_out2smooth = smoothe_over400(Ts, Flow_out2)
MetalFlow_insmooth = smoothe_over400(Ts, MetalFlow_in)
MetalFlow_outsmooth = smoothe_over400(Ts, MetalFlow_out)
MetalFlow_in2smooth = smoothe_over400(Ts, MetalFlow_in2)
MetalFlow_out2smooth = smoothe_over400(Ts, MetalFlow_out2)

#plt.plot(Ts, Flow_in, color='b', lw=2)
#plt.plot(Ts, MetalFlow_in*1000, color='k',  ls='--', lw=2)

#plt.plot(Ts, Flow_out, color='r', lw=2)
#plt.plot(Ts, MetalFlow_out*1000, color='g',  ls='--', lw=2)

#plt.plot(Ts, Flow_in2, color='b',  lw=2)
#plt.plot(Ts, MetalFlow_in2*1000, color='k', ls=':', lw=2)

#plt.plot(Ts, Flow_out2, color='r',  lw=2)
#plt.plot(Ts, MetalFlow_out2*1000, color='g', ls=':', lw=2)


plt.plot(Ts[OutFrac2cut], Flow_insmooth[OutFrac2cut], color='r', lw=2)
plt.plot(Ts[OutFrac2cut], MetalFlow_insmooth[OutFrac2cut]*1000, color='r',  ls='-.', lw=2)

plt.plot(Ts[OutFrac9cut], Flow_outsmooth[OutFrac9cut], color='g', lw=2)
plt.plot(Ts[OutFrac9cut], MetalFlow_outsmooth[OutFrac9cut]*1000, color='g',  ls='-.', lw=2)

plt.plot(Ts[InFrac2cut], -1*Flow_in2smooth[InFrac2cut], ls='--', color='b',  lw=2)
plt.plot(Ts[InFrac2cut], -1*MetalFlow_in2smooth[InFrac2cut]*1000, color='b', ls=':', lw=2)

plt.plot(Ts[InFrac9cut], -1*Flow_out2smooth[InFrac9cut], ls='--', color='k',  lw=2)
plt.plot(Ts[InFrac9cut], -1*MetalFlow_out2smooth[InFrac9cut]*1000, color='k', ls=':', lw=2)

plt.yscale('log')

plt.xlabel('cosmic time (Myr)',fontsize=20)
plt.ylabel(r'Outflow, Inflow  ($M_\odot$${\rm yr}^{-1}$)',fontsize=20)
#plt.yscale('log')
#ShellUseInner = '2'
#ShellUseOuter = '9'

thar22b =['0.'+ShellUseInner+r'5 $R_{\rm vir}$ out/inflow', '0.'+ShellUseInner+r'5 $R_{\rm vir}$ metal out/inflow $\times$ 1000', outerstring+r' $R_{\rm vir}$ out/inflow', outerstring+r' $R_{\rm vir}$ metal out/inflow $\times$ 1000']
thar22b =[innerstring+' out/inflow', innerstring+r' metal out/inflow $\times$ 1000', outerstring+r' out/inflow', outerstring+r' metal out/inflow $\times$ 1000',]
thar22b =[innerstring+' outflow', innerstring+r' metal outflow $\times$ 1000', outerstring+r' outflow', outerstring+r' metal outflow $\times$ 1000',  innerstring+' inflow', innerstring+r' metal inflow $\times$ 1000', outerstring+r' inflow', outerstring+r' metal inflow $\times$ 1000',  ]


plt.legend(thar22b, loc = 'best',  prop={'size':(legsize*0.75)})

ax = plt.gca()

x_range = ax.get_xlim()
y_range = ax.get_ylim()

ax2 = ax.twiny()


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

ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(T_wanted)
ax2.set_xticklabels(z_wanted_str)
ax2.set_xlabel('z', fontsize=20)
ticklabels = ax2.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(18)
    label.set_family('serif')
ticklabels = ax2.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(18)
    label.set_family('serif')
ax2.xaxis.set_tick_params(width=2.5, length=8)
ax2.yaxis.set_tick_params(width=2.5, length=8)    

plt.savefig(name+'/'+ str(today) + 'variousflowsC'+haloN+'.pdf')
plt.clf()








###################################################################
#writing text
#ML_integrated
###################################################################

total_expelled = np.trapz(Flow_in, Ts, 0.1)
total_formed =  np.trapz(SFR1, Ts, 0.1)
total_metal_expelled = np.trapz(MetalFlow_in, Ts, 0.1)
total_expelled9 = np.trapz(Flow_out, Ts, 0.1)
total_metal_expelled9 = np.trapz(MetalFlow_out, Ts, 0.1)
total_metal_gained = np.trapz(MetalFlow_in2, Ts, 0.1)
total_metal_gained9 = np.trapz(MetalFlow_out2, Ts, 0.1)
total_mass_gained = np.trapz(Flow_in2, Ts, 0.1)
total_mass_gained9 = np.trapz(Flow_out2, Ts, 0.1)


ISM_metal_fraction = (ISMShellMetal[0] + ISMShellMetal[1]) / (ISMShellMass[0] + ISMShellMass[1])
if (InnerISM_only):
	print 'what the fuckers '
	print ISMShellMetal[0], 'metals'
	print ISMShellMass[0], 'mass'
	ISM_metal_fraction = (ISMShellMetal[0]) / (ISMShellMass[0])

TotalISMMetal = np.sum( (ISMShellMetal[0] + ISMShellMetal[1]))
TotalISMMass = np.sum((ISMShellMass[0] + ISMShellMass[1]))

if (InnerISM_only):
	TotalISMMetal = np.sum( ISMShellMetal[0] )
	TotalISMMass = np.sum(ISMShellMass[0] )


ISM_metal_FS = (ISMShellMetal[0] + ISMShellMetal[1])
ISM_mass_FS = (ISMShellMass[0] + ISMShellMass[1])

if (InnerISM_only):
	ISM_metal_FS = ISMShellMetal[0]
	ISM_mass_FS = ISMShellMass[0]
	ISM_metal_FS = ISMMass
	ISM_mass_FS = ISM_TMass

AverageHaloZ = np.sum(ISMShellMetal[0:10]) / np.sum(ISMShellMass[0:10])
print Halo_metal_fraction

AverageISMZ = TotalISMMetal/TotalISMMass
#realcut = ISM_metal_fraction > 1e-10
#fakecut = np.invert(realcut)

#for theISM in AverageISMZ
	
realcut = np.isreal(ISM_metal_fraction)
poscut = ISM_metal_fraction > 1e-10 
realcut *= poscut

TimeAverage_ISMZ = np.trapz(ISM_metal_fraction[realcut], Ts[realcut], 0.1) / (Ts[-1] - Ts[0])

ML_from_total = total_expelled/total_formed
ML_from_total9 = total_expelled9/ total_formed
MetalLoad_from_total = total_metal_expelled / total_formed
MetalLoad_from_total9 = total_metal_expelled9 / total_formed

MG_from_total = total_mass_gained/total_formed
MG_from_total9 = total_mass_gained9/total_formed
MetalGain_from_total = total_metal_gained/total_formed
MetalGain_from_total9 = total_metal_gained9 / total_formed


print 'mass loading '
print 'Shell 2 flux cumulative ', ML_from_total
print 'Shell 2 metal flux cumulative *1000', MetalLoad_from_total*1000
print 'average ISM metal fraction ', AverageISMZ

MLline_list = [int(haloN), ML_from_total, MetalLoad_from_total, AverageISMZ, TimeAverage_ISMZ, AverageHaloZ, ML_from_total9, MetalLoad_from_total9, MG_from_total, MG_from_total9, MetalGain_from_total, MetalGain_from_total9,   zstart, zend]

#totalmasses, crossmasses, intotalmasses, incrosstotalmasses

print MLline_list

MLline = line_to_string(MLline_list)

g = open(name+'/'+'Z_integrated.txt', 'a')
g.write(MLline)
g.close()

fig2 = plt.figure(figsize=(10,8))
#plt.plot(Ts, ISM_metal_FS*1e10/little_h, '.k')
#plt.plot(Ts, ISM_mass_FS*1e10/little_h, '.b')
#plt.plot(Ts, Halo_metal_FS*1e10/little_h, '.r')
#plt.plot(Ts, Halo_mass_FS*1e10/little_h, '.g')

plt.plot(Ts, smoothe_over400(Ts,ISM_metal_FS)/ISM_metal_FS[-1], '-k', lw=3)
plt.plot(Ts, smoothe_over400(Ts,ISM_mass_FS)/ISM_mass_FS[-1], '--k', lw=3)
plt.plot(Ts, smoothe_over400(Ts,Halo_metal_FS)/Halo_metal_FS[-1], '-r', lw=3)
plt.plot(Ts, smoothe_over400(Ts,Halo_mass_FS)/Halo_mass_FS[-1], '--r', lw=3)




plt.legend(['ISM metals ('+"{0:.2g}".format(ISM_metal_FS[-1])+' at z=0)', 'ISM gas ('+"{0:.2g}".format(ISM_mass_FS[-1])+' at z=0)', 'CGM metals ('+"{0:.2g}".format(Halo_metal_FS[-1])+' at z=0)', 'CGM gas ('+"{0:.2g}".format(Halo_mass_FS[-1])+' at z=0)'], loc = 'best',  prop={'size':legsize})


plt.ylabel('fraction of z=0 value', fontsize=20)
plt.xlabel('cosmic time (Myr)', fontsize=20)

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
#plt.yscale('log')
plt.yscale('linear')

plt.savefig(name+'/'+ str(today) + 'buildup'+haloN+'.pdf')
plt.clf()

fig2b = plt.figure(figsize=(10,8))
#plt.plot(Ts, ISM_metal_FS*1e10/little_h, '.k')
#plt.plot(Ts, ISM_mass_FS*1e10/little_h, '.b')
#plt.plot(Ts, Halo_metal_FS*1e10/little_h, '.r')
#plt.plot(Ts, Halo_mass_FS*1e10/little_h, '.g')

plt.plot(Ts, smoothe_over400(Ts,ISM_metal_FS)/ISM_metal_FS[-1], '-k', lw=3)
plt.plot(Ts, smoothe_over400(Ts,ISM_mass_FS)/ISM_mass_FS[-1], '--k', lw=3)
plt.plot(Ts, smoothe_over400(Ts,Halo_metal_FS)/Halo_metal_FS[-1], '-r', lw=3)
plt.plot(Ts, smoothe_over400(Ts,Halo_mass_FS)/Halo_mass_FS[-1], '--r', lw=3)




plt.legend(['ISM metals ('+"{0:.2g}".format(ISM_metal_FS[-1])+' at z=0)', 'ISM gas ('+"{0:.2g}".format(ISM_mass_FS[-1])+' at z=0)', 'CGM metals ('+"{0:.2g}".format(Halo_metal_FS[-1])+' at z=0)', 'CGM gas ('+"{0:.2g}".format(Halo_mass_FS[-1])+' at z=0)'], loc = 'best',  prop={'size':legsize})


plt.ylabel('fraction of z=0 value', fontsize=20)
plt.xlabel('cosmic time (Myr)', fontsize=20)
plt.yscale('log')
plt.ylim([1e-2, 1e1])

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
#plt.yscale('log')
#plt.yscale('linear')

plt.savefig(name+'/'+ str(today) + 'buildupB'+haloN+'.pdf')
plt.clf()




fig3a = plt.figure(figsize=(10,8))

dT = Ts[1:] - Ts[:-1]
dT = np.append([dT[0]], dT)
IntegratedMetalFlux = np.cumsum(dT*MetalFlow_out)*1e6
IntegratedMetalGained = -1*np.cumsum(dT*MetalFlow_out2)*1e6
plt.plot(Ts/1000, IntegratedMetalFlux, '-r', lw=2)
plt.plot(Ts/1000, IntegratedMetalGained, '--b', lw=2)
plt.plot(Ts/1000, ISMMass, ':g',lw=2)
plt.plot(starTs/1000, starmetalcum, '-k',lw=2)
plt.ylabel('Mass of metals (Msun)',fontsize=20)
plt.xlabel('cosmic time (Gyr)', fontsize=20)
plt.legend(['Metals Lost', 'Metals Gained', 'ISM metals', 'stellar metals'], loc = 'best',  prop={'size':legsize})

ax = plt.gca()

x_range = ax.get_xlim()
y_range = ax.get_ylim()

ax2 = ax.twiny()


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

ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(T_wanted)
ax2.set_xticklabels(z_wanted_str)
ax2.set_xlabel('z', fontsize=20)
ticklabels = ax2.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(18)
    label.set_family('serif')
ticklabels = ax2.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(18)
    label.set_family('serif')
ax2.xaxis.set_tick_params(width=2.5, length=8)
ax2.yaxis.set_tick_params(width=2.5, length=8)   

plt.savefig(name+'/'+ str(today) + 'flower'+haloN+'.pdf')
plt.clf()

fig3b = plt.figure(figsize=(10,8))

#dT = Ts[1:] - Ts[:-1]
#dT = np.append([dT[0]], dT)
#IntegratedMetalFlux = np.cumsum(dT*MetalFlow_out)*1e6
#IntegratedMetalGained = -1*np.cumsum(dT*MetalFlow_out2)*1e6
plt.plot(Ts/1000, IntegratedMetalFlux, '-r', lw=2)
plt.plot(Ts/1000, IntegratedMetalGained, '--b', lw=2)
plt.plot(Ts/1000, ISMMass, ':g',lw=2)
plt.plot(starTs/1000, starmetalcum, '-k',lw=2)
plt.yscale('log')
plt.ylabel('Mass of metals (Msun)',fontsize=20)
plt.xlabel('cosmic time (Gyr)', fontsize=20)
plt.legend(['Metals Lost', 'Metals Gained', 'ISM metals', 'stellar metals'], loc = 'best',  prop={'size':legsize})

ax = plt.gca()

x_range = ax.get_xlim()
y_range = ax.get_ylim()
plt.ylim(y_range[1]/1e4, y_range[1])
ax2 = ax.twiny()


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

ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(T_wanted)
ax2.set_xticklabels(z_wanted_str)
ax2.set_xlabel('z', fontsize=20)
ticklabels = ax2.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(18)
    label.set_family('serif')
ticklabels = ax2.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(18)
    label.set_family('serif')
ax2.xaxis.set_tick_params(width=2.5, length=8)
ax2.yaxis.set_tick_params(width=2.5, length=8)   

plt.savefig(name+'/'+ str(today) + 'flowerlog'+haloN+'.pdf')
plt.clf()



print 'total metal expelled 9',total_metal_expelled9
print 'total expelled 9', total_expelled9

print 'now what the fuck ',ISM_metal_fraction

print' fuckuing ',  ISMShellMetal[0] 
print 'fuck ', ISMShellMass[0]

print 'kar ',len(SFRsmooth), len(SFR1), len(Ts)
