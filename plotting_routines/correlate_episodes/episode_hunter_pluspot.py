import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from datetime import date
import sys
import scipy.optimize as optimize

today = date.today()
print today

plot_markers = 'n'
naked = True
interactive_mode = True
jacket = False
specialshell = True #use specialval Extra Shell instead of whatever input was
doublespecial = True
specialval = 1 #for 1.0 Rvir, specialval = 5
Dspecialval = 1  #0 for 10 kpc, 1 for 25 kpc, 2 for 50 kpc - confusing but this is inner region
paramfac_flies = True
SFR_correct_factor = 1.15  #to correct for stars formed between snapshots 

doublespeciallabel = '10 kpc'
speciallabel = '25 kpc'

surf_den_boost=0
SFR_boost = 1

#these parameters should be set by command prompt
zstart =999
zend = 999
stagger = 999 #the stagger between start of integration for SFR and for outflow rate.

little_h = 0.7 
paramfac = 1e5
adaptive_paramfac = 'y'

#z_wanted = 


if (len(sys.argv) < 2):
	print 'syntax  blah.py Nhalo filename(outflow) filename2(inflow)  finname3(rhalf) finname4 (rhalf extra) [zstart zend stagger N inner Shell, N outer shell]'
	sys.exit()

haloN = str(sys.argv[1])
filename = str(sys.argv[2])
filename2 = str(sys.argv[3])
rhalf_finname = str(sys.argv[4])
pot_finname = str(sys.argv[5])


#sigh i suck at coding 
uselims = 'n'  #this is only for the y limits on the plots. Defaults set below 
subtract_delay = 'y'  #this is to "go back in time" after the SFR is computed. 
slack_time = 'n' #this one applies to whether the outflow episode can be prolonged by the "Outflow_1st_par" things. 
slack_delay = 'y'   #this is to prolong an SF episode after the fixed "delay" given. THe only meaning of "delay" is how long the SF episode is integrate for
stagger_outflow = 'y' # this is to officially start the integration for outflow total at a time after the SF integration. Given by "stagger"

the_y_lims_zoomed = (-15, 15) 
#default for 0.25 and 0.95 Rvir
ShellUseInner = '2'
ShellUseInnerOrig = ShellUseInner
ShellUseOuter = '9'

if (len(sys.argv) >= 7):
	zstart = float(sys.argv[6])
	zend = float(sys.argv[7])
	stagger = float(sys.argv[8])
	print 'got zstart, zend, stagger from command prompt'


if (len(sys.argv) >= 10):
	ShellUseInner = str(sys.argv[8])
	ShellUseOuter = str(sys.argv[9])
	print 'got Shell Inner and Outer from command prompt'

zmedpoint = 3

if (zstart > 3.5):
	z_wanted = np.array([4.5, 4.0, 3.5, 3.0, 2.5, 2.0])
	zmedpoint = 3.0
if (zstart > 1 and zstart < 3.5):
	z_wanted = np.array([2.0, 1.5, 1.25, 1.0, 0.75, 0.5])
	zmedpoint = 1.25
if (zstart < 1):
	z_wanted = np.array([0.5, 0.25, 0.1, 0.0])
	zmedpoint = 0.25	

#z_wanted = np.array([4.5, 4.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.25, 1.0, 0.75, 0.5, 0.25, 0.1, 0.0])

z_wanted_cut = (z_wanted <= zstart)
z_wanted_cut *= (z_wanted >= zend)
z_wanted = z_wanted[z_wanted_cut]
print 'using the following z_wanted array ', z_wanted

zmedpoint = max(zmedpoint, zend)

z_wanted_str = []
for q in z_wanted:
	z_wanted_str.append(str(q))
a_wanted = 1.0 / (1 + z_wanted)

print 'using shells ', ShellUseInner, ShellUseOuter
SFR_smoothing = 20 #Myr
delay = stagger # minimal length of SF interval
interval = stagger #minimal length of outflow interval lets keep interval = delay for now
max_interval = 250 + (stagger-60)

#delay = stagger
#interval = stagger
# i think this thing serves no function


#for initiating


def line_to_string(line, newline = 'y'):
	linestring = ''
	for q in line:
		linestring += "{0:.6g}".format(q) + '  '
	if (newline == 'y'):
		linestring+=' \n'
	return linestring

finname = 'tasz.txt'
f = open(finname)
dars = np.loadtxt(f)
theas = dars[:,0]
thets = dars[:,1]
f.close()

finname = filename
f = open(finname)
dars = np.loadtxt(f)
Zs = dars[:,1]
Ns = dars[:,0]
cut1 = Zs >= zend
cut2 = Zs <= zstart
cuts = cut1*cut2


Zs = Zs[cuts] 
Ns = Ns[cuts]
zend = max(zend, np.min(Zs))

Rvirused = dars[:,2][cuts]
SFRInner = dars[:,3][cuts]
SFRouter = dars[:,4][cuts]
#SFR1 = SFRInner + SFRouter
#for now excluding SF outside 0.2 Rvir
SFR1 = SFRInner * SFR_correct_factor
DenseGasMass = dars[:,5][cuts]
TotalMass = dars[:,6][cuts]
Mgas = dars[:,7][cuts]
Mstar = dars[:,8][cuts]

fbar = (Mgas + Mstar)/TotalMass
fgas = Mgas/TotalMass
fstar = Mstar/TotalMass




totalSmass = dars[:,9][cuts]
Flow1 = dars[:,10][cuts]
NewSFR_timescale = dars[:,11][cuts]


newtonG = 6.67384e-8
sSFR = SFR1/ (Mstar/little_h) 

Vcirc = np.sqrt((newtonG * TotalMass *  2e33 / little_h) / (Rvirused/little_h * 3.08e21))
Vcirc = Vcirc / 1e5
#Vcirc *= np.sqrt(1.0 + Zs)


getvcirc = interp1d(Zs[::-1], Vcirc[::-1], kind='linear')
getmass = interp1d(Zs[::-1], TotalMass[::-1], kind='linear')
getRvir  = interp1d(Zs[::-1], Rvirused[::-1], kind='linear')
getMstar = interp1d(Zs[::-1], Mstar[::-1], kind='linear')


print 'vcirc at z=',zend,getvcirc(zend),' km/s'
print 'mass at z=',zend,getmass(zend),' Msun/h'
print 'rvir at z=',zend,getRvir(zend),' kpc/h'
print 'Mstar at z=2',zend,getMstar(zend), ' Msun/h'

mend = float(getmass(zend))
vend = float(getvcirc(zend))
mstarend = float(getMstar(zend))

#mstart = float(getmass(zstart))
#mstarstart = float(getMstar(zstart))
#vstart = float(getvcirc(zstart))

mstart = TotalMass[0]
mstarstart = Mstar[0]
vstart = Vcirc[0]

if (adaptive_paramfac == 'y'):
	#paramfac =  20 *  (1e11/mend)
	#paramfac = 100 * 1e10/mstarend
	paramfac = max(1 * (1e10/mstarend)**2.0, 1.0)

SFR_1stpar = 1.0/paramfac
SFR_2ndpar = 1.0
SFR_3rdpar = 5.0/paramfac

#for slack
SFR_4thpar = 1.0/paramfac
SFR_5thpar = 1.0
SFR_6rhpar = 2.0/paramfac

def paramfac_fly(theMstar):
	myparamfac = max(1 * (1e10/theMstar)**2.0, 1.0)
	return myparamfac


massstring = "{0:.4g}".format(float(getmass(zend))/little_h)
#massstring = '%.3E' % Decimal(massstring)
#Flow1 = dars[:,7][cuts]
#Flow_Shells = []

#Cross_Shells = [] 


#needs to be one less than the actual number used - 0 counts.
NumShells = 20
Flow_Shells = []
Cross_Shells = []



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


starter = -1 
count = 0
while (count < NumShells):
	anticount = NumShells - count - 1
	riteindex = -6*anticount+starter
	TheCross = (dars[:,riteindex][cuts]* 1e10 / little_h) / NewSFR_timescale
	#print 'The Cross of Shell ', count,' is index ',riteindex
	Cross_Shells.append(TheCross)
	count+=1

count = 0
Extra_CrossShells = []
while (count < (NumShells + NumExtraShells)):
	anticount = NumShells + NumExtraShells - count - 1
	riteindex = -6*anticount+starter
	TheCross= (dars[:,riteindex][cuts]* 1e10 / little_h) / NewSFR_timescale
	Extra_CrossShells.append(TheCross)
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

NPartCross_ar=np.array(NPartCross_ar)
NPartFlux_ar=np.array(NPartFlux_ar)
NpartShell_ar=np.array(NpartShell_ar)


partmass = 2.7e4

Outflow_fraction = NPartFlux_ar/NpartShell_ar
Total_Outflow_Flux = (NPartFlux_ar* partmass / little_h) / NewSFR_timescale

TotalOutflow_fraction = sum(NPartFlux_ar[0:10])/sum(NpartShell_ar[0:10])

TheTotalOutflow = sum(Total_Outflow_Flux[1:10])
#TheTotalOutflow = np.copy(Total_Outflow_Flux[1])
#count = 2 
#while (count <= 9):
#	TheTotalOutflow += Total_Outflow_Flux[count]
#	count+=1 


f.close()

	

Flow_in = Flow_Shells[int(ShellUseInner)]
Flow_out = Flow_Shells[int(ShellUseOuter)]

Cross_in = Cross_Shells[int(ShellUseInner)]
Cross_out = Cross_Shells[int(ShellUseOuter)]

if (specialshell):
	Flow_out = Extra_FlowShells[int(specialval)]
	Cross_out = Extra_CrossShells[int(specialval)]
if (doublespecial):
	Flow_in = Extra_FlowShells[int(Dspecialval)]
	Cross_in = Extra_CrossShells[int(Dspecialval)]	

maxy = min((max(np.max(SFR1), np.max(Flow_in))), 100)

finname2 = filename2
f2 = open(finname2)
dars2 = np.loadtxt(f2)

InFlow_Shells = []
InCross_Shells = []


starter = -5 
count = 0
while (count < NumShells):
	anticount = NumShells - count  -1 
	riteindex = -6*anticount+starter
	TheFlow = dars2[:,riteindex][cuts]
	#print 'The InFlow of Shell ', count,' is index ',riteindex
	InFlow_Shells.append(TheFlow)
	count +=1

count = 0
Extra_InFlowShells = []
while (count < (NumShells + NumExtraShells)):
	anticount = NumShells + NumExtraShells - count - 1
	riteindex = -6*anticount+starter
	TheFlow= dars2[:,riteindex][cuts]
	Extra_InFlowShells.append(TheFlow)
	count += 1



starter = -1 
count = 0
while (count < NumShells):
	anticount = NumShells - count  -1
	riteindex = -6*anticount+starter
	TheCross = (dars2[:,riteindex][cuts]* 1e10 / little_h) / NewSFR_timescale * -1 
	#print 'The Cross of Shell ', count,' is index ',riteindex
	InCross_Shells.append(TheCross)
	count+=1
	
count = 0
Extra_InCrossShells = []
while (count < (NumShells + NumExtraShells)):
	anticount = NumShells + NumExtraShells - count - 1
	riteindex = -6*anticount+starter
	TheCross= (dars2[:,riteindex][cuts]* 1e10 / little_h) / NewSFR_timescale
	Extra_InCrossShells.append(TheCross)
	count += 1


Flow_in2 = InFlow_Shells[int(ShellUseInner)]
Flow_out2 = InFlow_Shells[int(ShellUseOuter)]

Cross_in2 = InCross_Shells[int(ShellUseInner)] 
Cross_out2 = InCross_Shells[int(ShellUseOuter)]


if (specialshell):
	Flow_out2 = Extra_InFlowShells[int(specialval)]
	Cross_out2 = Extra_InCrossShells[int(specialval)]
	ShellUseOuter = 'S'+str(specialval)
if (doublespecial):
	Flow_in2 = Extra_InFlowShells[int(Dspecialval)]
	Cross_in2 = Extra_InCrossShells[int(Dspecialval)]
	ShellUseInner = 'S'+str(Dspecialval)

count = 0
#print 'showing infall rates'
while (count < len(Flow_in2)):
	#print count, Zs[count], Cross_in2[count], Flow_in2[count] 
	count +=1

#Flow_in2 = dars2[:,-31][cuts]
#Flow_out2 = dars2[:,-3][cuts]
f2.close()


miny = max(-50, np.min(Flow_in2))
the_y_lims = (miny, maxy)


As = 1.0 / (1.0 + Zs)

tfora = interp1d(theas, thets)
afort = interp1d(thets, theas)
Ts = tfora(As)


Ts *= 1000.0

T_wanted = tfora(a_wanted)
T_wanted *= 1000


total_expelled = np.trapz(Flow_in, Ts, 0.1)
total_formed =  np.trapz(SFR1, Ts, 0.1)
total_expelled_cross = np.trapz(Cross_in, Ts, 0.1)
total_gained = np.trapz(Flow_in2, Ts, 0.1)
total_gained_cross = np.trapz(Cross_in2, Ts, 0.1)

totalmasses = []
cumtotalmasses = []

half_total_expelled = total_expelled/2.0

for shell in Flow_Shells:
	cumtotalmass = []
	if (len(shell) > 1): 
		thetotalmass = np.trapz(shell, Ts, 0.1)
		count = 0 
		while (count < len(shell)):
			mass_so_far =np.trapz(shell[0:count+1], Ts[0:count+1], 0.1)
			cumtotalmass.append(mass_so_far)			
			count+=1
	else: thetotalmass = 0
	cumtotalmass = np.array(cumtotalmass)
	totalmasses.append(thetotalmass)
	cumtotalmasses.append(cumtotalmass)

count = 0 
magiccount = 0 
for p in cumtotalmasses[int(ShellUseInnerOrig)]:
	#print 'doo be doo ',count, Ns[count], p, half_total_expelled
	if (p >= half_total_expelled):
		print 'half outflow point reached!' 
		print count, Ns[count], Zs[count]
		magiccount = count
		magicZ = Zs[count]
		break
	count +=1

incrosstotalmasses = []
incrosscumtotalmasses = []
for shell in InCross_Shells:
	cumtotalmass = []
	if (len(shell) > 1): 
		thetotalmass = np.trapz(shell, Ts, 0.1)
		count = 0 
		while (count < len(shell)):
			cumtotalmass.append(np.trapz(shell[0:count+1], Ts[0:count+1], 0.1))
			count+=1
	else: thetotalmass = 0
	cumtotalmass = np.array(cumtotalmass)
	incrosstotalmasses.append(thetotalmass)
	incrosscumtotalmasses.append(cumtotalmass)

intotalmasses = []
incumtotalmasses = []
for shell in InFlow_Shells:
	cumtotalmass = []
	if (len(shell) > 1): 
		thetotalmass = np.trapz(shell, Ts, 0.1)
		count = 0 
		while (count < len(shell)):
			cumtotalmass.append(np.trapz(shell[0:count+1], Ts[0:count+1], 0.1))
			count+=1
	else: thetotalmass = 0
	cumtotalmass = np.array(cumtotalmass)
	intotalmasses.append(thetotalmass)
	incumtotalmasses.append(cumtotalmass)



#for shell in Cross_Shells:
#	if (len(shell) > 1): 
#		thecrossmass = np.trapz(shell, Ts, 0.1)
#	else: thecrossmass = 0
#	crossmasses.append(thecrossmass)

crossmasses = [] 
crosscumtotalmasses = []
for shell in Cross_Shells:
	cumtotalmass = []
	if (len(shell) > 1): 
		thecrossmass = np.trapz(shell, Ts, 0.1)
		count = 0 
		while (count < len(shell)):
			cumtotalmass.append(np.trapz(shell[0:count+1], Ts[0:count+1], 0.1))
			count+=1
	else: thetotalmass = 0
	cumtotalmass = np.array(cumtotalmass)
	crossmasses.append(thecrossmass)
	crosscumtotalmasses.append(cumtotalmass)




ML_from_total = total_expelled/total_formed
ML_from_total_cross = total_expelled_cross/total_formed

SFRs = interp1d(Ts, SFR1, kind='linear')
SFRsOuter = interp1d(Ts, SFRouter, kind='linear')
Mstar_t = interp1d(Ts, Mstar, kind='linear')
Outflows_Inner = interp1d(Ts, Flow_in,  kind='linear')
Inflows_Inner = interp1d(Ts, Flow_in2,  kind='linear')

Outflows_Outer = interp1d(Ts, Flow_out, kind='linear')
Cross_Inner = interp1d(Ts, Cross_in, kind='linear')
Cross_Outer = interp1d(Ts, Cross_out, kind='linear')

Flow_Shells_Int = [] 

for shell in Flow_Shells:
	fs = interp1d(Ts, shell, kind='linear')
	Flow_Shells_Int.append(fs)

Cross_Shells_Int = []
for shell in Cross_Shells:
	fs = interp1d(Ts, shell, kind='linear')
	Cross_Shells_Int.append(fs)

#New_Outflow = interp1d(Ts, NewOutflowRate, kind='linear')


Inflows_Inner = interp1d(Ts, Flow_in2,  kind='linear')
Inflows_Outer = interp1d(Ts, Flow_out2, kind='linear')

count = 0
Tbegin = Ts[0]
Tend = Ts[len(Ts)-1]


#surface density stuff
f = open(rhalf_finname)
dars = np.loadtxt(f, ndmin=2)

f.close()

dtSF = dars[:,2][cuts]
rhalf_densegas = dars[:,11][cuts]
mhalf_densegas = dars[:,13][cuts]

rhalf_inner_densegas = dars[:,19][cuts]
mhalf_inner_densegas = dars[:,21][cuts]

rhalf_youngstars = dars[:,27][cuts]
mhalf_youngstars = dars[:,29][cuts]

rhalf_allstars = dars[:,35][cuts]
mhalf_allstars = dars[:,37][cuts]

rhalf_densegas_phys = rhalf_densegas * As /little_h
rhalf_inner_densegas_phys = rhalf_inner_densegas * As /little_h
rhalf_youngstars_phys = rhalf_youngstars * As / little_h
rhalf_allstars_phys =  rhalf_allstars * As / little_h

mhalf_densegas *= 1e10/little_h
mhalf_inner_densegas *= 1e10 /little_h
mhalf_youngstars *= 1e10 / little_h
mhalf_allstars *= 1e10 / little_h

youngstar_half_SFR = mhalf_youngstars / dtSF
totalSFR = youngstar_half_SFR*2.0

surfden_densegas = mhalf_densegas / (rhalf_densegas_phys * rhalf_densegas_phys * 1000 * 1000  * 3.14159)
surfden_inner_densegas = mhalf_inner_densegas / (rhalf_inner_densegas_phys * rhalf_inner_densegas_phys * 1000 * 1000  * 3.14159)
surfden_youngstars_SFR = youngstar_half_SFR / (rhalf_youngstars_phys * rhalf_youngstars_phys  * 3.14159)

realcut = np.isfinite(surfden_youngstars_SFR)
fakecut = np.invert(realcut)
surfden_youngstars_SFR[fakecut] = 0.0

surfden_youngstars_SFRs = interp1d(Ts, surfden_youngstars_SFR,  kind='linear')
weirdquant =  interp1d(Ts, surfden_youngstars_SFR*SFR1,  kind='linear')

f = open(pot_finname)
dars = np.loadtxt(f)
d_pot = dars[:,3][cuts] * As
M_pot = dars[:,7][cuts]
The_pot = M_pot/d_pot

realcut = np.isfinite(The_pot)
fakecut = np.invert(realcut)
The_pot[fakecut] = 0.0

weirdpotquant = interp1d(Ts, The_pot*SFR1,  kind='linear')


Mformed  = 0
Mexpelled = 0
Mcrossed = 0
Mnew = 0
AvgSFR =0 
AvgSFR_for_delay = 0
AvgOutflowInner = 0
AvgOutflowOuter = 0
AvgInflowInner = 0
AvgInflowOuter = 0
AvgNewOutflows = 0

Mformed_ar = []
Mejected_ar = []
Mgained_ar = []

interval_M_ar = []
interval_R_ar = []
interval_Z_ar =[]
interval_T_ar = []
interval_Mgas_ar = []
interval_Mstar_ar = []
interval_Inflow_start_ar = []
interval_vcirc_ar = []
interval_max_sigSFR_ar = []
interval_max_SFR_ar =[]
interval_max_Outflow_ar = []
interval_max_Inflow_ar = []
interval_SF_start_ar = []
interval_Outflow_start_ar = []
interval_SF_end_ar = []
interval_Outflow_end_ar = []
interval_start_ar =[]
interval_average_sigSFR_ar = []
interval_average_pot_ar = []
episode_adjustment_ar = []
zero_ar = []
interval_count_ar = []
tstep = 1

interval_count = 0

if (interval <= SFR_smoothing):
	print 'interval must be larger than SFR_smoothing!'
	sys.exit()

the_time = Tbegin
print Tbegin, Tend
SFRprev = 0
 
#delay = 60   - length of SF interval now default_len_SF
#interval = 60  - length of outflw itnerval now default_len_OUT
# stagger - variable 60 for hiz, 90 for medz 120 for loz leave as stagger
default_len_SF = 60
default_len_OUT = 60


def integrate_quant(start_time, end_time, Func, dt):
	total = 0
	time = start_time
	while (time <= end_time and time <Tend):
		total += dt*Func(time)
		time+=dt
	return total


T_practical_end = Tend - max(default_len_SF, (default_len_OUT+stagger) ) # after this, the episode cant meet minimum criteria for length

pastSFR = 0
while (the_time < T_practical_end):
	if (paramfac_flies):
		paramfac = paramfac_fly(Mstar_t(the_time))
		SFR_1stpar = 1.0/paramfac
		SFR_2ndpar = 1.0
		SFR_3rdpar = 5.0/paramfac

		#for slack
		SFR_4thpar = 1.0/paramfac
		SFR_5thpar = 1.0
		SFR_6rhpar = 2.0/paramfac

	episode_criterion =  (((SFRs(the_time) > SFR_1stpar) and (SFRs(the_time) > SFRprev*SFR_2ndpar)) or SFRs(the_time )> SFR_3rdpar) # conditions to start interval 
	if (episode_criterion):
		Flow_total_interval = np.zeros(10)
		interval_count += 1 
		interval_count_ar.append(interval_count)
		print '###############################################'
		print 'beginning interval ', interval_count, the_time, SFRs(the_time),  SFR_1stpar,  SFR_3rdpar
		interval_start_time = the_time #clock start time 
		
		inflow_start_time = max([ (the_time-default_len_OUT), Tbegin])
		
		trial_SF_interval_end = the_time + default_len_OUT

		
		
		trial_SF_interval_end = the_time + default_len_OUT		
		trial_outflow_interval_end = the_time + stagger + default_len_OUT
		trial_outflow_interval_duration = trial_outflow_interval_end - the_time - stagger
		trial_SF_interval_duration = trial_SF_interval_end - the_time

		#trial_interval_mass = integrate_quant(interval_start_time, trial_SF_interval_end, SFRs, tstep)
		
		finSFR = SFRs(trial_SF_interval_end)
		#the new way to do past SFR
		pastSFR = SFRs(trial_SF_interval_end-SFR_smoothing)
		recentSFR = SFRs(trial_SF_interval_end-tstep)
		
		continue_SF_trial = (trial_outflow_interval_end < Tend and ((finSFR > SFR_4thpar and finSFR > SFR_5thpar * pastSFR) or finSFR>recentSFR or finSFR > SFR_6rhpar) and (trial_SF_interval_duration < max_interval))
		
		time_addition = 0
		
		
		while (continue_SF_trial):
			time_addition+=tstep
			trial_SF_interval_end = the_time + default_len_SF + time_addition
			trial_outflow_interval_end = the_time + stagger + default_len_OUT + time_addition
			trial_outflow_interval_duration = trial_outflow_interval_end - the_time - stagger 
			trial_SF_interval_duration = trial_SF_interval_end - the_time
			finSFR = SFRs(trial_SF_interval_end)
			#new way
			pastSFR = SFRs(trial_SF_interval_end-SFR_smoothing)
			recentSFR = SFRs(trial_SF_interval_end-tstep)
			#print 'diagnostic: ',trial_SF_interval_end, finSFR, SFR_4thpar,   SFR_5thpar,  pastSFR, recentSFR, SFR_6rhpar, time_addition
 			continue_SF_trial = (trial_outflow_interval_end < Tend and ((finSFR > SFR_4thpar and finSFR > SFR_5thpar * pastSFR) or finSFR>recentSFR or finSFR > SFR_6rhpar) and (trial_SF_interval_duration < max_interval))

			# the old way: uses time_addition instead of duration
 			#continue_SF_trial = (trial_outflow_interval_end < Tend and ((finSFR > SFR_4thpar and finSFR > SFR_5thpar * pastSFR) or finSFR>recentSFR or finSFR > SFR_6rhpar) and (time_addition < max_interval))

#while (the_time < Tend and ((SFRs(the_time) > SFR_4thpar and SFRs(the_time) > SFR_5thpar * SFRprev) or SFRs(the_time)>SFRs(the_time-tstep) or SFRs(the_time) > SFR_6rhpar) and (the_delay_slack < max_interval)):
	
		the_SF_interval_end = trial_SF_interval_end
		the_outflow_interval_end = trial_outflow_interval_end
		the_SF_interval_duration = trial_SF_interval_duration
		outflow_start_time = the_time + stagger #should be equal to: the_outflow_interval_end - time_addition 
		
		print 'inflow started at ',inflow_start_time
		print 'episode started at ',interval_start_time
		print 'SF integration ended at ',the_SF_interval_end
		print 'episode ended at ',the_outflow_interval_end
		
		
		
		if (interactive_mode):
			#plot intervals
			Tcut1 = Ts > 0.9 * interval_start_time
			Tcut2 = Ts < 1.1 * the_outflow_interval_end
			Tcut = Tcut1*Tcut2
			fig1 = plt.figure(figsize=(10,8))
			plt.plot(Ts[Tcut], SFR_boost*SFR1[Tcut], color='k')
			plt.plot(Ts[Tcut], Flow_in[Tcut], '--r')
			plt.plot(Ts[Tcut], Flow_in2[Tcut], '--b')
			plt.plot(Ts[Tcut], surfden_youngstars_SFR[Tcut]*surf_den_boost, ':m')

			plt.plot(Ts[Tcut], Flow_out[Tcut], ':r')
			plt.plot([interval_start_time], [0], 'vr', ms=13)
			plt.plot([outflow_start_time], [0], 'vg', ms=13)
			plt.plot([the_outflow_interval_end], [0], '^g', ms=13)
			plt.plot([the_SF_interval_end], [0], '^r', ms=13)
			fig2 = plt.figure(figsize=(10,8))
			plt.plot(Ts, SFR_boost*SFR1, color='k')
			plt.plot(Ts, Flow_in, '--r')
			plt.plot(Ts, Flow_in2, '--b')
			plt.plot(Ts, surfden_youngstars_SFR*surf_den_boost, ':m')
			plt.plot(Ts, Flow_out, ':r')
			plt.plot([interval_start_time], [0], 'vr', ms=13)
			plt.plot([outflow_start_time], [0], 'vg', ms=13)
			plt.plot([the_outflow_interval_end], [0], '^g', ms=13)
			plt.plot([the_SF_interval_end], [0], '^r', ms=13)		
			plt.show()

		
			#ask for confirmation
			override = raw_input('do you agree?')
			if (override == 'n'):
				print 'enterring correction mode '
				while (True):
					print 'input change in interval'
					delta = float(raw_input('delta t(Myr) = '))
					the_SF_interval_end = trial_SF_interval_end + delta
					the_outflow_interval_end = trial_outflow_interval_end + delta
					the_SF_interval_duration = the_SF_interval_end - the_time
					outflow_start_time = the_time + stagger 
	
					print 'trial '
					print 'episode started at ',interval_start_time
					print 'SF integration ended at ',the_SF_interval_end
					print 'episode ended at ',the_outflow_interval_end
	
					#plot intervals
					Tcut1 = Ts > 0.9 * interval_start_time
					Tcut2 = Ts < 1.1 * the_outflow_interval_end
					Tcut = Tcut1*Tcut2
					fig1 = plt.figure(figsize=(10,8))
					plt.plot(Ts[Tcut], SFR_boost*SFR1[Tcut], color='k')
					plt.plot(Ts[Tcut], Flow_in[Tcut], '--r')
					plt.plot(Ts[Tcut], Flow_in2[Tcut], '--b')
					plt.plot(Ts[Tcut], surfden_youngstars_SFR[Tcut]*surf_den_boost, ':m')

					plt.plot(Ts[Tcut], Flow_out[Tcut], ':r')
					plt.plot([interval_start_time], [0], 'vr', ms=13)
					plt.plot([outflow_start_time], [0], 'vg', ms=13)
					plt.plot([the_outflow_interval_end], [0], '^g', ms=13)
					plt.plot([the_SF_interval_end], [0], '^r', ms=13)
					fig2 = plt.figure(figsize=(10,8))
					plt.plot(Ts, SFR_boost*SFR1, color='k')
					plt.plot(Ts, surfden_youngstars_SFR*surf_den_boost, ':m')

					plt.plot(Ts, Flow_in, '--r')
					plt.plot(Ts, Flow_in2, '--b')
					plt.plot(Ts, Flow_out, ':r')
					plt.plot([interval_start_time], [0], 'vr', ms=13)
					plt.plot([outflow_start_time], [0], 'vg', ms=13)
					plt.plot([the_outflow_interval_end], [0], '^g', ms=13)
					plt.plot([the_SF_interval_end], [0], '^r', ms=13)		
					plt.show()		
				
					override = raw_input( 'do you like it now?')
					if (override == 'y'):
						episode_adjustment_ar.append(delta)
						break
			else:
				episode_adjustment_ar.append(0.0)
		else:
			episode_adjustment_ar.append(0.0)
		#integrate masses
		interval_Mformed = integrate_quant(interval_start_time, the_SF_interval_end, SFRs, tstep)
		interval_Mejected = integrate_quant(outflow_start_time, the_outflow_interval_end, Outflows_Inner, tstep)		
		interval_average_sigSFR = integrate_quant(interval_start_time, the_SF_interval_end, weirdquant, tstep) / interval_Mformed
		interval_average_pot = integrate_quant(interval_start_time, the_SF_interval_end, weirdpotquant, tstep) / interval_Mformed
		
		interval_Mgained =integrate_quant(inflow_start_time, the_SF_interval_end, Inflows_Inner, tstep)	

		
		Tcut1 = Ts > 1.0 * interval_start_time
		Tcut2 = Ts < 1.0 * the_outflow_interval_end		
		Tcut = Tcut1*Tcut2

		interval_median = len(Ts[Tcut])/2 
		interval_median_M = TotalMass[Tcut][interval_median]
		interval_median_R = Rvirused[Tcut][interval_median]
		interval_median_Z = Zs[Tcut][interval_median]
		interval_median_T = Ts[Tcut][interval_median]
		interval_median_Mgas = Mgas[Tcut][interval_median]
		interval_median_Mstar = Mstar[Tcut][interval_median]
		interval_median_vcirc = Vcirc[Tcut][interval_median]
		interval_median_vcirc *= (1.0 + interval_median_Z)**0.5
		Xcut = (SFR1>0) * np.isfinite(surfden_youngstars_SFR)
		if (len(surfden_youngstars_SFR[Tcut*Xcut]) > 0):
			interval_max_sigSFR = max(surfden_youngstars_SFR[Tcut*Xcut])
		else: 
			interval_max_sigSFR = 0
		interval_max_SFR = max(SFR1[Tcut])
		
		TcutO1 = Ts > 1.0 * outflow_start_time
		TcutO2 = Ts < 1.0 * the_outflow_interval_end		
		TcutO = TcutO1*TcutO2
		if (len(Flow_in[TcutO])>1):
			interval_max_Outflow = max(Flow_in[TcutO])
		else:
			interval_max_Outflow = 0
		
		TcutI1 = Ts > 1.0 * inflow_start_time
		TcutI2 = Ts < 1.0 * the_SF_interval_end		
		TcutI = TcutI1*TcutI2
		#interval_max_Inflow = min(Flow_in2[TcutI])

		if (len(Flow_in2[TcutI])>1):
			interval_max_Inflow = min(Flow_in2[TcutI])
		else:
			interval_max_Inflow = 0

		interval_M_ar.append(interval_median_M)
		interval_R_ar.append(interval_median_R)
		interval_Z_ar.append(interval_median_Z)
		interval_T_ar.append(interval_median_T)
		interval_Mgas_ar.append(interval_median_Mgas)
		interval_Mstar_ar.append(interval_median_Mstar)
		interval_vcirc_ar.append(interval_median_vcirc)
		interval_max_sigSFR_ar.append(interval_max_sigSFR)
		interval_max_SFR_ar.append(interval_max_SFR)
		interval_max_Outflow_ar.append(interval_max_Outflow)
		interval_max_Inflow_ar.append(interval_max_Inflow)
		interval_average_sigSFR_ar.append(interval_average_sigSFR)
		interval_average_pot_ar.append(interval_average_pot)
		Mformed_ar.append(interval_Mformed)
		Mejected_ar.append(interval_Mejected)
		Mgained_ar.append(interval_Mgained)
		
		interval_Inflow_start_ar.append(inflow_start_time)
		interval_SF_start_ar.append(interval_start_time)
		interval_Outflow_start_ar.append(outflow_start_time)
		interval_SF_end_ar.append(the_SF_interval_end)
		interval_Outflow_end_ar.append(the_outflow_interval_end)
		
		
		the_time = the_SF_interval_end
	else:
		the_time+= tstep
		
		

Mformed_ar = np.array(Mformed_ar)
Mejected_ar = np.array(Mejected_ar)
Mgained_ar = np.array(Mgained_ar)
interval_M_ar = np.array(interval_M_ar)
interval_R_ar = np.array(interval_R_ar)
interval_Z_ar = np.array(interval_Z_ar)
interval_T_ar = np.array(interval_T_ar)
interval_Mgas_ar = np.array(interval_Mgas_ar)
interval_Mstar_ar = np.array(interval_Mstar_ar)
interval_vcirc_ar = np.array(interval_vcirc_ar)
interval_max_sigSFR_ar = np.array(interval_max_sigSFR_ar)
interval_max_SFR_ar = np.array(interval_max_SFR_ar)
interval_max_Outflow_ar = np.array(interval_max_Outflow_ar) 
interval_max_Inflow_ar = np.array(interval_max_Inflow_ar)
interval_SF_start_ar = np.array(interval_SF_start_ar)
interval_Outflow_start_ar = np.array(interval_Outflow_start_ar)
interval_SF_end_ar = np.array(interval_SF_end_ar)
interval_Outflow_end_ar = np.array(interval_Outflow_end_ar)
interval_count_ar = np.array(interval_count_ar)
interval_average_sigSFR_ar = np.array(interval_average_sigSFR_ar)
interval_average_pot_ar = np.array(interval_average_pot_ar)
episode_adjustment_ar = np.array(episode_adjustment_ar)

eta_ar = Mejected_ar/Mformed_ar
print eta_ar
Total_ejected = np.sum(Mejected_ar)
Total_formed = np.sum(Mformed_ar)
episode_eject_fraction = Mejected_ar/Total_ejected
eta_by_fraction = np.sum(eta_ar * episode_eject_fraction)
episode_form_fraction = Mformed_ar/Total_formed
eta_by_formed = np.sum(eta_ar * episode_form_fraction)

print 'eta by outflow fraction ',eta_by_fraction
print 'eta by SF fraction ',eta_by_formed

masscolumns = np.column_stack([interval_count_ar, interval_Z_ar, interval_T_ar, Mformed_ar, Mejected_ar, interval_M_ar, interval_R_ar, interval_Mgas_ar, interval_Mstar_ar,interval_vcirc_ar, interval_max_sigSFR_ar,  interval_max_SFR_ar, interval_max_Outflow_ar, interval_SF_start_ar,interval_Outflow_start_ar, interval_SF_end_ar,  interval_Outflow_end_ar, eta_ar, episode_eject_fraction,  episode_form_fraction, interval_average_sigSFR_ar, interval_average_pot_ar,  Mgained_ar, interval_Inflow_start_ar, interval_max_Inflow_ar, episode_adjustment_ar])

jorb = np.column_stack([surfden_youngstars_SFR, SFR1])
print 'jorb ',jorb

foutname = ('./' + str(today) + 'episodes'+ShellUseInner+ ShellUseOuter +'N'+haloN+'.txt')
np.savetxt(foutname, masscolumns, fmt='%1.6g', header='0 N interval, 1 interval_Z_ar, 2 interval_T_ar, 3 Mformed_ar, 4 Mejected_ar, 5 interval_M_ar, 6 interval_R_ar, 7 interval_Mgas_ar, 8 interval_Mstar_ar, 9 interval_vcirc_ar, 10 interval_max_sigSFR_ar, 11 interval_max_SFR_ar, 12 interval_max_Outflow_ar,  13 interval_SF_start_ar, 14 interval_Outflow_start_ar, 15 interval_SF_end_ar, 16 interval_Outflow_end_ar, 17 eta_ar, 18 episode_eject_fraction, 19 episode_form_fraction, 20 interval_average_sigSFR_ar 21 interval_average_pot_ar 22 Mgained_ar 23 interval_Inflow_start_ar 24 interval_max_Outflow_ar 25 episode_adjustment_ar')

print 'finished'
