import sys
import os
import numpy as np

#0 for hi z,  1for med z, 2 for loz

setting = 0
prefix = '2015-05-05'


if  (len(sys.argv) < 3):
	print 'syntax  blah.py setting name'
	sys.exit()
if (len(sys.argv) > 3):
	'using extra weirdness '
	weirdness = str(sys.argv[3])
else:
	weirdness = ''

setting = int(sys.argv[1])
name = str(sys.argv[2])

if (setting == 0):
	finname = '/Users/Muratov/research/2014A/62714/'+name+'/results-hiz'+weirdness+'/good_halos.txt'
	force_zend = False
	forced_zstart = True
	outflowdir = 'outflows_50515'
	zender = 0.499
	zstarter = 4.0
	Nstarter = '90'
	Nender = ''
	z_tolerance = 1.0
	stagger = 60


if (setting == 1):
	finname = '/Users/Muratov/research/2014A/62714/'+name+'/results-loz'+weirdness+'/good_halos.txt'
	force_zend = True
	forced_zstart = False
	outflowdir = 'outflows_50515-loz'
	zender = 0.499
	zstarter = 0.5
	Nstarter = ''
	Nender = '340'
	z_tolerance = 0.75
	stagger = 90
	
if (setting == 2):
	finname = '/Users/Muratov/research/2014A/62714/'+name+'/results-loz'+weirdness+'/good_halos.txt'
	force_zend = False
	forced_zstart = True
	outflowdir = 'outflows_50515-loz'
	zender = 0.499
	zstarter = 0.5
	Nstarter = '340'
	Nender = ''
	z_tolerance = 0.25
	stagger = 120

dars = np.loadtxt(finname, ndmin=2)
haloN = dars[:,0]
zstarts = dars[:,4]
zends = dars[:,5]
Nstarts = np.array(map(str, map(int, dars[:,6])))
Nends = np.array(map(str, map(int, dars[:,7])))

print Nstarts

count = 0
if (force_zend):
	for q in zends:
		if (q < zender):
			print 'adjusting end ',count,' from ',q,' to ',zender
			zends[count] = zender
			Nends[count] = Nender
		count +=1

count = 0
if (forced_zstart):
	for q in zstarts:
		if (q > zstarter):
			print 'adjusting start ',count,' from ',q,' to ',zstarter
			zstarts[count] = zstarter
			Nstarts[count] = Nstarter
		count +=1


zint = zstarts - zends
zcut = zint >= z_tolerance

print len(zends[zcut]),' out of ',len(zends),' meet tolerance criteria'

haloN = haloN[zcut]
zstarts = zstarts[zcut]
zends = zends[zcut]
Nstarts = Nstarts[zcut]
Nends = Nends[zcut]

halo_to_do = []
for q in haloN:
	halo_to_do.append(int(q))

print 'halos on agenda ',haloN
print 'zstarts on agenda ',zstarts
print 'zends on agenda ',zends
print 'Nstarts on agenda ',Nstarts
print 'Nends on agenda ',Nends

#lim1 = '80'
#lim2 = '190'
count = 0
for haloN in halo_to_do:
	print 'doing ',haloN
	bigstring = name+'/'+outflowdir+'/'+prefix+'histogram'
	#str1 = outflowdir+'/'+prefix+'OutflowN'+str(haloN)+'.txt'
	#str2 = outflowdir+'/'+prefix+'InflowN'+str(haloN)+'.txt'
	mycommand = 'python '+name+'/Acreate_table.py '+str(haloN)+' '+Nstarts[count]+' '+Nends[count]+' '+bigstring+' '+ name
	os.system(mycommand)
	count +=1

