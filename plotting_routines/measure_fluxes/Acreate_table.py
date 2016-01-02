import sys
import os

if (len(sys.argv) < 4 ):
    print 'syntax: blah.py haloN snapmin snapmax [rootstring]'
    sys.exit()

haloN = str(sys.argv[1])
snapmin = int(sys.argv[2])
snapmax = int(sys.argv[3])
rootstring = 'outflows_90514/2014-09-05histogram_mass'
rootstring = str(sys.argv[4])
name = str(sys.argv[5])

count = snapmin

while(count <= snapmax):
	if (count < 10):
		Nsnapstring = '00'+str(count)
	elif (count < 100):
		Nsnapstring = '0'+str(count)
	else:
		Nsnapstring = str(count)
	print 'doing ',count
        bigstring = rootstring+Nsnapstring+'N'+haloN+'.txt'
	mycommand = 'python '+name+'/Acreate_inner_histo.py '+Nsnapstring+' '+haloN+' '+bigstring+' '+name
	os.system(mycommand)
	count = count +1
