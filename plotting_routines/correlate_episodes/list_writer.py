import numpy as np
import matplotlib.pyplot as plt
from datetime import date
import sys
import os
today = str(date.today())
print today

if (len(sys.argv) < 2):
	print 'syntax  blah.py  relevant_snaps_file'
	sys.exit()

finname1 = str(sys.argv[1])
f = open(finname1)
dars = np.loadtxt(f)
f.close()


epN = dars[:,0].astype(int)
zs = dars[:,1].astype(str)
Nstart = dars[:,2].astype(int)
Nend =dars[:,3].astype(int)

print zs
foutname = today+'list_of_eps.txt'
#f = open(foutname, 'a')

count = 0
while (count < len(epN)):
	myN = Nstart[count]
	myep = str(epN[count])
	mycommand = 'echo \# z='+zs[count]+' episode number '+myep+' >> '+foutname
	#f.write('#z='+zs[count]+' episode number '+myep)
	os.system(mycommand)
	mycommand2 = 'echo \# z='+zs[count]+' episode number '+myep
	os.system(mycommand2)
	while (myN <= Nend[count]):
		mycommand = 'ls EscSurf/*'+str(myN)+'N*.txt >> '+foutname
		os.system(mycommand)
		myN+=1
	count+=1
mycommand = 'echo \# done >> '+foutname
os.system(mycommand)
