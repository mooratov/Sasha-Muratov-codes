import sys
import os
import re

finname = 'foullist.txt'
f = open(finname)
dars = f.readlines()
p = re.compile("\d+")
q = re.compile("AHF_\D+")
for a in dars:
    jo = p.findall(a)
    bo = q.findall(a)
    print jo
    print bo
    if (len(bo)>0):
        far = 'snap'+jo[0]+'RPep.z'+jo[1]+'.'+jo[2]+'.'+bo[0]
        print 'far ',far.strip()
        print 'a ', a.strip()
        mycommand = 'cp '+a.strip()+' '+far.strip()
        print 'command ',mycommand
        os.system(mycommand)
