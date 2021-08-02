#!usr/bin/env python3
import numpy
import math

###############################################################################################
#            This code is to calculate the MSD from QE cp-MD output file .pos                 #
#                                                                                             #
# Usage: python3 MSDQE.py.                                                                    #
# This code need the .pos file in the directory.                                              #
# The variable you need modify will be labeled below:                                         #
#    1. 'dt' is the time step in ps, which is set in varaible 'dt' in cp.in file,but the units#
#       is different. In the cp.in file, the time unit is atomic units:                       #
#                 1 atomic unit time = 2.418e-5 ps                                            #
#    2. 'skip' is the configuration interval in the .pos file.                                #
#    3. 'confs' is the total configurations in the .pos file.                                 #
#    4. 'ntype' is the atomic types in the system.                                            #
#    5. 'nelement' is the element label which must correspond to atomic positions in the .pos #
#       file.                                                                                 #
#    6. 'ntyp' is the number of each type of atom in the system.                              #
#    7. You must modify the input file 'nh432h2o.pos' in line 36.                             #
#    8. The MSD is in the 'MSD.out' file.                                                     #
###############################################################################################

dt = 0.0001209442  #time step ps
skip = 10  #skip of MD steps
confs = 5001 #numbers of configurations
ntype = 4   #numbers of atom types
nelement = ['N','Cl','O','H']
ntyp = [1,1,32,68]  #numbers of each type of atom

ntime = dt * skip
n = sum(ntyp)   #number of atoms in the system
pos = []  #atoms position
msd = [[]]*ntype  #MSD data
first_line = 't/ps'
with open('nh432h2o.pos') as file:
  for line in file:
    ll = list(line.strip(' ').split(' '))
    if len(ll)==2: continue
    po = []
    for item in ll:
      if item:
        po.append(float(item))
    pos.append(po)
apos = numpy.array(pos)
epos, ename = [], []
for e in nelement:
  epos.append(numpy.empty((confs,ntyp[nelement.index(e)],3)))
  ename.append(e+'.pos')
  first_line = first_line+'            MSD_'+e

for i in range(ntype):
  with open(ename[i], "w") as file:
    for line in range(confs):
      if i == 0:
        down = 0; up = sum(ntyp[:i+1])
      else:
        down = sum(ntyp[:i]); up = sum(ntyp[:i+1])
      epos[i][line,:,:] = apos[line*n+down:line*n+up,:]
      for sle in epos[i][line,:,:].tolist():
        file.write(str(sle[0])+" "+str(sle[1])+" "+str(sle[2])+"\n")

for k in range(1,confs):
  msdt = []
  if k%10 == 0:
    print(k)
  for i in range(ntype):
    msdtt=[]
    for x1 in range(0,confs-k):
      tmp = (epos[i][x1+k,:,:]-epos[i][x1,:,:])**2; msdtt.append(tmp)
    msdt.append(msdtt)
  lmsd=[ntime*k]
  for j in range(ntype):
    lmsd.append(sum(msdt[j]).sum()/ntyp[j]/(confs-k))
  msd.append(lmsd)

file_name = 'MSD.out'
with open(file_name, "w") as file_object:
  file_object.write(first_line+"\n")
  for line in msd:
    if line:
      for j in range(0,len(line)):
        file_object.write(str(line[j])+" ")
      file_object.write("\n")
