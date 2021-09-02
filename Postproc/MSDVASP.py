#!usr/bin/env python3
import numpy
import math

###############################################################################################
#            This code is to calculate the MSD from VASP MD output file XDATCAR               #
#                                                                                             #
# Usage: python3 MSDVASP.py.                                                                  #
# This code need the XDATCAR in the directory.                                                #
# The variable you need modify is 'dt' and 'skip',note that the latter is the skip of the     #
# configurarions in XDATCAR, and the former is the variable of 'POTIM' set in INCAR.          #
###############################################################################################

dt = 0.0005  #time step ps. You need modify
skip = 10  #skip of MD steps. You need modify

pos = []  #atoms position
msd = []  #MSD data
msd_first_line = ['t/ps'] 
ntime = dt * skip
atoms = []  #atoms figures
out_file_name = {}   #position file names to be writen
n_files = {}   
msds = {}
ofns = []    
ntyp = []  #numbers of each type of atom
po = []
system = ''
with open('XDATCAR') as file:
  for line in file:
    ll = list(line.strip(' ').split())
    po1 = []
    for item in ll:
      if item:
        po1.append(item)
    po.append(po1)
    
ntype = len(po[5])                      #numbers of atom types  
cell = float(po[2][0]); ce2 = cell**2
n = sum(int(po[6][k]) for k in range(ntype))  #number of atoms in the system
confs = int(len(po[7:])/(n+1))                      #numbers of configurations
print(ntype);print(cell);print(n);print(confs)
for j in range(len(po[5])):
  system = system+po[5][j]+'_'+str(po[6][j])+' '
  atoms.append(po[5][j])                #atoms figures
  ntyp.append(int(po[6][j]))           #numbers of each type of atom
  ofn = (po[5][j]+'.pos'); ofns.append(ofn)
  out_file_name[ofn] = numpy.empty((confs,ntyp[j],3))
  n_files[ofn] = 'file_'+po[5][j]
  msd_first_line.append('        MSD_'+po[5][j])
msd.append(msd_first_line)
print('MD system is: '+system)
print('Cell parameter is: '+str(cell)+' A')
print('Number of atom types is: ' + str(ntype))
print('Number of configurations is: '+str(confs))
print('There are '+str(n)+' atoms in each configuration')
for jtem in po[7:]:
  if len(jtem)==3 and jtem[0] != 'Direct': 
    pos.append([float(jtem[0])*cell, float(jtem[1])*cell, float(jtem[2])*cell])

apos = numpy.array(pos);print(numpy.shape(apos))

for k, v in n_files.items():
  site = ofns.index(k)
  with open(k, "w") as v:
    for line in range(0,confs):
      for il in range(0,n):
        sle = apos[line*n+il,:]
        if sum(ntyp[:site]) <= il < sum(ntyp[:site+1]):
          out_file_name[k][line,il-sum(ntyp[:site]),:] = sle
          v.write(atoms[site]+" "+str(sle[0])+" "+str(sle[1])+" "+str(sle[2])+"\n")
        else: continue

for k in range(1,confs):
  for i in ofns:
    msds[i] = []
  if k%10 == 0: print(k)
  x1m = confs-k; list1 = [ntime*k]
  for x1 in range(0,x1m):
    for key, v in n_files.items():
      tmp = (out_file_name[key][x1+k,:,:] - out_file_name[key][x1,:,:])**2; msds[key].append(tmp)
  for key, v in n_files.items():
    site = ofns.index(key)
    list1.append(sum(msds[key]).sum()/ntyp[site]/len(msds[key]))
  msd.append(list1)
# print(msd)

file_name = 'MSD.out'
with open(file_name, "w") as file_object:
  for line in msd:
    for j in range(len(line)):
      file_object.write(str(line[j])+" ")
    file_object.write("\n")

    
