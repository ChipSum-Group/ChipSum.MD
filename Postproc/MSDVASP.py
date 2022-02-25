#!usr/bin/env python3
import numpy

###############################################################################################
#            This code is to calculate the MSD from VASP MD output file XDATCAR               #
#                                                                                             #
# Usage: python3 MSDVASP.py.                                                                  #
# This code need the XDATCAR in the directory.                                                #
# The variable you need modify is 'potim' and 'skip',note that the latter is the skip of the     #
# configurarions in XDATCAR, and the former is the variable of 'POTIM' set in INCAR.          #
###############################################################################################

potim = 0.0005  #time step ps. You need modify
skip = 1  #skip of MD steps. You need modify

pos = []  #atoms position
msd = []  #MSD data
msd_first_line = ['t/ps'] 
dt = potim * skip
out_file_name = {}   #position file names to be writen
n_files = {}   
msds = {}
fns = []    
po = []
system = ''
with open('XDATCAR') as file:
    for line in file:
        ll = list(line.strip(' ').split())
        po.append(ll)
    
Ntype = len(po[5])                      #numbers of atom types  
cell = [[float(po[2][0]),float(po[2][1]),float(po[2][2])],
        [float(po[3][0]),float(po[3][1]),float(po[3][2])],
        [float(po[4][0]),float(po[4][1]),float(po[4][2])]]      #cell paramater
N = sum(int(po[6][k]) for k in range(Ntype))  #number of atoms in the system
confs = int(len(po[7:])/(N+1))                      #numbers of configurations
atom = po[5][:]                          #atoms figures
natom = []                        #numbers of each type of atom
print(Ntype);print(cell);print(N);print(confs)

for j in range(Ntype):
    natom.append(int(po[6][j]))                  #numbers of each type of atom
    system = system+po[5][j]+'_'+str(po[6][j])+' '
    file_name = (po[5][j]+'.pos'); fns.append(file_name)
    out_file_name[file_name] = numpy.empty((confs,natom[j],3))
    n_files[file_name] = 'file_'+po[5][j]
    msd_first_line.append('        MSD_'+po[5][j])
msd.append(msd_first_line)

print('MD system is: '+system)
print('Cell parameter is:)
print(cell)
print('Number of atom types is: ' + str(Ntype))
print('Number of configurations is: '+str(confs))
print('There are '+str(N)+' atoms in each configuration')

for jtem in po[7:]:
    if jtem[0] != 'Direct': 
        pos.append([float(jtem[0]), float(jtem[1]), float(jtem[2])])

apos = numpy.dot(numpy.array(pos),numpy.array(cell));print(numpy.shape(apos))

for k, v in n_files.items():
    site = fns.index(k)
    with open(k, "w") as v:
        for cn in range(0,confs):
            for ia in range(0,N):
                sle = apos[cn*N+ia,:]
                if sum(natom[:site]) <= ia < sum(natom[:site+1]):
                    out_file_name[k][cn,ia-sum(natom[:site]),:] = sle
                    v.write(atom[site]+" "+str(sle[0])+" "+str(sle[1])+" "+str(sle[2])+"\n")
                else: continue

for k in range(1,confs):
    for element in fns:
        msds[element] = []
    if k%100 == 0: print(k)
    x1m = confs-k; list1 = [dt*k]
    for x1 in range(0,x1m):
        for key in fns:
            tmp = (out_file_name[key][x1+k,:,:] - out_file_name[key][x1,:,:])**2; msds[key].append(tmp)
    for key, v in n_files.items():
        site = fns.index(key)
        list1.append(sum(msds[key]).sum()/natom[site]/len(msds[key]))
    msd.append(list1)

file_name = 'MSD.out'
with open(file_name, "w") as file_object:
    for line in msd:
        for j in range(len(line)):
            file_object.write(str(line[j])+" ")
        file_object.write("\n")
