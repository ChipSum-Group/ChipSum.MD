#!usr/bin/env python3
import numpy as np
import Eckart

###############################################################################################
#            This code is to calculate the MSD from VASP MD output file XDATCAR               #
#                                                                                             #
# Usage: python3 MSDVASP.py.                                                                  #
# This code need the XDATCAR in the directory.                                                #
# The variable you need modify is 'dt' and 'skip',note that the latter is the skip of the     #
# configurarions in XDATCAR, and the former is the variable of 'POTIM' set in INCAR.          #
###############################################################################################

potim = 0.001  #time step ps. You need modify
skip = 1  #skip of MD steps. You need modify

pos = []  #atoms position
dt = potim * skip
out_file_name = {}   #position file names to be writen
n_files, ne_files = {}, {}
fns, fns_e = [], []    
po = []
system = ''
with open('XDATCAR') as file:
    for line in file:
        ll = list(line.strip(' ').split())
        po.append(ll)
    
Ntype = len(po[5])                      #numbers of atom types  
cell = float(po[2][0])                 #cell paramater
N = sum(int(po[6][k]) for k in range(Ntype))  #number of atoms in the system
confs = int(len(po[7:])/(N+1))                      #numbers of configurations
atom = po[5][:]                          #atoms figures
natom = []                        #numbers of each type of atom
print(Ntype);print(cell);print(N);print(confs)

for j in range(Ntype):
    natom.append(int(po[6][j]))                  #numbers of each type of atom
    system = system+po[5][j]+'_'+str(po[6][j])+' '
    file_name = (po[5][j]+'.pos'); fns.append(file_name)
    file_e_name = (po[5][j]+'_e.pos'); fns_e.append(file_e_name)
    n_files[file_name] = 'file_'+po[5][j]
    ne_files[file_e_name] = 'file_e'+po[5][j]

print('MD system is: '+system)
print('Cell parameter is: '+str(cell)+' A')
print('Number of atom types is: ' + str(Ntype))
print('Number of configurations is: '+str(confs))
print('There are '+str(N)+' atoms in each configuration')

for jtem in po[7:]:
    if jtem[0] != 'Direct': 
        pos.append([float(jtem[0])*cell, float(jtem[1])*cell, float(jtem[2])*cell])

apos = np.array(pos);print(np.shape(apos))

# for k, v in n_files.items():
#     site = fns.index(k)
#     with open(k, "w") as v:
#         for cn in range(0,confs):
#             for ia in range(0,N):
#                 sle = apos[cn*N+ia,:]
#                 if sum(natom[:site]) <= ia < sum(natom[:site+1]):
#                     out_file_name[k][cn,ia-sum(natom[:site]),:] = sle
#                     v.write(atom[site]+str(ia-sum(natom[:site])+1)+" "+str(sle[0])+" "+str(sle[1])+" "+str(sle[2])+"\n")
#                 else: continue

list2 = [atom, natom]
ar_list2 = np.array(list2)
ar_list = Eckart.div(apos, N)
ar_list1 = Eckart.tran_mcenter(ar_list, ar_list2)
ju = Eckart.cal_U(ar_list1, ar_list2)
print(np.shape(ar_list1))

e_pos = []
for i in range(len(ju)):
    e_pos.append(np.dot(ar_list1[i], ju[i].T).tolist())
e_apos = np.array(e_pos)    
print(np.shape(e_apos))
for k, v in ne_files.items():
    site = fns_e.index(k)
    print(site)
    with open(k, "w") as v:
        down = sum(natom[:site])
        up = down + natom[site]
        for i in range(confs-1):
            for j in range(down, up):
                sle = e_apos[i,j,:]
                v.write(atom[site]+str(j-down+1)+" "+str(sle[0])+" "+str(sle[1])+" "+str(sle[2])+"\n")
        