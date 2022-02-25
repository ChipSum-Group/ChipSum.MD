#user/bin/python3
#
# This script is used to calculate the Eckart Frame Algorithm

import numpy as np

periodic_mass = {'H': 1.00794, 'He': 4.002602, 'Li': 6.941, 'Be': 9.0121831, 'B': 10.811, 'C': 12.0107, 'N': 14.0067, 'O': 15.9994, 'F': 18.99840316, 'Ne': 20.1797, 'Na': 22.98976928, 'Mg': 24.305, 'Al': 26.9815385, 'Si': 28.0855, 'P': 30.973762, 'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.0983, 'Ca': 40.078, 'Sc': 44.955908, 'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938044, 'Fe': 55.845, 'Co': 58.933194, 'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.38, 'Ga': 69.723, 'Ge': 72.64, 'As': 74.921595, 'Se': 78.971, 'Br': 79.904, 'Kr': 83.798, 'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.90584, 'Zr': 91.224, 'Nb': 92.90637, 'Mo': 95.95, 'Tc': 98.9072, 'Ru': 101.07, 'Rh': 102.9055, 'Pd': 106.42, 'Ag': 107.8682, 'Cd': 112.414, 'In': 114.818, 'Sn': 118.71, 'Sb': 121.76, 'Te': 127.6, 'I': 126.90447, 'Xe': 131.293, 'Cs': 132.905452, 'Ba': 137.327, 'La': 138.90547, 'Ce': 140.116, 'Pr': 140.90766, 'Nd': 144.242, 'Pm': 144.9, 'Sm': 150.36, 'Eu': 151.964, 'Gd': 157.25, 'Tb': 158.92535, 'Dy': 162.5, 'Ho': 164.93033, 'Er': 167.259, 'Tm': 168.93422, 'Yb': 173.054, 'Lu': 174.9668, 'Hf': 178.49, 'Ta': 180.94788, 'W': 183.84, 'Re': 186.207, 'Os': 190.23, 'Ir': 192.217, 'Pt': 195.084, 'Au': 196.966569, 'Hg': 200.59, 'Tl': 204.3833, 'Pb': 207.2, 'Bi': 208.9804, 'Po': 208.9824, 'At': 209.9871, 'Rn': 222.0176, 'Fr': 223.0197, 'Ra': 226.0245, 'Ac': 227.0277, 'Th': 232.0377, 'Pa': 231.03588, 'U': 238.02891, 'Np': 237.0482, 'Pu': 239.0642, 'Am': 243.0614, 'Cm': 247.0704, 'Bk': 247.0703, 'Cf': 251.0796, 'Es': 252.083, 'Fm': 257.0591, 'Md': 258.0984, 'No': 259.101, 'Lr': 262.1097, 'Rf': 267.1218, 'Db': 268.1257, 'Sg': 269.1286, 'Bh': 274.1436, 'Hs': 277.1519, 'Mt': 278, 'Ds': 281, 'Rg': 282, 'Cn': 285, 'Nh': 284, 'Fl': 289, 'Mc': 288, 'Lv': 292}

def tran_mcenter(ary1, ary2):
    sh = np.shape(ary1)
    new_pos, m = [], []
    total_mass = 0
    for el in ary2.T:
        m = m + [periodic_mass[el[0]]]*int(el[1])
        total_mass = total_mass + periodic_mass[el[0]]*int(el[1])
    am3 = np.array([m, m, m]).T
    for i in range(sh[0]):
        mcenter = (ary1[i]*am3).sum(axis=0)/total_mass
        new_pos.append(ary1[i]-mcenter)
    return np.array(new_pos)

def cal_U(arry1, arry2):
    s = np.shape(arry1)
    U_all, mass = [], []
    for e in arry2.T:
        mass = mass + [periodic_mass[e[0]]]*int(e[1])
    for i in range(1,s[0]):
        c = np.zeros((4,4), dtype=float)
        U = np.zeros((3,3), dtype=float)
        
        t1 = ((arry1[0,:,:]-arry1[i,:,:])*(arry1[0,:,:]-arry1[i,:,:])).sum(axis=1)
        t2 = (arry1[0,:,1]+arry1[i,:,1])*(arry1[0,:,2]-arry1[i,:,2])-(arry1[0,:,2]+arry1[i,:,2])*(arry1[0,:,1]-arry1[i,:,1])
        t3 = (arry1[0,:,2]+arry1[i,:,2])*(arry1[0,:,0]-arry1[i,:,0])-(arry1[0,:,0]+arry1[i,:,0])*(arry1[0,:,2]-arry1[i,:,2])
        t4 = (arry1[0,:,0]+arry1[i,:,0])*(arry1[0,:,1]-arry1[i,:,1])-(arry1[0,:,1]+arry1[i,:,1])*(arry1[0,:,0]-arry1[i,:,0])
        t5 = (arry1[0,:,0]-arry1[i,:,0])*(arry1[0,:,0]-arry1[i,:,0])+(arry1[0,:,1]+arry1[i,:,1])*(arry1[0,:,1]+arry1[i,:,1])+(arry1[0,:,2]+arry1[i,:,2])*(arry1[0,:,2]+arry1[i,:,2])
        t6 = (arry1[0,:,0]-arry1[i,:,0])*(arry1[0,:,1]-arry1[i,:,1])-(arry1[0,:,1]+arry1[i,:,1])*(arry1[0,:,0]+arry1[i,:,0])
        t7 = (arry1[0,:,0]-arry1[i,:,0])*(arry1[0,:,2]-arry1[i,:,2])-(arry1[0,:,2]+arry1[i,:,2])*(arry1[0,:,0]+arry1[i,:,0])
        t8 = (arry1[0,:,0]+arry1[i,:,0])*(arry1[0,:,0]+arry1[i,:,0])+(arry1[0,:,1]-arry1[i,:,1])*(arry1[0,:,1]-arry1[i,:,1])+(arry1[0,:,2]+arry1[i,:,2])*(arry1[0,:,2]+arry1[i,:,2])
        t9 = (arry1[0,:,1]-arry1[i,:,1])*(arry1[0,:,2]-arry1[i,:,2])-(arry1[0,:,2]+arry1[i,:,2])*(arry1[0,:,1]+arry1[i,:,1])
        t10 = (arry1[0,:,0]+arry1[i,:,0])*(arry1[0,:,0]+arry1[i,:,0])+(arry1[0,:,1]+arry1[i,:,1])*(arry1[0,:,1]+arry1[i,:,1])+(arry1[0,:,2]-arry1[i,:,2])*(arry1[0,:,2]-arry1[i,:,2])
        
        c[0,0] = (t1*mass).sum()
        c[0,1] = (t2*mass).sum()
        c[0,2] = (t3*mass).sum()
        c[0,3] = (t4*mass).sum()
        c[1,0] = (t2*mass).sum()
        c[1,1] = (t5*mass).sum()
        c[1,2] = (t6*mass).sum()
        c[1,3] = (t7*mass).sum()
        c[2,0] = (t3*mass).sum()
        c[2,1] = (t6*mass).sum()
        c[2,2] = (t8*mass).sum()
        c[2,3] = (t9*mass).sum()
        c[3,0] = (t4*mass).sum()
        c[3,1] = (t7*mass).sum()
        c[3,2] = (t9*mass).sum()
        c[3,3] = (t10*mass).sum()
        
        value, vector = np.linalg.eig(c)
        site = value.tolist().index(min(value))
        Q = vector[site]
        
        U[0,0] = Q[0]*Q[0]+Q[1]*Q[1]-Q[2]*Q[2]-Q[3]*Q[3]
        U[0,1] = 2*(Q[1]*Q[2]+Q[0]*Q[3])
        U[0,2] = 2*(Q[1]*Q[3]-Q[0]*Q[2])
        U[1,0] = 2*(Q[1]*Q[2]-Q[0]*Q[3])
        U[1,1] = Q[0]*Q[0]-Q[1]*Q[1]+Q[2]*Q[2]-Q[3]*Q[3]
        U[1,2] = 2*(Q[1]*Q[0]+Q[2]*Q[3])
        U[2,0] = 2*(Q[1]*Q[3]+Q[0]*Q[2])
        U[2,1] = 2*(Q[2]*Q[3]-Q[1]*Q[0])
        U[2,2] = Q[0]*Q[0]-Q[1]*Q[1]-Q[2]*Q[2]+Q[3]*Q[3]
        
        U_all.append(U)
        
    return U_all

def div(arry3, n):
    lst3 = []
    for i in range(int(np.shape(arry3)[0]/n)):
        lst3.append(arry3[i:i+n].tolist())
    return np.array(lst3)


