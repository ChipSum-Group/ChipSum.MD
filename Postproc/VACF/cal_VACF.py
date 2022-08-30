#user/bin/python3
import numpy as np
import linecache
import statsmodels.tsa.api as smt
from scipy.fftpack import fft
from scipy.optimize import curve_fit
import math, sys, Eckart

#############################################################################################################
# This script is used to calculate the ACF of velocity(position) based on VASP MD trajectory file "XDATCAR".#
#Note:                                                                                                      #
# 1. The time step "dt" you should modify based on your MD simulation.                                      #
# 2. Used after "output_pos.py".                                                                            #
# 3. The POSCAR should be in the directory.                                                                 #
# 4. Usage: python3 acf_fft.py dev=0/1/2/3 num=10000/20000/30000...                                         #
#    "dev" is the i-th derivative of trajectory position.                                                   #
#    "num" is the last num MD steps used to calculate VACF.                                                 #
#############################################################################################################

dev = int(sys.argv[1].split("=")[1])
print(dev)
num = int(sys.argv[2].split("=")[1])
species = linecache.getline('POSCAR', 6).strip(' ').split()
num_species = linecache.getline('POSCAR', 7).strip(' ').split()

time, f = [], []
paras = [["element","xc","2sd"]]
force_constant = [["element","FC_sum","2sd"]]
dt = 1e-15         #取样时间间隔s
Fs = (1/dt)*1e-12    #抽样频率THz
C = 0.06555545
periodic_mass = {'H': 1.00794, 'He': 4.002602, 'Li': 6.941, 'Be': 9.0121831, 'B': 10.811, 'C': 12.0107, 'N': 14.0067, 'O': 15.9994, 'F': 18.99840316, 'Ne': 20.1797, 'Na': 22.98976928, 'Mg': 24.305, 'Al': 26.9815385, 'Si': 28.0855, 'P': 30.973762, 'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.0983, 'Ca': 40.078, 'Sc': 44.955908, 'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938044, 'Fe': 55.845, 'Co': 58.933194, 'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.38, 'Ga': 69.723, 'Ge': 72.64, 'As': 74.921595, 'Se': 78.971, 'Br': 79.904, 'Kr': 83.798, 'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.90584, 'Zr': 91.224, 'Nb': 92.90637, 'Mo': 95.95, 'Tc': 98.9072, 'Ru': 101.07, 'Rh': 102.9055, 'Pd': 106.42, 'Ag': 107.8682, 'Cd': 112.414, 'In': 114.818, 'Sn': 118.71, 'Sb': 121.76, 'Te': 127.6, 'I': 126.90447, 'Xe': 131.293, 'Cs': 132.905452, 'Ba': 137.327, 'La': 138.90547, 'Ce': 140.116, 'Pr': 140.90766, 'Nd': 144.242, 'Pm': 144.9, 'Sm': 150.36, 'Eu': 151.964, 'Gd': 157.25, 'Tb': 158.92535, 'Dy': 162.5, 'Ho': 164.93033, 'Er': 167.259, 'Tm': 168.93422, 'Yb': 173.054, 'Lu': 174.9668, 'Hf': 178.49, 'Ta': 180.94788, 'W': 183.84, 'Re': 186.207, 'Os': 190.23, 'Ir': 192.217, 'Pt': 195.084, 'Au': 196.966569, 'Hg': 200.59, 'Tl': 204.3833, 'Pb': 207.2, 'Bi': 208.9804, 'Po': 208.9824, 'At': 209.9871, 'Rn': 222.0176, 'Fr': 223.0197, 'Ra': 226.0245, 'Ac': 227.0277, 'Th': 232.0377, 'Pa': 231.03588, 'U': 238.02891, 'Np': 237.0482, 'Pu': 239.0642, 'Am': 243.0614, 'Cm': 247.0704, 'Bk': 247.0703, 'Cf': 251.0796, 'Es': 252.083, 'Fm': 257.0591, 'Md': 258.0984, 'No': 259.101, 'Lr': 262.1097, 'Rf': 267.1218, 'Db': 268.1257, 'Sg': 269.1286, 'Bh': 274.1436, 'Hs': 277.1519, 'Mt': 278, 'Ds': 281, 'Rg': 282, 'Cn': 285, 'Nh': 284, 'Fl': 289, 'Mc': 288, 'Lv': 292}

def read_pos(file_name):
    pos =[]
    with open(file_name) as fi:
        for line in fi:
            ll = line.strip(' ').split()
            pos.append(list(map(float, ll[1:])))
    return np.array(pos)

def divide(array, n):
    result = []
    for i in range(n):
        x = array[i::n][:,0].tolist()
        y = array[i::n][:,1].tolist()
        z = array[i::n][:,2].tolist()
        result.append(x);result.append(y);result.append(z);
    return np.array(result)

def cal_vel(array1):
    velocity = []
    for i in range(np.shape(array1)[1]-1):
        v = (array1[:,i+1]-array1[:,i])*1e-10/dt
        velocity.append(v.tolist())
    return np.array(velocity).T

def cal_acf(lis):
    N = len(lis)
#   acf = smt.stattools.acf(lis[:N], nlags=N-1,unbiased=True)
    acf = smt.stattools.acf(lis[:], nlags=N-1)
    return (acf/acf[0]).tolist()

def cal_fft(lis1):
    fft_am = np.abs(fft(lis1)[:int(len(lis1)/2)])
    return fft_am.tolist()

def write_result(ary, name):
    s = np.shape(ary)
    with open(name, "w") as wfe:
        for i in range(s[1]):
            for j in range(s[0]):
                wfe.write(str(ary[j,i])+' ')
            wfe.write("\n")

def func(x,y0,A,xc,w):
    return y0+2*A/math.pi*w/(4*np.power((x-xc),2)+np.power(w,2))

for ele in species:
    mass = periodic_mass[ele]
    acf_vel_all, fft_acf_all = [], []
    f_name = ele + '.pos'
    n_ele = int(num_species[species.index(ele)])
    position = read_pos(f_name)[-(num+1)*n_ele:]
    pos_div = divide(position,n_ele)
    if dev == 0:
        vel = pos_div[:,:]
    elif dev == 1:
        vel = cal_vel(pos_div)*1e5
    elif dev == 2:
        vel1 = cal_vel(pos_div)
        vel = cal_vel(vel1)*1e5
    elif dev == 3:
        vel2 = cal_vel(pos_div)
        vel1 = cal_vel(vel2)
        vel = cal_vel(vel1)*1e5
    else:
        print("You should give the variable 'dev' as 'python3 script.py dev=0,1,2,3'")
        exit()
    row = np.shape(vel)[0]
    for i in range(row):
        acf_vel = cal_acf(vel[i,:].tolist())
        fft_acf = cal_fft(acf_vel[:])
        acf_vel_all.append(acf_vel)
        fft_acf_all.append(fft_acf)
    acf_file_name = ele + '_acf.out'
    fft_file_name = ele + '_fft.out'
    print(np.shape(fft_acf_all))
    n_data = np.shape(acf_vel_all)[1]
    f_data = np.shape(fft_acf_all)[1]
    time = [dt*i for i in range(n_data)]
    f = [Fs*i/(n_data) for i in range(int(f_data))]
    acf_vel_ave = np.average(acf_vel_all, axis=0)
    fft_acf_ave = np.average(fft_acf_all, axis=0)
    acf_vel_all.insert(0,time)
    a_acf = np.array(acf_vel_all)
    fft_acf_all.insert(0,f)
    a_fft = np.array(fft_acf_all)
    para = []
    for i in range(1,len(fft_acf_all)):
        popt, pcov = curve_fit(func, a_fft[0], a_fft[i])
        fc = C*mass*np.power(popt[2],2)
        err = C*mass*np.power(np.diag(pcov)[2],2)
        para.append([ele+str(i),fc,err])
    paras = paras + para
    for i in range(int(len(para)/3)):
        s = i*3
        force_sum = (para[s][1]+para[s+1][1]+para[s+2][1])
        sd = np.sqrt(np.power(para[s][2],2)+np.power(para[s+1][2],2)+np.power(para[s+2][2],2))
        force_constant.append([ele+str(i),force_sum,sd])
    
    write_result(a_acf, acf_file_name)
    write_result(a_fft, fft_file_name)
    
    with open(ele+"_acf_ave.out", "w") as wfe:
        for i in range(len(acf_vel_ave)):
            wfe.write(str(time[i])+' '+str(acf_vel_ave[i])+'\n')
    
    with open(ele+"_fft_ave.out", "w") as wfe:
        for i in range(len(fft_acf_ave)):
            wfe.write(str(f[i])+' '+str(fft_acf_ave[i])+'\n')

print(np.shape(paras))
print(np.shape(force_constant))

with open("force_sonst.out", "w") as peak:
    for data in paras[:]:
        peak.write(data[0]+' '+str(data[1])+' '+str(data[2])+'\n')

with open("FC_sum.out","w") as fc:
    for f in force_constant:
        fc.write(f[0]+' '+str(f[1])+' '+str(f[2])+'\n')