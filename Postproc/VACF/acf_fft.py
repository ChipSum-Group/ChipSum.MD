#user/bin/python3
import numpy as np
import linecache
import statsmodels.tsa.api as smt
from scipy.fftpack import fft
import math, sys

#############################################################################################################
# This script is used to calculate the ACF of velocity(position) based on VASP MD trajectory file "XDATCAR".#
#Note:                                                                                                      #
# 1. The time step "dt" you should modify based on your MD simulation.                                      #
# 2. Used after "output_pos.py".                                                                            #
# 3. Variable "f_name" you could choose ".pos" file and "_e.pos" file, which is similar for "acf_file_name" #
#    and "fft_file_name".                                                                                   #
# 4. Usage: python3 acf_fft.py dev=0,1,2,3; "dev" is the i-th derivative of trajectory position.            #
#############################################################################################################

dt = 1e-15         #time steps /s
time, f = [], []
Fs = (1/dt)*1e-12    #sampling frequency /THz

dev = int(sys.argv[1].split("=")[1])
print(dev)
species = linecache.getline('XDATCAR', 6).strip(' ').split()
num_species = linecache.getline('XDATCAR', 7).strip(' ').split()

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
        v = (array1[:,i+1]-array1[:,i])
        velocity.append(v.tolist())
    return np.array(velocity).T

def cal_acf(lis):
    N = len(lis)
#   acf = smt.stattools.acf(lis[:N], nlags=N-1,unbiased=True)     #unbiased acf
    acf = smt.stattools.acf(lis[:], nlags=N-1)                    #biased acf
    return (acf/acf[0]).tolist()

def cal_fft(lis1):
    n1 = int(math.log(len(lis1),2))
    if np.power(2,n1) == len(lis1):
        N = np.power(2,n1)
        lis3 = lis1[:]
    else:
        N = np.power(2,n1+1)
        n2 = N - len(lis1)
        lis2 = [0.0]*n2
        lis3 = lis1[:] + lis2[:]
#   fft_am = np.abs(fft(lis1)[:int(N/2)]/N)   #normaled fft
    fft_am = np.abs(fft(lis3)[:int(N/2)])  #unnormaled fft
    return fft_am.tolist()

def write_result(ary, name):
    s = np.shape(ary)
    with open(name, "w") as wfe:
        for i in range(s[1]):
            for j in range(s[0]):
                wfe.write(str(ary[j,i])+' ')
            wfe.write("\n")
        
for ele in species:
    acf_vel_all, fft_acf_all = [], []
    f_name = ele + '.pos'
    position = read_pos(f_name)
    pos_div = divide(position,int(num_species[species.index(ele)]))
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
    fft_acf_all.insert(0,f)
    write_result(np.array(acf_vel_all), acf_file_name)
    write_result(np.array(fft_acf_all), fft_file_name)
    
    
    with open(ele+"_acf_ave.out", "w") as wfe:
        for i in range(len(acf_vel_ave)):
            wfe.write(str(time[i])+' '+str(acf_vel_ave[i])+'\n')
    
    with open(ele+"_fft_ave.out", "w") as wfe:
        for i in range(len(fft_acf_ave)):
            wfe.write(str(time[i])+' '+str(fft_acf_ave[i])+'\n')
