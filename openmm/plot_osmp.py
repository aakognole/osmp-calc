import numpy as np
import matplotlib.pyplot as plt
import glob
from sys import argv
from os.path import exists as file_exists

methods = ['drude', 'c36']
mol1, mol2 = str(argv[1]), str(argv[2])
sysname = mol1+'_'+mol2

def blockavg(x,nblocks=30):
    lblock = int(len(x)/nblocks)
    m = []
    for i in range(nblocks):
        start = i*lblock
        end = (i+1)*lblock
        m.append(np.mean(x[start:end]))
    m = np.array(m)
    return np.mean(m), np.std(m)

for method in methods:
    dirs = sorted(glob.glob('%s_at_*'%(method)))
    if len(dirs) == 0:
        continue
    print(method.upper(),':',mol1.upper(),'-',mol2.upper())
    osmp = []
    f = open('OSMP_%s_%s_%s.dat'%(mol1,mol2,method), 'w')
    f.write('# %8s %10s %10s\n'%('Conc (M)','OsmP (bar)','Error'))
    print('# %8s %10s %10s'%('Conc (M)','OsmP (bar)','Error'))
    for d in dirs:
        c = d.split("_")[2]
        r1 = np.loadtxt('%s/osmp.%s_%s_%s.1.dat'%(d,mol1,mol2,c))
        r2 = np.loadtxt('%s/osmp.%s_%s_%s.2.dat'%(d,mol1,mol2,c))
        r3 = np.loadtxt('%s/osmp.%s_%s_%s.3.dat'%(d,mol1,mol2,c))    
        r = np.concatenate((r1,r2,r3))/100000.0
        m,s = blockavg(r[:,1])
        print("%10.1f %10.3f %10.3f"%(float(c),m,s))
        f.write("%10.1f %10.3f %10.3f\n"%(float(c),m,s))
        osmp.append((float(c),m,s))
    
    osmp = np.array(osmp)
    f.close()
    # plot
    plt.figure()
    plt.title(method.upper()+': '+mol1.upper()+' - '+mol2.upper())
    plt.errorbar(osmp[:,0],osmp[:,1],yerr=osmp[:,2],marker='o',markersize=5,capsize=3)
    plt.xlabel('Concentration (M)')
    plt.ylabel('Osmotic Pressure (bar)')
    plt.tight_layout()
    plt.savefig('OSMP_%s_%s_%s.png'%(mol1,mol2,method))
    plt.close()

if file_exists('OSMP_%s_%s_drude.dat'%(mol1,mol2)) and file_exists('OSMP_%s_%s_c36.dat'%(mol1,mol2)):
    osmp_drude = np.loadtxt('OSMP_%s_%s_drude.dat'%(mol1,mol2))
    osmp_c36 = np.loadtxt('OSMP_%s_%s_c36.dat'%(mol1,mol2))
    plt.figure()
    plt.title(mol1.upper()+' - '+mol2.upper())
    plt.errorbar(osmp_drude[:,0],osmp_drude[:,1],yerr=osmp_drude[:,2],marker='o',markersize=5,capsize=3,label='drude')
    plt.errorbar(osmp_c36[:,0],osmp_c36[:,1],yerr=osmp_c36[:,2],marker='o',markersize=5,capsize=3,label='c36')
    plt.xlabel('Concentration (M)')
    plt.ylabel('Osmotic Pressure (bar)')
    plt.legend()
    plt.tight_layout()
    plt.savefig('OSMP_%s_%s_both.png'%(mol1,mol2))
    plt.close()
