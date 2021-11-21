import numpy as np
import MDAnalysis as mda
from sys import argv
from os.path import exists as file_exists

mol1, mol2, conc = str(argv[1]), str(argv[2]), str(argv[3])
sysname = mol1+'_'+mol2+'_'+conc

f = open('../atom_selection.txt','r')
sele = str(f.read())
f.close()

boxx, boxy = 48.0, 48.0; area = boxx * boxy

k = 1 # kcal/mol/Å^2 ( 1 kcal/mol/Å^2 = 418.4 kJ/mol/nm^2 )
exp = 2 
factor = ( 1000 * 4.184 ) / ( 6.0221367 * ( 10**(23) ) ) / ( 10**(-30) )

upperwall, lowerwall = 72.0, 24.0
ftot, itot = 0, 0

for run in range(1,4):
    psf = '../'+sysname+'.c36.nowat.psf'
    dcd = sysname+'.'+str(run)+'.nowat.dcd'
    u = mda.Universe(psf,dcd)
    atms = u.select_atoms('(%s) and resid 1'%(sele))
    
    file1 = open('osmp.'+sysname+'.'+str(run)+'.dat', 'w')
    file1.write('# %15s %15s\n'%('time','instant_osmp'))
    
    sel1 = u.select_atoms('resname %s and name %s'%(mol1.upper(),atms.atoms[0].name))
    sel2 = u.select_atoms('resname %s and name %s'%(mol2.upper(),atms.atoms[1].name))
    
    nmol1 = sel1.n_atoms
    nmol2 = sel2.n_atoms
    
    for ts in u.trajectory:
        fts = 0
        for r1 in sel1.resids:
            sel1_1 = sel1.select_atoms('resid %s'%(r1))
            if abs(sel1_1.positions[0][2]) > upperwall :
                dz = abs(abs(sel1_1.positions[0][2])-upperwall)
                fts = fts + ( k * ( dz ** exp ) )
                ftot = ftot + ( k * ( dz ** exp ) )
            elif abs(sel1_1.positions[0][2]) < lowerwall :
                dz = abs(abs(sel1_1.positions[0][2])-lowerwall)
                fts = fts + ( k * ( dz ** exp ) )
                ftot = ftot + ( k * ( dz ** exp ) )
        for r2 in sel2.resids:
            sel2_1 = sel2.select_atoms('resid %s'%(r2))
            if abs(sel2_1.positions[0][2]) > upperwall :
                dz = abs(abs(sel2_1.positions[0][2])-upperwall)
                fts = fts + ( k * ( dz ** exp ) )
                ftot = ftot + ( k * ( dz ** exp ) )
            elif abs(sel2_1.positions[0][2]) < lowerwall :
                dz = abs(abs(sel2_1.positions[0][2])-lowerwall)
                fts = fts + ( k * ( dz ** exp ) )
                ftot = ftot + ( k * ( dz ** exp ) )
        time = ts.frame / 2000
        force = fts / 2 / area    
        pressure = force * factor
        file1.write('  %15.5f %15.5f\n'%(time,pressure))
        itot += 1
    file1.close()

force = ftot / 2 / area / itot
pressure = force * factor / 100000.0
print('# Average Osmotic Pressure = %.3f bar'%(pressure))

exit()
