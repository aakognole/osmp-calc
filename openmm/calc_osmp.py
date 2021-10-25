import numpy as np
import MDAnalysis as mda
from sys import argv

mol1, mol2, conc = str(argv[1]), str(argv[2]), str(argv[3])
sysname = mol1+'_'+mol2+'_'+conc

f = open('../atom_selection.txt','r')
sele = str(f.read())
f.close()

boxx, boxy = 48.0, 48.0; area = boxx * boxy

k = 10 # kcal/mol/Å^2
exp = 1 
factor = ( 1000 * 4.184 ) / ( 6.0221367 * ( 10**(23) ) ) / ( 10**(-30) )

upperwall, lowerwall = 72.0, 24.0

psf = '../'+sysname+'.psf'
dcd = sysname+'.1.dcd'
u = mda.Universe(psf,dcd)
atms = u.select_atoms('(%s) and resid 1'%(sele))

file1 = open('osmp.'+sysname+'.dat', 'w')
file1.write('# %15s %15s\n'%('time','instant_osmp'))

sel1 = u.select_atoms('resname %s and name %s'%(mol1.upper(),atms.atoms[0].name))
sel2 = u.select_atoms('resname %s and name %s'%(mol2.upper(),atms.atoms[1].name))

nmol1 = sel1.n_atoms
nmol2 = sel2.n_atoms

ftot = 0

for ts in u.trajectory:
    fts = 0
    for r1 in sel1.resids:
        sel1_1 = sel1.select_atoms('resid %s'%(r1))
        if sel1_1.positions[0][2] > 72.0 :
            dz = abs(sel1_1.positions[0][2]-72.0)
            fts = fts + ( k * ( dz ** exp ) )
            ftot = ftot + ( k * ( dz ** exp ) )
        elif sel1_1.positions[0][2] < 24.0 :
            dz = abs(sel1_1.positions[0][2]-24.0)
            fts = fts + ( k * ( dz ** exp ) )
            ftot = ftot + ( k * ( dz ** exp ) )
    for r2 in sel2.resids:
        sel2_1 = sel2.select_atoms('resid %s'%(r2))
        if sel2_1.positions[0][2] > 72.0 :
            dz = abs(sel2_1.positions[0][2]-72.0)
            fts = fts + ( k * ( dz ** exp ) )
            ftot = ftot + ( k * ( dz ** exp ) )
        elif sel2_1.positions[0][2] < 24.0 :
            dz = abs(sel2_1.positions[0][2]-24.0)
            fts = fts + ( k * ( dz ** exp ) )
            ftot = ftot + ( k * ( dz ** exp ) )
    time = ts.frame / 1000
    force = fts / 2 / area    
    pressure = force * factor
    file1.write('  %15.5f %15.5f\n'%(time,pressure))

force = ftot / 2 / area / u.trajectory.n_frames
pressure = force * factor
file1.write('# Average Osmotic Pressure = %.3f bar'%(pressure))
file1.close()

exit()
