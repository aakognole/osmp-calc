from __future__ import print_function
import MDAnalysis as mda
from sys import argv
from os.path import exists as file_exists

mol1, mol2, conc = str(argv[1]), str(argv[2]), str(argv[3])
pdb='../'+mol1+'_'+mol2+'_'+conc+'.mini.c36.pdb'
u = mda.Universe(pdb)

if file_exists('../atom_selection.txt'):
    f = open('../atom_selection.txt','r')
    sele = str(f.read())
    f.close()
else:
    print("------")
    print("Select one atom from each of the two molecules to add restraints:")
    print("Atoms in Mol1:", u.select_atoms('resid 1 and resname %s'%(mol1.upper())).atoms.names)
    print("Atoms in Mol2:", u.select_atoms('resid 1 and resname %s'%(mol2.upper())).atoms.names)
    print("selection of atoms (e.g. 'name MG or name CLA')")
    print(">>> ")
    sele = str(input())
    f = open('../atom_selection.txt','w')
    f.write(sele)
    f.close()
print('atoms for plumed =', sele)
g1 = u.select_atoms(sele)

f = open("plumed.py","w")
f.write('def plumedscript():\n')
f.close()
f = open("plumed.py","a")
f.write('    script = """\n')
for i in g1.atoms.indices:
    f.write("POSITION ATOM=%s LABEL=d%s\n" % (i+1,i+1))
f.write('\n')
        
upper="UPPER_WALLS "

upper=upper+'ARG='
for i in g1.atoms.indices:
    upper=upper+"d%s.z" % (i+1)+','

upper=upper[:-1]+' '

upper=upper+'KAPPA='
for i in g1.atoms.indices:
    upper=upper+"4184,"

upper=upper[:-1]+' '

upper=upper+'EXP='
for i in g1.atoms.indices:
    upper=upper+"1,"

upper=upper[:-1]+' '

upper=upper+'EPS='
for i in g1.atoms.indices:
    upper=upper+"1,"

upper=upper[:-1]+' '

upper=upper+'OFFSET='
for i in g1.atoms.indices:
    upper=upper+"0,"

upper=upper[:-1]+' '

upper=upper+'AT='
for i in g1.atoms.indices:
    upper=upper+"7.2,"

upper=upper[:-1]+' LABEL=uwall'

f.write(upper)
f.write('\n')
f.write('\n')

lower='LOWER_WALLS '

lower=lower+'ARG='
for i in g1.atoms.indices:
    lower=lower+"d%s.z" % (i+1)+','

lower=lower[:-1]+' '

lower=lower+'KAPPA='
for i in g1.atoms.indices:
    lower=lower+"4184,"

lower=lower[:-1]+' '

lower=lower+'EXP='
for i in g1.atoms.indices:
    lower=lower+"1,"

lower=lower[:-1]+' '

lower=lower+'EPS='
for i in g1.atoms.indices:
    lower=lower+"1,"

lower=lower[:-1]+' '

lower=lower+'OFFSET='
for i in g1.atoms.indices:
    lower=lower+"0,"

lower=lower[:-1]+' '

lower=lower+'AT='
for i in g1.atoms.indices:
    lower=lower+"2.4,"

lower=lower[:-1]+' LABEL=lwall'

f.write(lower)
f.write('\n')
f.write('\n')
f.write('PRINT ARG=uwall.bias,lwall.bias STRIDE=500 FILE=cv.0.dat """\n')
f.write('    return script\n')
f.close()
