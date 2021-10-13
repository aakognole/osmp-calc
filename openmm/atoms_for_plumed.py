from __future__ import print_function
import MDAnalysis as mda
from sys import argv

ion=str(argv[1])
mol=str(argv[2])
conc=str(argv[3])
pdb='../'+ion+'_'+mol+'_'+conc+'.mini.pdb'
u = mda.Universe(pdb)

#sele = str(input("selection of atoms (e.g. 'name POT or name P'):"))
sele = 'name MG or name P'
print('atoms for plumed =', sele)
g1 = u.select_atoms(sele)

#print(g1)

f = open("./template/plumed.py","w")
f.write('def plumed():\n')
f.close()
f = open("./template/plumed.py","a")
f.write('    script = """\n')
for i in g1.atoms.indices:
    f.write("POSITION ATOM=%s NOPBC LABEL=d%s\n" % (i+1,i+1))
f.write('\n')
        
#print()

upper="UPPER_WALLS "

#print('ARG', end='=')
upper=upper+'ARG='
for i in g1.atoms.indices:
    #print("d%s.z" % (i+1), end=',')
    upper=upper+"d%s.z" % (i+1)+','
#print()

upper=upper[:-1]+' '

#print()

#print('KAPPA', end='=')
upper=upper+'KAPPA='
for i in g1.atoms.indices:
    #print("150", end=',')
    upper=upper+"150,"
#print()

upper=upper[:-1]+' '

#print()

#print('EXP', end='=')
upper=upper+'EXP='
for i in g1.atoms.indices:
    #print("2", end=',')
    upper=upper+"2,"
#print()

upper=upper[:-1]+' '

#print()

#print('EPS', end='=')
upper=upper+'EPS='
for i in g1.atoms.indices:
    #print("1", end=',')
    upper=upper+"1,"
#print()

upper=upper[:-1]+' '

#print()

#print('OFFSET', end='=')
upper=upper+'OFFSET='
for i in g1.atoms.indices:
    #print("0", end=',')
    upper=upper+"0,"
#print()

upper=upper[:-1]+' '

#print()

#print('AT', end='=')
upper=upper+'AT='
for i in g1.atoms.indices:
    #print("7.2", end=',')
    upper=upper+"7.2,"
#print()

upper=upper[:-1]+' LABEL=uwall'

f.write(upper)
f.write('\n')
f.write('\n')
#print()

lower='LOWER_WALLS '

lower=lower+'ARG='
for i in g1.atoms.indices:
    #print("d%s.z" % (i+1), end=',')
    lower=lower+"d%s.z" % (i+1)+','
#print()

lower=lower[:-1]+' '

#print()

#print('KAPPA', end='=')
lower=lower+'KAPPA='
for i in g1.atoms.indices:
    #print("150", end=',')
    lower=lower+"150,"
#print()

lower=lower[:-1]+' '

#print()

#print('EXP', end='=')
lower=lower+'EXP='
for i in g1.atoms.indices:
    #print("2", end=',')
    lower=lower+"2,"
#print()

lower=lower[:-1]+' '

#print()

#print('EPS', end='=')
lower=lower+'EPS='
for i in g1.atoms.indices:
    #print("1", end=',')
    lower=lower+"1,"
#print()

lower=lower[:-1]+' '

#print()

#print('OFFSET', end='=')
lower=lower+'OFFSET='
for i in g1.atoms.indices:
    #print("0", end=',')
    lower=lower+"0,"
#print()

lower=lower[:-1]+' '

#print()

#print('AT', end='=')
lower=lower+'AT='
for i in g1.atoms.indices:
    #print("2.4", end=',')
    lower=lower+"2.4,"
#print()

lower=lower[:-1]+' LABEL=lwall'

f.write(lower)
f.write('\n')
f.write('\n')
f.write('PRINT ARG=uwall.bias,lwall.bias STRIDE=500 FILE=cv.0.dat """\n')
f.write('    return script\n')
f.close()
