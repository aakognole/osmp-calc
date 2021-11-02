from simtk import openmm as mm
from simtk.openmm import app
from simtk.openmm.app import (forcefield as ff)
from simtk.unit import *
from sys import stdout, exit, stderr, argv

from simtk.openmm import *
from simtk.openmm.app import *

#Platform.loadPluginsFromDirectory('/opt/openmm/7.4.0/plugins/plumed/plumed_in_openmm/platforms/cuda')
from openmmplumed import PlumedForce
from plumed import plumedscript

from glob import glob
import MDAnalysis as mda

jobname = str(argv[1])+'_'+str(argv[2])+'_'+str(argv[3])
run = str(argv[4])

psf = app.CharmmPsfFile('../'+jobname+'.xplor.psf')
crd = app.CharmmCrdFile('./'+jobname+'.omm.crd')

psf.setBox(4.8*nanometer,4.8*nanometer,9.6*nanometer)

toppar = glob('../../toppar/toppar_c36/*')
params = app.CharmmParameterSet(*toppar)

system = psf.createSystem(params,  nonbondedMethod=ff.PME,
                          nonbondedCutoff=1.2*nanometer,
                          switchDistance=1.0*nanometer,
                          ewaldErrorTolerance = 0.0001,
                          constraints=ff.HBonds)

#system.addForce(mm.MonteCarloBarostat(1*bar, 298*kelvin)) # This is NVT ensemble, no pressure control

script = plumedscript()
system.addForce(PlumedForce(script))

integrator = mm.DrudeLangevinIntegrator(298.15*kelvin, 5/picosecond, 0.001*picoseconds)

simulation = app.Simulation(psf.topology, system, integrator)
simulation.context.setPositions(crd.positions)

simulation.minimizeEnergy(maxIterations=20)

nsavcrd = 500       # save frames every 0.5 ps
nstep   = 10000000  # simulate every 10 ns
nprint  = 10000     # report every 10 ps

dcd = app.DCDReporter(jobname+'.'+run+'.dcd', nsavcrd)
firstdcdstep = (0)*nstep + nsavcrd
dcd._dcd = app.DCDFile(dcd._out, simulation.topology, simulation.integrator.getStepSize(),
                       firstdcdstep, nsavcrd) # charmm doesn't like first step to be 0
simulation.reporters.append(dcd)
simulation.reporters.append(app.StateDataReporter(jobname+'.'+run+'.out', nprint, step=True,
                                                  kineticEnergy=True, potentialEnergy=True,
                                                  totalEnergy=True, temperature=True, volume=True,
                                                  speed=True, remainingTime=True,
                                                  totalSteps=nstep, separator='\t')) 
simulation.step(nstep)
simulation.reporters.pop()
simulation.reporters.pop()

# write restart file
state = simulation.context.getState( getPositions=True, getVelocities=True )
with open(jobname+'.'+run+'.rst', 'w') as f:
    f.write(mm.XmlSerializer.serialize(state))
with open(jobname+'.'+run+'.chk', 'wb') as f:
    f.write(simulation.context.createCheckpoint())

# write trajectory stripping water
u = mda.Universe('../'+jobname+'.psf',jobname+'.'+run+'.dcd')
nowat = u.select_atoms('all and not resname SWM4 TIP3')
nowatdcd = jobname+'.'+run+'.nowat.dcd'
with MDAnalysis.Writer(nowatdcd, nowat.n_atoms) as W:
    for ts in u.trajectory:
        W.write(nowat)

exit()
