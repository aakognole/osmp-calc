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

jobname = str(argv[1])+'_'+str(argv[2])+'_'+str(argv[3])

psf = app.CharmmPsfFile('../'+jobname+'.xplor.psf')
crd = app.CharmmCrdFile('./'+jobname+'.omm.crd')

psf.setBox(4.8*nanometer,4.8*nanometer,9.6*nanometer)

toppar = glob('../../toppar/toppar_drude/toppar_drude*str')
params = app.CharmmParameterSet(*toppar)

system = psf.createSystem(params,  nonbondedMethod=ff.PME,
                          nonbondedCutoff=1.2*nanometer,
                          switchDistance=1.0*nanometer,
                          ewaldErrorTolerance = 0.0001,
                          constraints=ff.HBonds)

#system.addForce(mm.MonteCarloBarostat(1*bar, 298*kelvin)) # This is NVT ensemble, no pressure control

script = plumedscript()
system.addForce(PlumedForce(script))

integrator = mm.DrudeLangevinIntegrator(303.15*kelvin, 5/picosecond, 1*kelvin, 20/picosecond,
                                        0.001*picoseconds)
integrator.setMaxDrudeDistance(0.02) # Drude Hardwall

simulation = app.Simulation(psf.topology, system, integrator)
simulation.context.setPositions(crd.positions)
simulation.context.computeVirtualSites()

start, end  = 1, 2
if start == 1:
    start = 1 #do a very short minimization    
    simulation.minimizeEnergy(maxIterations=20)
else:
    # start from the positions and velocities
    with open(jobname+'.'+str(start-1)+'.rst', 'r') as f:
        simulation.context.setState(mm.XmlSerializer.deserialize(f.read()))
    with open(jobname+'.'+str(start-1)+'.chk', 'rb') as f:
        simulation.context.loadCheckpoint(f.read())

nsavcrd = 500       # save frames every 0.5 ps
nstep   = 10000000  # simulate every 10 ns
nprint  = 10000     # report every 10 ps

for j in range(start,end):
    dcd = app.DCDReporter(jobname+'.'+str(j)+'.dcd', nsavcrd)
    firstdcdstep = (j-1)*nstep + nsavcrd
    dcd._dcd = app.DCDFile(dcd._out, simulation.topology, simulation.integrator.getStepSize(),
                           firstdcdstep, nsavcrd) # charmm doesn't like first step to be 0
    simulation.reporters.append(dcd)
    simulation.reporters.append(app.StateDataReporter(jobname+'.'+str(j)+'.out', nprint, step=True,
                                                      kineticEnergy=True, potentialEnergy=True,
                                                      totalEnergy=True, temperature=True, volume=True,
                                                      speed=True, remainingTime=True,
                                                      totalSteps=nstep, separator='\t')) 
    simulation.step(nstep)
    simulation.reporters.pop()
    simulation.reporters.pop()

    # write restart file
    state = simulation.context.getState( getPositions=True, getVelocities=True )
    with open(jobname+'.'+str(j)+'.rst', 'w') as f:
        f.write(mm.XmlSerializer.serialize(state))
    with open(jobname+'.'+str(j)+'.chk', 'wb') as f:
        f.write(simulation.context.createCheckpoint())


