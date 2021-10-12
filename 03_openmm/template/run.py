from simtk import openmm as mm
from simtk.openmm import app
from simtk.unit import *
from sys import stdout, exit, stderr, argv
from drudeutility import *

from simtk.openmm import *
from simtk.openmm.app import *

Platform.loadPluginsFromDirectory('/opt/openmm/7.4.0/plugins/plumed/plumed_in_openmm/platforms/cuda')
from openmmplumed import PlumedForce

from plumed import plumed

ion=str(argv[1])
mol=str(argv[2])
conc=str(argv[3])
start, end = 1, 2

jobname=ion+'_'+mol+'_'+conc

psf = CharmmDrudePsfFile('./'+jobname+'.xplor.psf')
crd = app.CharmmCrdFile('./'+jobname+'.omm.crd')

psf.setBox(4.8*nanometer,4.8*nanometer,9.6*nanometer)

params = app.CharmmParameterSet('./toppar_drude/toppar_drude_master_protein_2020b.str', './toppar_drude/toppar_drude_nucleic_acid_2020b.str', './toppar_drude/toppar_drude_model_2020b.str')

system = psf.createSystem(params,  nonbondedMethod=ff.PME, nonbondedCutoff=1.2*nanometer, switchDistance=1.0*nanometer, ewaldErrorTolerance = 0.0001, constraints=ff.HBonds)

nbforce = [system.getForce(i) for i in range(system.getNumForces()) if isinstance(system.getForce(i), mm.NonbondedForce)][0]
nbforce.setNonbondedMethod(mm.NonbondedForce.PME)
nbforce.setEwaldErrorTolerance(0.0001)
nbforce.setCutoffDistance(1.2*nanometer)
nbforce.setUseSwitchingFunction(True)
nbforce.setSwitchingDistance(1.0*nanometer)

# if there is nbfix
cstnb = [system.getForce(i) for i in range(system.getNumForces()) if isinstance(system.getForce(i), mm.CustomNonbondedForce)]
if cstnb:
    nbfix = cstnb[0]
    nbfix.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    nbfix.setCutoffDistance(1.2*nanometer)
    nbfix.setUseSwitchingFunction(True)
    nbfix.setSwitchingDistance(1.0*nanometer)

system.addForce(mm.MonteCarloBarostat(1*bar, 298*kelvin))

script = plumed()
system.addForce(PlumedForce(script))

integrator = mm.DrudeLangevinIntegrator(298*kelvin, 5/picosecond, 1*kelvin, 20/picosecond, 0.001*picoseconds)
integrator.setMaxDrudeDistance(0.02) # Drude Hardwall
simulation = app.Simulation(psf.topology, system, integrator)
simulation.context.setPositions(crd.positions)
simulation.context.computeVirtualSites()

# do a very short minimization
if start == 1:
    #do a very short minimization
    #simulation.minimizeEnergy(maxIterations=20)
    start = 1 
else:
    # start from the positions and velocities
    with open(jobname+'.'+str(start-1)+'.rst', 'r') as f:
        simulation.context.setState(mm.XmlSerializer.deserialize(f.read()))
    with open(jobname+'.'+str(start-1)+'.chk', 'rb') as f:
        simulation.context.loadCheckpoint(f.read())

nsavcrd=500 # save every 10 ps
nstep=10000000 # simulate every 10 ns
nprint=10000 # report every 10 ps

for ii in range(start,end):
    dcd=app.DCDReporter(jobname+'.'+str(ii)+'.dcd', nsavcrd)
    firstdcdstep = (ii-1)*nstep + nsavcrd
    dcd._dcd = app.DCDFile(dcd._out, simulation.topology, simulation.integrator.getStepSize(), firstdcdstep, nsavcrd) # charmm doesn't like first step to be 0
    simulation.reporters.append(dcd)
    simulation.reporters.append(app.StateDataReporter(jobname+'.'+str(ii)+'.out', nprint, step=True, kineticEnergy=True, potentialEnergy=True, totalEnergy=True, temperature=True, volume=True, speed=True, remainingTime=True, totalSteps=nstep, separator='\t')) 
    simulation.step(nstep)
    simulation.reporters.pop()
    simulation.reporters.pop()

    # write restart file
    state = simulation.context.getState( getPositions=True, getVelocities=True )
    with open(jobname+'.'+str(ii)+'.rst', 'w') as f:
        f.write(mm.XmlSerializer.serialize(state))
    with open(jobname+'.'+str(ii)+'.chk', 'wb') as f:
        f.write(simulation.context.createCheckpoint())


