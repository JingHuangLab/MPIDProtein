from openmm import app
from openmm.app import *
from openmm import *
from openmm.unit import *
import openmm.app.element as elem
import openmm.app.forcefield as forcefield
import openmm as mm
from sys import stdout, argv
import mpidplugin
import numpy as np
import argparse

jobname=sys.argv[1]
ffname=sys.argv[2]
start=int(sys.argv[3])
stop=int(sys.argv[4])

prePath='./'

pdb = app.PDBFile(jobname+'_'+ffname+'.pdb')
forcefield = ForceField('./para/mpid_2019g_220628_pme.xml')


system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, ewaldErrorTolerance=1e-4, nonbondedCutoff=12*angstrom, switchDistance=10*angstrom, polarization='extrapolated', constraints=HBonds, useDispersionCorrection=False )

for force in system.getForces():
     if type(force).__name__ == "NonbondedForce":
        print("nonbondForce",force.getCutoffDistance())
     if type(force).__name__ == "CustomNonbondedForce":
        print("customnonbondForce",force.getCutoffDistance())
        print(force.getSwitchingDistance())
        print(force.getUseSwitchingFunction())
        print("custom_lrc",force.getUseLongRangeCorrection())
     if type(force).__name__ == "Force":
        mpid_force = mpidplugin.MPIDForce.cast(force)
        print(mpid_force)
        print("MPID:",mpid_force.getCutoffDistance())

platform = mm.Platform.getPlatformByName('CUDA')
system.addForce(mm.MonteCarloBarostat(1*bar, 300*kelvin))
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtoseconds)
simulation = app.Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)

#context = simulation.context
if pdb.topology.getPeriodicBoxVectors():
    simulation.context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())

box=simulation.context.getState().getPeriodicBoxVectors()
print ('Box Size ')
print(box[0],"    ", box[1],"      ", box[2] ,"    ")

#print out energy and forces
state=simulation.context.getState(getPositions=True, getForces=True, getEnergy=True)
print ("energy:",state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))
print ("Force:")
for f in state.getForces().value_in_unit(kilocalories_per_mole/angstroms):
    print (f)
        
sys.stdout.flush()

if start > 0: # restart from the positions and velocities
    restart=str(start-1)
    #with open('dcd/'+jobname+'.'+restart+'.rst', 'r') as f:
    #    simulation.context.setState(mm.XmlSerializer.deserialize(f.read()))
    with open(prePath+'openmm_output/'+jobname+'.'+ffname+restart+'.chk', 'rb') as f:
        simulation.context.loadCheckpoint(f.read())

nsavcrd=10000 # save every 10 ps
nstep=10000000 # simulate every 10 ns
nprint=100000 # report every 100 ps


for ii in range(start,stop+1):
    dcd=app.DCDReporter(prePath+'openmm_output/'+jobname+'.'+ffname+str(ii)+'.dcd', nsavcrd)

    firstdcdstep = ii*nstep + nsavcrd
    while (firstdcdstep > 2000000000):
        firstdcdstep -= 2000000000
    dcd._dcd = app.DCDFile(dcd._out, simulation.topology, simulation.integrator.getStepSize(), firstdcdstep, nsavcrd) # charmm doesn't like first step to be 0
    simulation.reporters.append(dcd)
    simulation.reporters.append(app.StateDataReporter(prePath+'openmm_output/'+jobname+'.'+ffname+str(ii)+'.out', nprint, step=True, kineticEnergy=True, potentialEnergy=True, totalEnergy=True, temperature=True, volume=True, speed=True))
    simulation.step(nstep)
    simulation.reporters.pop()
    simulation.reporters.pop()
    dcd._out.close() # close the dcd file to make sure the buffers are written.

    # write restart file
    state = simulation.context.getState( getPositions=True, getVelocities=True )
    with open(prePath+'openmm_output/'+jobname+'.'+ffname+str(ii)+'.rst', 'w') as f:
        f.write(mm.XmlSerializer.serialize(state))
    with open(prePath+'openmm_output/'+jobname+'.'+ffname+str(ii)+'.chk', 'wb') as f:
        f.write(simulation.context.createCheckpoint())
