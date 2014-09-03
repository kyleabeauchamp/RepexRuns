import mdtraj as md
import os
import numpy as np
import simtk.unit as u
from simtk.openmm import app
import simtk.openmm as mm
import pdbfixer

name = "1vii"
pdb_filename = os.path.join("./pdbs/%s.pdb" % name)
equil_pdb_filename = os.path.join("./equil/%s.pdb" % name)
equil_dcd_filename = os.path.join("./equil/%s.dcd" % name)

temperature = 300 * u.kelvin
friction = 1.0 / u.picosecond
pressure = 1.0 * u.atmospheres
equil_timestep = 1.0 * u.femtosecond
barostat_frequency = 25
n_equil_steps = 100000
equil_output_frequency = 1000

forcefield = app.ForceField("amber10.xml", "tip3p.xml")

pdb = app.PDBFile(pdb_filename)
model = app.modeller.Modeller(pdb.topology, pdb.positions)
model.addSolvent(forcefield, padding=0.5 * u.nanometer)

topology = model.getTopology()
positions = model.getPositions()

system = forcefield.createSystem(model.topology, nonbondedMethod=app.CutoffPeriodic, nonbondedCutoff=1.0 * u.nanometers, constraints=app.HBonds)

integrator = mm.LangevinIntegrator(temperature, friction, equil_timestep)
system.addForce(mm.MonteCarloBarostat(pressure, temperature, barostat_frequency))

simulation = app.Simulation(topology, system, integrator)
simulation.context.setPositions(positions)

print('Minimizing.')
simulation.minimizeEnergy()

simulation.context.setVelocitiesToTemperature(temperature)
print('Equilibrating.')

simulation.reporters.append(app.PDBReporter(equil_pdb_filename, n_equil_steps - 1))
simulation.reporters.append(app.DCDReporter(equil_dcd_filename, equil_output_frequency))
simulation.step(n_equil_steps)
del simulation
del system
traj = md.load(equil_dcd_filename, top=equil_pdb_filename)[-1]
traj.save(equil_pdb_filename)
