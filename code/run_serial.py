import numpy as np
import simtk.unit as u
from simtk.openmm import app
import simtk.openmm as mm
from repex.thermodynamics import ThermodynamicState
from repex.parallel_tempering import ParallelTempering
import repex.netcdf_io

name = "1vii"
equil_filename = "./equil/%s.pdb" % name

nc_filename = "./%s.nc" % name

n_replicas = 50 # number of temperature replicas
T_min = 300.0 * u.kelvin # minimum temperature
T_max = 450.0 * u.kelvin # maximum temperature

parameters = dict(pressure=1.0 * u.atmospheres, number_of_iterations=1000)

T_i = [ T_min + (T_max - T_min) * (np.exp(float(i) / float(n_replicas-1)) - 1.0) / (np.e - 1.0) for i in range(n_replicas)]

forcefield = app.ForceField("amber10.xml", "tip3p.xml")
pdb = app.PDBFile(equil_filename)

system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.CutoffPeriodic, nonbondedCutoff=1.0 * u.nanometers, constraints=app.HBonds)

coordinates = [pdb.positions] * n_replicas

rex = ParallelTempering.create(system, coordinates, nc_filename, T_min=T_min, T_max=T_max, n_temps=n_replicas, parameters=parameters)
rex.run()
