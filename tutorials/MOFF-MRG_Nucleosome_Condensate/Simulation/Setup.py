import numpy as np
import pandas as pd
import sys
import os
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit
import mdtraj


from openabc.forcefields.parsers import MOFFParser, MRGdsDNAParser
from openabc.forcefields import MOFFMRGModel
from openabc.utils.insert import insert_molecules
from openabc.forcefields.rigid import createRigidBodies
from nucl_utils import get_MOFF_MRG_chromatin_rigid_bodies

_amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS',
                'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                'LEU', 'LYS', 'MET', 'PHE', 'PRO',
                'SER', 'THR', 'TRP', 'TYR', 'VAL']

_nucleotides = ['DA', 'DT', 'DC', 'DG']


# set slab simulation box size
box_a = 49
box_b = 49
box_c = 300

# we have to rebuild the system as this time there is no MonteCarloBarostat in it
histone = MOFFParser.from_atomistic_pdb('starting_structure/histone.pdb', 'single_histone_CA.pdb', default_parse=False)
histone.parse_mol(get_native_pairs=False)
dna = MRGdsDNAParser.from_atomistic_pdb('starting_structure/dna.pdb', 'cg_single_nucl_dna.pdb', default_parse=False)
dna.parse_mol(bonded_energy_scale=0.8)

single_nucl = MOFFMRGModel()
single_nucl.append_mol(histone)
single_nucl.append_mol(dna)
single_nucl.atoms_to_pdb('cg_single_nucl.pdb')

n_atoms_per_single_nucl = len(single_nucl.atoms.index)

# Number of nucleosomes
n_nucl = 100

insert_molecules('cg_single_nucl.pdb', 'cg_n_nucl.pdb', n_mol=n_nucl, box=[box_a, box_b, box_c],reset_serial=False)

multi_nucl = MOFFMRGModel()
for i in range(n_nucl):
    multi_nucl.append_mol(histone)
    multi_nucl.append_mol(dna)

top = app.PDBFile('cg_n_nucl.pdb').getTopology()

# Create a system to start the NVT simulation
multi_nucl.create_system(top, box_a=box_a, box_b=box_b, box_c=box_c)

init_coord = app.PDBFile('./cg_n_nucl.pdb').getPositions()
# Create the rigid bodies;
single_nucl_rigid_body = np.array(get_MOFF_MRG_chromatin_rigid_bodies(n_nucl=1, nrl=147, n_rigid_bp_per_nucl=73)[0])
rigid_bodies = []
for i in range(n_nucl):
    rigid_bodies.append((single_nucl_rigid_body + i*n_atoms_per_single_nucl).tolist())
createRigidBodies(multi_nucl.system, init_coord, rigid_bodies)

# load trajectory and get the compressed configuration
# for easier visualization, we move the geometric center of all the atoms to the box center
npt_traj = mdtraj.load_dcd('../../../../equil-200mM/output_NPT.dcd', top='cg_n_nucl.pdb')
#npt_traj = mdtraj.load_dcd('./NPT_compress.dcd', top='cg_n_nucl.pdb')
init_coord = npt_traj.xyz[-1]
init_coord -= np.mean(init_coord, axis=0)
init_coord += 0.5*np.array([box_a, box_b, box_c])

# Start to run simulation;
salt_conc = 200*unit.millimolar
temperature = 300*unit.kelvin
manning_scale = 1.0
multi_nucl.add_protein_bonds(force_group=1)
multi_nucl.add_protein_angles(force_group=2)
#multi_nucl.add_protein_dihedrals(force_group=3)
#multi_nucl.add_native_pairs(force_group=4)
multi_nucl.add_dna_bonds(force_group=5)
multi_nucl.add_dna_angles(force_group=6)
multi_nucl.add_dna_fan_bonds(force_group=7)
multi_nucl.add_contacts(force_group=8)

multi_nucl.add_elec_switch_map(salt_conc, manning_scale, temperature, cutoff1=4.0*unit.nanometer, cutoff2=4.3*unit.nanometer, force_group=9)

# Save the system into system.xml
multi_nucl.save_system('system.xml')


# Begin the simulation
platform_name = 'CUDA'
#platform_name = 'CPU'

# Use LangevinMiddleIntegrator with friction coefficient=0.01/ps to speed it up; Langevin middle integrator is required for simulations with rigid bodies, see discussion with Peter Eastman: https://github.com/openmm/openmm/issues/3993;
friction_coeff = 0.01/unit.picosecond
timestep = 10*unit.femtosecond
integrator = mm.LangevinMiddleIntegrator(temperature, friction_coeff, timestep)
multi_nucl.set_simulation(integrator, platform_name, init_coord=init_coord)
simulation = multi_nucl.simulation

simulation.minimizeEnergy()
output_interval = 50000
output_dcd = './traj_temp300.dcd'
output_state = './state_report.txt'
dcd_reporter = app.DCDReporter(output_dcd, output_interval, enforcePeriodicBox=True)
state_reporter = app.StateDataReporter(output_state, output_interval, step=True, time=True, potentialEnergy=True,
                                       kineticEnergy=True, totalEnergy=True, temperature=True, speed=True)

simulation.reporters.append(dcd_reporter)
simulation.reporters.append(state_reporter)
simulation.context.setVelocitiesToTemperature(temperature)

# run NVT compression
print('Start NVT run.')
simulation.step(5000000)
print('NVT simulation is finished.')

# print final box vectors
state = simulation.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True,
                                    getParameters=True, enforcePeriodicBox=True)
box_vec = state.getPeriodicBoxVectors(asNumpy=True)

print('Final box vectors:')
print(box_vec)

# save the final state
with open('./state_cpt.xml', 'w') as f:
    f.write(mm.XmlSerializer.serialize(state))
~                                                        
