import numpy as np
import pandas as pd
try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
except ImportError:
    import simtk.openmm as mm
    import simtk.openmm.app as app
    import simtk.unit as unit
import mdtraj
import sys
import os

#__location__ = os.path.dirname(os.path.abspath(__file__))
#sys.path.append(f'{__location__}/../..')

from openabc.forcefields import MOFFMRGModel
from openabc.forcefields.parsers import MOFFParser, MRGdsDNAParser
from openabc.utils.shadow_map import load_ca_pairs_from_gmx_top

"""
Compare energy with GROMACS output. 
Note the native pairs in GROMACS topology file are produced by SMOG, which may be slightly different from the native pairs found by our shadow map algorithm code. 
To keep consistency, we directly load native pairs from GROMACS topology file. 
"""

salt = 200
print(f"Processing salt concentration: {salt}")
salt_string = str(salt) + 'mM'

histone_parser = MOFFParser.from_atomistic_pdb('histone.pdb', 'single_histone_CA.pdb', default_parse=False)
histone_parser.parse_mol(get_native_pairs=False)
histone_parser.parse_exclusions()

dsDNA_parser = MRGdsDNAParser.from_atomistic_pdb('dna.pdb', 'cg_single_nucl_dna.pdb')

protein_dna = MOFFMRGModel()
protein_dna.append_mol(histone_parser)
protein_dna.append_mol(dsDNA_parser)

## This parser belongs to the nucleosome that is the focus of the energy analysis
histone_parser_altered = MOFFParser.from_atomistic_pdb('histone.pdb', 'single_histone_CA.pdb', default_parse=False)
histone_parser_altered.parse_mol(get_native_pairs=False)
histone_parser_altered.parse_exclusions()

dsDNA_parser_altered = MRGdsDNAParser.from_atomistic_pdb('dna.pdb', 'cg_single_nucl_dna.pdb')

# Protein parser adjustments
histone_parser_altered.atoms['name'] = histone_parser_altered.atoms['name'].replace({'CA': 'NP'})

# DNA parser adjustments
dsDNA_parser_altered.atoms['name'] = dsDNA_parser_altered.atoms['name'].replace({'DN': 'ND'})


# we need to write all the atoms to pdb file
protein_dna.atoms_to_pdb('cg_protein_dna.pdb')

# Number of nucleosomes
n_nucl = 100
if not os.path.exists('cg_n_nucl.pdb'):
    insert_molecules('cg_protein_dna.pdb', 'cg_n_nucl.pdb', n_mol=n_nucl, box=[49, 49, 300],reset_serial=False)

## This will go through each of the 100 nucleosomes by index
openmm_energies = []
for i in range(0,100):
    protein_dna = MOFFMRGModel()
    if i == 0:
        protein_dna.append_mol(histone_parser_altered)
        protein_dna.append_mol(dsDNA_parser_altered)
        for j in range(1,n_nucl):
            protein_dna.append_mol(histone_parser)
            protein_dna.append_mol(dsDNA_parser)
    else:
        for j in range(0,i):
            protein_dna.append_mol(histone_parser)
            protein_dna.append_mol(dsDNA_parser)
        protein_dna.append_mol(histone_parser_altered)
        protein_dna.append_mol(dsDNA_parser_altered)
        for j in range(i+1,100):
            protein_dna.append_mol(histone_parser)
            protein_dna.append_mol(dsDNA_parser)

    top = app.PDBFile('cg_protein_dna.pdb').getTopology()
    protein_dna.create_system(top, box_a=49, box_b=49, box_c=300)


    manning_scale=1.0
    salt_conc = salt*unit.millimolar
    temperature = 300*unit.kelvin
    protein_dna.add_protein_bonds(force_group=1)
    protein_dna.add_protein_angles(force_group=2)
    protein_dna.add_dna_bonds(force_group=5)
    protein_dna.add_dna_angles(force_group=6)
    protein_dna.add_dna_fan_bonds(force_group=7)
    protein_dna.add_contacts_test(force_group=8)

    protein_dna.add_elec_switch_map_test(salt_conc, manning_scale, temperature, cutoff1=8.3*unit.nanometer, force_group=9)

    protein_dna.save_system('system.xml')
    collision = 1/unit.picosecond
    timestep = 10*unit.femtosecond
    integrator = mm.NoseHooverIntegrator(temperature, collision, timestep)
    platform_name = 'CUDA'
    protein_dna.set_simulation(integrator, platform_name, init_coord=None)
    simulation = protein_dna.simulation

    #The path of the trajectory for analysis
    traj_path = '100_nucleosome_condensate_200mM.xtc'

    traj = mdtraj.load_xtc(traj_path,top='cg_n_nucl.pdb')
    n_frames = traj.xyz.shape[0]
    nucleosome_energies = []
    for k in range(n_frames):
        row = []
        simulation.context.setPositions(traj.xyz[k])
        for j in range(9,10):
            state = protein_dna.simulation.context.getState(getEnergy=True, groups={j})
            energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            row.append(energy)
        nucleosome_energies.append(row)
    openmm_energies.append(nucleosome_energies)

openmm_energies = np.array(openmm_energies)

df_openmm_energies = pd.DataFrame(np.squeeze(openmm_energies)).round(6)
df_openmm_energies.round(2).to_csv('openmm_energies.csv', index=False)


