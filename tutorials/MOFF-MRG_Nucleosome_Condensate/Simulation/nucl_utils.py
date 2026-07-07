import numpy as np

# histone tails (atom index starts from 1): 1-43, 136-159, 238-257, 353-400, 488-530, 623-646, 725-744, 840-887
# in openmm, we use index starting from 0
_histone_tail_start_atoms = np.array([1, 136, 238, 353, 488, 623, 725, 840]) - 1
_histone_tail_end_atoms = np.array([43, 159, 257, 400, 530, 646, 744, 887]) - 1
_histone_tail_atoms = []
for i in range(8):
    _histone_tail_atoms += list(range(_histone_tail_start_atoms[i], _histone_tail_end_atoms[i] + 1))
_histone_tail_atoms = np.array(_histone_tail_atoms)

_n_bp_per_nucl = 147

_n_CA_atoms_per_histone = 974

_histone_core_atoms = np.array([x for x in range(_n_CA_atoms_per_histone) if x not in _histone_tail_atoms])

def get_MOFF_MRG_chromatin_rigid_bodies(n_nucl, nrl, n_rigid_bp_per_nucl=73):
    """
    Get chromatin rigid bodies for MOFF and MRG models. 
    The chromatin should possess uniform linker length without additional linkers on both ends. 
    
    Parameters
    ----------
    n_nucl : int
        Nucleosome number. 
    
    nrl : int
        Nucleosome repeat length. 
    
    n_flexible_bp_per_nucl : int
        The number of flexible nucleosomal base pairs for each nucleosome. 
    
   Returns
    -------
    rigid_bodies : list
        List of rigid bodies. 
    
    """
    n_bp = nrl*(n_nucl - 1) + _n_bp_per_nucl
    assert n_rigid_bp_per_nucl > 0
    n_CA_atoms = n_nucl*_n_CA_atoms_per_histone
    n_dna_atoms = 2*n_bp
    n_atoms = n_CA_atoms + n_dna_atoms
    bp_id_to_atom_id_dict = {}
    for i in range(n_bp):
        bp_id_to_atom_id_dict[i] = []
        # first ssDNA chain
        bp_id_to_atom_id_dict[i] += [n_CA_atoms + i]
        # second ssDNA chain
        bp_id_to_atom_id_dict[i] += [n_atoms - i - 1]
    rigid_bodies = []
    for i in range(n_nucl):
        rigid_bodies.append([])
        rigid_bodies[i] += (_histone_core_atoms + i*_n_CA_atoms_per_histone).tolist()
        start_bp_id = int((_n_bp_per_nucl - n_rigid_bp_per_nucl)/2) + i*nrl
        end_bp_id = start_bp_id + n_rigid_bp_per_nucl - 1
        for j in range(start_bp_id, end_bp_id + 1):
            rigid_bodies[i] += bp_id_to_atom_id_dict[j]
        rigid_bodies[i] = sorted(rigid_bodies[i])
    return rigid_bodies


                                                                                                          63,0-1        Bot
