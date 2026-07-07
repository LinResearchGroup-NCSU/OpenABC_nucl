import MDAnalysis as mda
import numpy as np
import itertools

# ----------------------------
# USER INPUT: update these
# ----------------------------

import sys

salt = 
kappa = 

cutoff = float(kappa) * 10
#cutoff = 10.0
topology_file = "cg_n_nucl.pdb"
trajectory_file = "100_nucleosome_condensate_200mM.xtc"

# ----------------------------
# LOAD TRAJECTORY
# ----------------------------
u = mda.Universe(topology_file, trajectory_file)
n_frames = len(u.trajectory)
## ----------------------------
# Define ranges
# ----------------------------
histone_tail_segments = []
histone_tail_sets = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 623, 624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642, 643, 644, 645, 646, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 739, 740, 741, 742, 743, 744, 840, 841, 842, 843, 844, 845, 846, 847, 848, 849, 850, 851, 852, 853, 854, 855, 856, 857, 858, 859, 860, 861, 862, 863, 864, 865, 866, 867, 868, 869, 870, 871, 872, 873, 874, 875, 876, 877, 878, 879, 880, 881, 882, 883, 884, 885, 886, 887], dtype=int) - 1


DNA_ends_sets = []
original_DNA = list(itertools.chain(range(974,1004),range(1090,1150),range(1237,1267)))

total = 1268

for i in range(100):
    temp = [((i*total) + r) for r in histone_tail_sets]
    histone_tail_segments.append(temp)

    temp = [((i*total) + r) for r in original_DNA]
    DNA_ends_sets.append(temp)

## ----------------------------
# Compute coordinates
# ----------------------------
frame = 0
for ts in u.trajectory:
    network = np.zeros((100,100))
    ### rows are tails, columns are DNA
    for r in range(100):
        for c in range(100):
            if (r != c):
                groupA_indices = histone_tail_A_sets[r]
                groupA = u.atoms[groupA_indices]

                groupB_indices = DNA_ends_sets[c]
                groupB = u.atoms[groupB_indices]

                contacts = 0

                posA = groupA.positions
                posB = groupB.positions

                distances = np.linalg.norm(posA[:, None, :] - posB[None, :, :], axis=2)
                cn = np.sum(distances < cutoff)
                contacts = contacts + cn

                network[r,c] = contacts

    np.savetxt("contact_map_" + str(frame) + "_frame.csv", network, delimiter=",")
    frame = frame + 1


                                                                        
