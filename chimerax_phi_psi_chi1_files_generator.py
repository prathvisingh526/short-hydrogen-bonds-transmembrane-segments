import os
from chimerax.core.commands import run

path_molprobity_files = "D:\\files_molprobity\\tm\\" #chimeraX throws error if the path contains whitespaces
path_destination      = "D:\\phi_psi_chi1_tm\\"

for pdb_file in os.listdir(path_molprobity_files):
	run(session, "open " + path_molprobity_files + pdb_file)
	for dihedral_angle in ["phi", "psi", "chi1"]:
		run(session, "save " + path_destination + pdb_file[0:4].upper() + "_" + dihedral_angle + ".defattr attrName " + dihedral_angle + " format defattr model #1")
	run(session, "close")
run(session, "exit")