import re
import os
import pprint
import json
import math
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.PDB import *

def rescounter_tm(structure, chain_auth_asym_id, segments_relevant_chain, pdb_id_with_entity_id_cdhit, pdb_id_molprobity, dictdata, list_disulfides):
	chains    =  structure.get_chains()
	for chain in chains:
		parent_chain_id = chain.get_id()
		if parent_chain_id == chain_auth_asym_id and "HBCIL" in segments_relevant_chain.keys():
			for segment_range in segments_relevant_chain["HBCIL"]:
				lower_limit_of_segment_range = int(segment_range.split("-")[0])
				upper_limit_of_segment_range = int(segment_range.split("-")[1])
				residues                     = chain.get_residues()
				for residue in residues:
					if residue.get_resname() != "CYS":
						is_cys_disulfide_bonded = "NA"
					if residue.get_resname() == "CYS":
						is_cys_disulfide_bonded  = False
						reschain_plus_resno = str(residue.get_full_id()[2]) + "." + str(residue.get_full_id()[3][1]) # "reschain" plus "resno"
						if reschain_plus_resno in list_disulfides:
							is_cys_disulfide_bonded = True

					if residue.get_resname() in dictdata.keys() and residue.get_full_id()[3][1] in range(lower_limit_of_segment_range, upper_limit_of_segment_range + 1) and residue.get_full_id()[3][0] == " " and residue.get_full_id()[3][2] == " ": #if the relevant residue is non-heteroatom entity with no insertion code
						dictdata[residue.get_resname()].setdefault("Resname",              []).append(residue.get_resname())
						dictdata[residue.get_resname()].setdefault("Resno",                []).append(residue.get_full_id()[3][1])
						dictdata[residue.get_resname()].setdefault("Chain ID",             []).append(residue.get_full_id()[2])
						dictdata[residue.get_resname()].setdefault("Entity ID",            []).append(pdb_id_with_entity_id_cdhit)
						dictdata[residue.get_resname()].setdefault("PDB ID",               []).append(pdb_id_molprobity[0:4].upper())
						dictdata[residue.get_resname()].setdefault("CYS disulfide bonded", []).append(is_cys_disulfide_bonded)
