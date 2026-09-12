import pandas as pd
import numpy as np
import re
import os
import pprint
from Bio.PDB import *
import statistics
import sys
sys.path.append("G:\\My Drive\\python scripts and parameters\\")
import range_to_table
import importlib
importlib.reload(range_to_table)

def undesired_bifurcates_remover(df1, excel_filename, reskey, donor_atoms_ON, donor_atoms_S, acceptor_atoms_ON, acceptor_atoms_S, molprobity_files):
    df1.sort_values("Donor PDB ID", inplace = True)
    df2                       = df1.copy() #this line causes the code to throw the warning: "A value is trying to be set on a copy of a slice from a DataFrame". It doesn't mean that something is wrong with your code. It simply means that any changes made in "df2" will not affect "df1"
    dict_common_acceptor_info = dict()
    rows_indices_unwanted_bif = list()
    for i in df1.itertuples():
        acceptor_atom     = df1.at[i[0], "Acceptor atom"]
        acceptor_resname  = df1.at[i[0], "Acceptor resname"]
        acceptor_resno    = df1.at[i[0], "Acceptor resno"]
        acceptor_reschain = df1.at[i[0], "Acceptor reschain"]
        acceptor_pdb_id   = df1.at[i[0], "Acceptor PDB ID"]
        tup1              = tuple() + (acceptor_atom, acceptor_resname, acceptor_resno, acceptor_reschain, acceptor_pdb_id) #using tuple and not list because a list cannot become key of a dictionary
        dict_common_acceptor_info.setdefault(tup1, []).append(i[0])                                                         #it will look like: {('O', 'ALA', 52, 'A', '3DB2') [1299, 1300]} for an acceptor which is mentioned in >1 row

    for key, value in dict_common_acceptor_info.items():
        if len(value) > 1: #if len(value) > 1, it means that there are >1 rows with the same values of acceptor info.
            dict_donor = dict()
            for i in df1.itertuples():
                if i[0] in value:
                    donor_resno    = df1.at[i[0], "Donor resno"]
                    donor_reschain = df1.at[i[0], "Donor reschain"]
                    donor_pdb_id   = df1.at[i[0], "Donor PDB ID"]
                    ha_distance    = df1.at[i[0], "HA distance (Å)"]
                    tup2           = tuple() + (donor_resno, donor_reschain, donor_pdb_id)
                    dict_donor.setdefault(tup2, []).append(ha_distance) #if two donors atoms are donating to the same residue, they can either belong to the same residue or different residues. If they belong to same residue, then their HA distances will append to the same key

            for key2, value2 in dict_donor.items():
                if len(value2) > 1:                    #if "len(value2) > 1", it means that >1 donor atoms, belonging to the same residue, are donating to the same acceptor atom
                    for hadistance in value2:          #in such a case, only keep the instance which is closer to the acceptor atom
                        if hadistance != min(value2):  #if the value in the current iteration is not the smallest value (I am doing this because I want the index number of those rows whose HA distance isn't minimum so that I can delete them later)
                            for j in df2.itertuples(): #then search for it corresponding row in the dataframe & grab its row number
                                donor_resno    = df2.at[j[0], "Donor resno"]
                                donor_reschain = df2.at[j[0], "Donor reschain"]
                                donor_pdb_id   = df2.at[j[0], "Donor PDB ID"]
                                ha_distance    = df2.at[j[0], "HA distance (Å)"]
                                if donor_resno == key2[0] and donor_reschain == key2[1] and donor_pdb_id == key2[2] and ha_distance == hadistance:
                                    rows_indices_unwanted_bif.append(j[0])

    if len(rows_indices_unwanted_bif) == 0:
        print("No unwanted instances of bifurcated H-bonds found.")
        df3 = pd.DataFrame() #empty dataframe

    if len(rows_indices_unwanted_bif) > 0:
        print("Unwanted instances of bifurcated H-bonds found. Their info will be written in a separate tab of the excel file and then their info from the main dataframe will be deleted.")
        df3 = df2.loc[rows_indices_unwanted_bif]
        df3.sort_values("Donor PDB ID", inplace = True)

        for i in rows_indices_unwanted_bif[::-1]: #this loop will do the same thing as "for i in range(len(df1)-1, -1, -1)"; the loop being used right now requires one step less
            df1 = df1.drop([i], axis = 0)
        df1 = df1.reset_index(drop = True)

##################################################################################################################
#Create a datafame where DA occupancy is 1:
##################################################################################################################

    df4 = df1[(df1["Occupancy of donor atom"] == 1) & (df1["Occupancy of acceptor atom"] == 1)]
    df4 = df4.reset_index(drop = True)
    df4.sort_values("Donor PDB ID", inplace = True)

##################################################################################################################
#Create a datafame for self contacts:
##################################################################################################################

    df5 = df4[(df4["Donor resname"] == df4["Acceptor resname"]) & (df4["Donor resno"] == df4["Acceptor resno"]) & (df4["Donor reschain"] == df4["Acceptor reschain"])]
    df5 = df5.reset_index(drop = True) #if you don't reset index, then if a self contact was in the 7th row of df4 and it gets added to the 1st row of df5, its index will still stay 7
    df5.sort_values("Donor PDB ID", inplace = True)

#################################################################################################################################################################################################################
#Create a datafames which stores only those instances where the donor is mainchain atom (N-atom):
#################################################################################################################################################################################################################

    df6 = df4[df4["Donor atom"].isin(["N"])]
    if len(df6) != 0:
        df6.sort_values("Donor PDB ID", inplace = True)

#################################################################################################################################################################################################################
#Create a datafames which stores only those instances where the acceptor is mainchain atom (N-atom or O-atom):
#################################################################################################################################################################################################################
    
    df7 = df4[df4["Acceptor atom"].isin(["O", "N"])]
    if len(df7) != 0:
        df7.sort_values("Donor PDB ID", inplace = True)

#################################################################################################################################################################################################################
#Create a datafames which stores only those instances where the donor is sidechain atom:
#################################################################################################################################################################################################################

    df8 = df4[~df4["Donor atom"].isin(["N"])]
    if len(df8) != 0:
        df8.sort_values("Donor PDB ID", inplace = True)

#################################################################################################################################################################################################################
#Create a datafames which stores only those instances where the acceptor is sidechain atom:
#################################################################################################################################################################################################################

    df9 = df4[~df4["Acceptor atom"].isin(["O", "N"])]
    if len(df9) != 0:
        df9.sort_values("Donor PDB ID", inplace = True)

#################################################################################################################################################################################################################
#Create a dataframe which stores the total INSTANCES where the donor is mainchain atom (N-atom):
#################################################################################################################################################################################################################

    df_filtered       = df4[df4["Donor atom"].isin(["N"])]
    dict_value_counts = df_filtered["Donor atom"].value_counts().to_dict()

    if len(dict_value_counts) == 0:
        df10 = pd.DataFrame()
    if len(dict_value_counts) != 0:
        dict_instances_counts = dict()
        for key, value in dict_value_counts.items():
            dict_instances_counts.setdefault("Atom1", []).append("Total instances where mainchain " + key + "-atom is the donor")
            dict_instances_counts.setdefault("Atom2", []).append(key)
            dict_instances_counts.setdefault("Count", []).append(value)
        df10 = pd.DataFrame(dict_instances_counts)
        df10.sort_values("Count", inplace = True)

#################################################################################################################################################################################################################
#Create a dataframe which stores the total INSTANCES where the acceptor is mainchain atom (N-atom or mainchain O-atom):
#################################################################################################################################################################################################################

    df_filtered       = df4[df4["Acceptor atom"].isin(["O", "N"])] #this dataframe's length will be 0 if the objective is to find Hbond instances where MET side chain is acceptor
    dict_value_counts = df_filtered["Acceptor atom"].value_counts().to_dict()

    if len(dict_value_counts) == 0:
        df11 = pd.DataFrame()
    if len(dict_value_counts) != 0:
        dict_instances_counts = dict()
        for key, value in dict_value_counts.items():
            dict_instances_counts.setdefault("Atom1", []).append("Total instances where mainchain " + key + "-atom is the acceptor")
            dict_instances_counts.setdefault("Atom2", []).append(key)
            dict_instances_counts.setdefault("Count", []).append(value)
        df11 = pd.DataFrame(dict_instances_counts)
        df11.sort_values("Count", inplace = True)

#################################################################################################################################################################################################################
#Create a dataframe which stores the total INSTANCES where ASN side chain is the donor, the total INSTANCES where GLN side chain is the donor and so on:
#################################################################################################################################################################################################################

    df_filtered       = df4[~df4["Donor atom"].isin(["N"])]
    dict_value_counts = df_filtered["Donor resname"].value_counts().to_dict()

    if len(dict_value_counts) == 0:
        df12 = pd.DataFrame()
    if len(dict_value_counts) != 0:
        dict_instances_counts = dict()
        for key, value in dict_value_counts.items():
            dict_instances_counts.setdefault("Resname1", []).append("Total instances where " + key + " sidechain is the donor")
            dict_instances_counts.setdefault("Resname2", []).append(key)
            dict_instances_counts.setdefault("Count",    []).append(value)
        df12 = pd.DataFrame(dict_instances_counts)
        df12.sort_values("Count", inplace = True)

#################################################################################################################################################################################################################
#Create a dataframe which stores the total INSTANCES where ASN side chain is the acceptor, the total INSTANCES where GLN side chain is the acceptor and so on:
#################################################################################################################################################################################################################

    df_filtered       = df4[~df4["Acceptor atom"].isin(["O", "N"])] #when the acceptor IS NOT mainchain "O" or "N"
    dict_value_counts = df_filtered["Acceptor resname"].value_counts().to_dict()

    if len(dict_value_counts) == 0:
        df13 = pd.DataFrame()
    if len(dict_value_counts) != 0:
        dict_instances_counts = dict()
        for key, value in dict_value_counts.items():
            dict_instances_counts.setdefault("Resname1", []).append("Total instances where " + key + " sidechain is the acceptor")
            dict_instances_counts.setdefault("Resname2", []).append(key)
            dict_instances_counts.setdefault("Count",    []).append(value)
        df13 = pd.DataFrame(dict_instances_counts)
        df13.sort_values("Count", inplace = True)

#################################################################################################################################################################################################################
#Create dataframes which stores only unique donor residues, donating to SON, ON, O and N:
#################################################################################################################################################################################################################
    # For the "df14_SON", I have used the drop_duplicates() method to remove duplicates. But if a donor forms >1 interaction, then this method will cause retention of the first instance it encountered.
    # As a result, if a donor residue was donating to O as well as N and if the method encountered the "O" row first, it will retain the "O" row's donor info only.
    # As a result, it will give the wrong impression that the unique donor is donating only to "O" atom. This is why I will retain only the donor info, else in future, an unsuspecting me can draw wrong conclusions from it.
    # I am also not using "H-atom involved" because if two rows corresponding to the same residue differ by "H-atom involved", the code will identify them as two different residues"

    df14_SON = df4.drop_duplicates(subset = ["Donor resname", "Donor resno", "Donor reschain", "Donor PDB ID", "Donor insertion code"])
    df14_SON = df14_SON[["Donor resname", "Donor resno", "Donor reschain", "Donor PDB ID", "Donor insertion code", "Donor sec. str.", "Chi1 angle donor (degrees)", "Phi angle donor (degrees)", "Psi angle donor (degrees)", "SASA donor (Å)^2"]]
    df14_SON = df14_SON.reset_index(drop = True)
    df14_SON.sort_values("Donor PDB ID", inplace = True)

    list_indices_ON = list()
    list_indices_O  = list()
    list_indices_N  = list()

    for i in df14_SON.itertuples():
        donor_resname  = df14_SON.at[i[0], "Donor resname"]
        donor_resno    = df14_SON.at[i[0], "Donor resno"]
        donor_reschain = df14_SON.at[i[0], "Donor reschain"]
        donor_pdb_id   = df14_SON.at[i[0], "Donor PDB ID"]
        list_indices_temp = list()
        for j in df4.itertuples():
            donor_resname_copy  = df4.at[j[0], "Donor resname"]
            donor_resno_copy    = df4.at[j[0], "Donor resno"]
            donor_reschain_copy = df4.at[j[0], "Donor reschain"]
            donor_pdb_id_copy   = df4.at[j[0], "Donor PDB ID"]
            if donor_resname == donor_resname_copy and donor_resno == donor_resno_copy and donor_reschain == donor_reschain_copy and donor_pdb_id == donor_pdb_id_copy:
                acc_atom_type = df4.at[j[0], "Acceptor atom"][0] # grab only the 1st alphabet of the acceptor atom
                list_indices_temp.append(acc_atom_type) # suppose a given donor residue donates to "O" as well as "N", the list will contain Os and Ns

        if "O" in list_indices_temp:
            list_indices_O.append(i[0])
        if "N" in list_indices_temp:
            list_indices_N.append(i[0])
        if "O" in list_indices_temp or "N" in list_indices_temp:
            list_indices_ON.append(i[0])

    df14_O  = df14_SON.loc[list_indices_O] # it will contain the info of all those unique donor residues which donated to O-atom (they may or may not have donated to S or N)
    df14_N  = df14_SON.loc[list_indices_N] # similar logic as above as above
    df14_ON = df14_SON.loc[list_indices_ON]

    df14_O  = df14_O.reset_index(drop = True)
    df14_N  = df14_N.reset_index(drop = True)
    df14_ON = df14_ON.reset_index(drop = True)

#################################################################################################################################################################################################################
#Create a dataframe which stores only unique acceptor residues:
#################################################################################################################################################################################################################

    df17 = df4.drop_duplicates(subset = ["Acceptor resname", "Acceptor resno", "Acceptor reschain", "Acceptor PDB ID", "Acceptor insertion code"])
    df17 = df17[["Acceptor resname", "Acceptor resno", "Acceptor reschain", "Acceptor PDB ID", "Acceptor insertion code", "Acceptor sec. str.", "Chi1 angle acceptor (degrees)", "Phi angle acceptor (degrees)", "Psi angle acceptor (degrees)", "SASA acceptor (Å)^2"]]
    df17.sort_values("Acceptor PDB ID", inplace = True) #If an acceptor forms >1 interaction, this dataframe will contain the very 1st instance among all

#################################################################################################################################################################################################################
#Create a datafame which stores the total number of acceptors per UNIQUE donor RESIDUE:
#################################################################################################################################################################################################################

    dict_acceptor_count_for_excel = dict()
    divider_symbol                = "."
    for i in df14_SON.itertuples():
        acceptor_info_all = list()
        unique_don_info = str(df14_SON.at[i[0], "Donor PDB ID"]) + divider_symbol + df14_SON.at[i[0], "Donor resname"] + divider_symbol + str(df14_SON.at[i[0], "Donor resno"]) + divider_symbol + df14_SON.at[i[0], "Donor reschain"]
        for j in df4.itertuples():
            don_info = str(df4.at[j[0], "Donor PDB ID"]) + divider_symbol + df4.at[j[0], "Donor resname"]  + divider_symbol + str(df4.at[j[0], "Donor resno"])  + divider_symbol + df4.at[j[0], "Donor reschain"]
            if unique_don_info == don_info:
                common_accep_info = str(df4.at[j[0], "Acceptor PDB ID"]) + divider_symbol + df4.at[j[0], "Acceptor atom"] + divider_symbol + df4.at[j[0], "Acceptor resname"] + divider_symbol + str(df4.at[j[0], "Acceptor resno"]) + divider_symbol +  df4.at[j[0], "Acceptor reschain"]
                acceptor_info_all.append(common_accep_info)

        acceptor_atoms_O    = [acceptor_atom for acceptor_atom in acceptor_atoms_ON if "O" in acceptor_atom]
        acceptor_atoms_N    = [acceptor_atom for acceptor_atom in acceptor_atoms_ON if "N" in acceptor_atom]

        ctr_mainchain_n     = 0
        ctr_mainchain_o     = 0
        ctr_mainchain_total = 0

        ctr_sidechain_n     = 0
        ctr_sidechain_o     = 0
        ctr_sidechain_s     = 0
        ctr_sidechain_total = 0

        ctr_n_total         = 0
        ctr_o_total         = 0

        amino_acids         = {'ALA': False, 'ARG': False, 'ASN': False, 'ASP': False, 'CYS': False, 'GLN': False, 'GLU': False, 'GLY': False, 'HIS': False, 'ILE': False, 'LEU': False, 'LYS': False, 'MET': False, 'PHE': False, 'PRO': False, 'SER': False, 'THR': False, 'TYR': False, 'TRP': False,  'VAL': False}
    
        for value in acceptor_info_all:
            match1 = re.search(r"(\w+)(\.)(\w+)(\.)(\w+)(\.)(\d+)(\.)(\w+)", value) #this regex is different from the one used for "donor_info_all" as acceptor info doesn't include the "H-atom involved" info
            if match1:
                if match1.group(3) in acceptor_atoms_S:
                    ctr_sidechain_s += 1
                    if match1.group(5) in amino_acids.keys():
                        amino_acids[match1.group(5)] = True

                if match1.group(3) in acceptor_atoms_O:
                    if match1.group(3) == "O":
                        ctr_mainchain_o += 1
                    if match1.group(3) != "O":
                        ctr_sidechain_o += 1
                        if match1.group(5) in amino_acids.keys():
                            amino_acids[match1.group(5)] = True

                if match1.group(3) in acceptor_atoms_N:
                    if match1.group(3) == "N":
                        ctr_mainchain_n += 1
                    if match1.group(3) != "N":
                        ctr_sidechain_n += 1
                        if match1.group(5) in amino_acids.keys():
                            amino_acids[match1.group(5)] = True

        ctr_sidechain_total = ctr_sidechain_n + ctr_sidechain_o + ctr_sidechain_s
        ctr_mainchain_total = ctr_mainchain_n + ctr_mainchain_o
        ctr_n_total         = ctr_sidechain_n + ctr_mainchain_n
        ctr_o_total         = ctr_sidechain_o + ctr_mainchain_o

        dict_acceptor_count_for_excel.setdefault("Donor info",                                   []).append(unique_don_info)
        dict_acceptor_count_for_excel.setdefault("Acceptor info",                                []).append(",".join(acceptor_info_all))
        dict_acceptor_count_for_excel.setdefault("Acceptor count total",                         []).append(len(acceptor_info_all))
        dict_acceptor_count_for_excel.setdefault("Acceptor count sidechain S",                   []).append(ctr_sidechain_s)
        dict_acceptor_count_for_excel.setdefault("Acceptor count sidechain O",                   []).append(ctr_sidechain_o)
        dict_acceptor_count_for_excel.setdefault("Acceptor count sidechain N",                   []).append(ctr_sidechain_n)
        dict_acceptor_count_for_excel.setdefault("Acceptor count sidechain total",               []).append(ctr_sidechain_total)
        dict_acceptor_count_for_excel.setdefault("Acceptor count mainchain O",                   []).append(ctr_mainchain_o)
        dict_acceptor_count_for_excel.setdefault("Acceptor count mainchain N",                   []).append(ctr_mainchain_n)
        dict_acceptor_count_for_excel.setdefault("Acceptor count mainchain total",               []).append(ctr_mainchain_total)
        dict_acceptor_count_for_excel.setdefault("Acceptor count mainchain N plus sidechain N",  []).append(ctr_n_total)
        dict_acceptor_count_for_excel.setdefault("Acceptor count mainchain O plus sidechain O",  []).append(ctr_o_total)
        for key, value in amino_acids.items():
            dict_acceptor_count_for_excel.setdefault(key + " sidechain only",                    []).append(value)

    df15 = pd.DataFrame(dict_acceptor_count_for_excel)
    df15.sort_values("Donor info", inplace = True)

#################################################################################################################################################################################################################
#Create a datafame which stores the total number of donors per UNIQUE acceptor RESIDUE:
#################################################################################################################################################################################################################

    dict_donor_count_for_excel = dict()
    divider_symbol             = "."
    for i in df17.itertuples():
        donor_info_all    = list()
        unique_accep_info = str(df17.at[i[0], "Acceptor PDB ID"]) + divider_symbol + df17.at[i[0], "Acceptor resname"] + divider_symbol + str(df17.at[i[0], "Acceptor resno"]) + divider_symbol + df17.at[i[0], "Acceptor reschain"]

        for j in df4.itertuples():
            accep_info = str(df4.at[j[0], "Acceptor PDB ID"]) +     divider_symbol + df4.at[j[0], "Acceptor resname"] + divider_symbol + str(df4.at[j[0], "Acceptor resno"]) + divider_symbol + df4.at[j[0], "Acceptor reschain"]

            if unique_accep_info == accep_info:
                common_donor_info = str(df4.at[j[0], "Donor PDB ID"]) + divider_symbol + df4.at[j[0], "Donor atom"] + divider_symbol + df4.at[j[0], "H-atom involved"] + divider_symbol + df4.at[j[0], "Donor resname"] + divider_symbol + str(df4.at[j[0], "Donor resno"]) + divider_symbol +  df4.at[j[0], "Donor reschain"]
                donor_info_all.append(common_donor_info)

        donor_atoms_O       = [donor_atom for donor_atom in donor_atoms_ON if "O" in donor_atom]
        donor_atoms_N       = [donor_atom for donor_atom in donor_atoms_ON if "N" in donor_atom]

        ctr_n_mainchain     = 0 #there is no "ctr_mainchain_o" variable because the mainchain O-atom cannot donate a H-bond

        ctr_o_sidechain     = 0
        ctr_n_sidechain     = 0
        ctr_s_sidechain     = 0
        ctr_sidechain_total = 0

        ctr_n_total         = 0

        amino_acids         = {'ALA': False, 'ARG': False, 'ASN': False, 'ASP': False, 'CYS': False, 'GLN': False, 'GLU': False, 'GLY': False, 'HIS': False, 'ILE': False, 'LEU': False, 'LYS': False, 'MET': False, 'PHE': False, 'PRO': False, 'SER': False, 'THR': False, 'TYR': False, 'TRP': False,  'VAL': False}
        
        for value in donor_info_all:
            match1 = re.search(r"(\w+)(\.)(\w+)(\.)(\w+)(\.)(\w+)(\.)(\d+)(\.)(\w+)", value)
            if match1:
                if match1.group(3) in donor_atoms_S:
                    ctr_s_sidechain += 1
                    if match1.group(7) in amino_acids.keys():
                        amino_acids[match1.group(7)] = True

                if match1.group(3) in donor_atoms_O:
                    ctr_o_sidechain += 1
                    if match1.group(7) in amino_acids.keys():
                        amino_acids[match1.group(7)] = True

                if match1.group(3) in donor_atoms_N:
                    if match1.group(3) == "N":
                        ctr_n_mainchain += 1
                    if match1.group(3) != "N":
                        ctr_n_sidechain += 1
                        if match1.group(7) in amino_acids.keys():
                            amino_acids[match1.group(7)] = True

        ctr_sidechain_total = ctr_o_sidechain + ctr_n_sidechain + ctr_s_sidechain
        ctr_n_total         = ctr_n_sidechain + ctr_n_mainchain

        dict_donor_count_for_excel.setdefault("Acceptor info",                            []).append(unique_accep_info)
        dict_donor_count_for_excel.setdefault("Donor info",                               []).append(",".join(donor_info_all))
        dict_donor_count_for_excel.setdefault("Donor count total",                        []).append(len(donor_info_all))
        dict_donor_count_for_excel.setdefault("Donor count sidechain S",                  []).append(ctr_s_sidechain)
        dict_donor_count_for_excel.setdefault("Donor count sidechain O",                  []).append(ctr_o_sidechain)
        dict_donor_count_for_excel.setdefault("Donor count sidechain N",                  []).append(ctr_n_sidechain)
        dict_donor_count_for_excel.setdefault("Donor count sidechain total",              []).append(ctr_sidechain_total)
        dict_donor_count_for_excel.setdefault("Donor count mainchain N",                  []).append(ctr_n_mainchain)
        dict_donor_count_for_excel.setdefault("Donor count mainchain N plus sidechain N", []).append(ctr_n_total)
        for key, value in amino_acids.items():
            dict_donor_count_for_excel.setdefault(key + " sidechain only",                []).append(value)

    df18 = pd.DataFrame(dict_donor_count_for_excel)
    df18.sort_values("Acceptor info", inplace = True)

##############################################################################################################################################################
# Create dataframes containing normalized B-factors of unique donor residues which are donating to ON, which are donating to O and which are donating to N:
##############################################################################################################################################################
    def normalized_bfactor_saver(df_unique_residue, molprobity_files):
        pdbparser = PDBParser(PERMISSIVE = 1, QUIET = True)
        dictzzz   = dict()
        list_aas  = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

        for i in df_unique_residue.itertuples():
            bfact_resname  = df_unique_residue.at[i[0], "Donor resname"]
            bfact_resno    = df_unique_residue.at[i[0], "Donor resno"]
            bfact_reschain = df_unique_residue.at[i[0], "Donor reschain"]
            bfact_pdb_id   = df_unique_residue.at[i[0], "Donor PDB ID"]

            for molprob_file in os.listdir(molprobity_files):
                if molprob_file[0:4].upper() == bfact_pdb_id:
                    bfact_structure = pdbparser.get_structure("protein", molprobity_files + molprob_file)
                    ##################################################################################################################################################################
                    # first find the average and standard deviation of B-factor of all the CA atoms in donor residue's chain (and not from all the chains; I confirmed it with sir):
                    ##################################################################################################################################################################
                    bfact_CA        = list()
                    bfact_residues  = bfact_structure.get_residues()

                    for bfact_residue in bfact_residues:
                        if (bfact_residue.get_full_id()[2] == bfact_reschain and # if the residue belongs to the same chain as the unique donor
                            bfact_residue.get_full_id()[3][0] == " " and         # if the residue is a non-heteroatom entity
                            bfact_residue.get_full_id()[3][2] == " " and         # if the residue has no insertion code
                            bfact_residue.get_resname() in list_aas):            # if the residue is one of the 20 aacids
                            bfact_atoms = bfact_residue.get_atoms()

                            for bfact_atom in bfact_atoms:
                                if bfact_atom.get_name() == "CA":
                                    bfact_CA.append(bfact_atom.get_bfactor())

                    average_bfact_CA_all      = round((sum(bfact_CA) / len(bfact_CA)), 2)
                    # standard_dev_bfact_CA_all = round(np.std(bfact_CA), 2)
                    standard_dev_bfact_CA_all = round(statistics.stdev(bfact_CA), 2) # "statistics.stdev()" calculates SD using the sample formula, unlike "np.std()" which calculates SD using the population formula
                    ########################################################################################################################
                    # then find the normalized B-factor (using 2 ways) of the unique donor residue in the iteration:
                    ########################################################################################################################
                    bfact_residues = bfact_structure.get_residues() #maybe it is also a generator thats why I initialized it again
                    for bfact_residue in bfact_residues:
                        if bfact_residue.get_resname() == bfact_resname and bfact_residue.get_full_id()[3][1] == bfact_resno and bfact_residue.get_full_id()[2] == bfact_reschain:
                            bfact_atoms       = bfact_residue.get_atoms()
                            bfact_heavy_atoms = list()

                            for bfact_atom in bfact_atoms:
                                if bfact_atom.get_name()[0] != "H":
                                    bfact_heavy_atoms.append(bfact_atom.get_bfactor())

                            avg_b_factor_heavy_atoms = round((sum(bfact_heavy_atoms) / len(bfact_heavy_atoms)), 2)
                            sum_b_factor_heavy_atoms = sum(bfact_heavy_atoms)
                            norm_b_factor1           = round(((avg_b_factor_heavy_atoms - average_bfact_CA_all) / standard_dev_bfact_CA_all), 2)
                            norm_b_factor2           = round(((sum_b_factor_heavy_atoms - average_bfact_CA_all) / standard_dev_bfact_CA_all), 2) 

                            dictzzz.setdefault("Donor resname",                                                                 []).append(bfact_resname)
                            dictzzz.setdefault("Donor resno",                                                                   []).append(bfact_resno)
                            dictzzz.setdefault("Donor reschain",                                                                []).append(bfact_reschain)
                            dictzzz.setdefault("Donor PDB ID",                                                                  []).append(bfact_pdb_id)
                            dictzzz.setdefault("Average B-factor of all heavy atoms of donor residue",                          []).append(avg_b_factor_heavy_atoms)
                            dictzzz.setdefault("Sum B-factor of all heavy atoms of donor residue",                              []).append(sum_b_factor_heavy_atoms)
                            dictzzz.setdefault("Average B-factor of all CA atoms in donor residue's chain",                     []).append(average_bfact_CA_all)
                            dictzzz.setdefault("Stdev B-factor of all CA atoms in donor residue's chain",                       []).append(standard_dev_bfact_CA_all)
                            dictzzz.setdefault("Normalized B-factor based on average B-factor of heavy atoms of donor residue", []).append(norm_b_factor1)
                            dictzzz.setdefault("Normalized B-factor based on sum of B-factor of heavy atoms of donor residue",  []).append(norm_b_factor2)
        unique_res_df = pd.DataFrame(dictzzz)
        if len(unique_res_df) != 0:
            unique_res_df.sort_values("Donor PDB ID", inplace = True)
        return unique_res_df

    df_B_factors_of_unique_don_res_accep_is_O = normalized_bfactor_saver(df14_O, molprobity_files)
    df_B_factors_of_unique_don_res_accep_is_N = normalized_bfactor_saver(df14_N, molprobity_files)
    df_B_factors_of_unique_don_res_accep_is_ON = normalized_bfactor_saver(df14_ON, molprobity_files)

#############################################################################################################################################
# Create dataframes containing normalized B-factors of donor residues engaged in SELF CONTACTS, donating to mainchain O and mainchain N:
#############################################################################################################################################

    unique_unique_self_cont_don_res_accep_is_mainchain_O = df5[(df5["Acceptor atom"] == "O")]
    unique_unique_self_cont_don_res_accep_is_mainchain_O = unique_unique_self_cont_don_res_accep_is_mainchain_O.reset_index(drop = True)
    df_unique_self_cont_don_res_accep_is_mainchain_O     = normalized_bfactor_saver(unique_unique_self_cont_don_res_accep_is_mainchain_O, molprobity_files)

    unique_unique_self_cont_don_res_accep_is_mainchain_N = df5[(df5["Acceptor atom"] == "N")]
    unique_unique_self_cont_don_res_accep_is_mainchain_N = unique_unique_self_cont_don_res_accep_is_mainchain_N.reset_index(drop = True)
    df_unique_self_cont_don_res_accep_is_mainchain_N     = normalized_bfactor_saver(unique_unique_self_cont_don_res_accep_is_mainchain_N, molprobity_files)

    unique_unique_self_cont_don_res_accep_is_mainchain_ON = df5
    df_unique_self_cont_don_res_accep_is_mainchain_ON     = normalized_bfactor_saver(unique_unique_self_cont_don_res_accep_is_mainchain_ON, molprobity_files)

# #####################################################################################################################################################################################################
# # Create a dataframe which stores, out of all the UNIQUE donor residues, how many are donating to SER sidechain, how many to THR sidechain and so on. I will not be counting for mainchain:
# #####################################################################################################################################################################################################
#     #The code takes a lot of time to execute (for CYS/SER/THR globular dataset, not for MET globular) if I don't use "if" statments for the below codes
#     if "MET is the acceptor" in excel_filename:
#         df16 = pd.DataFrame()
#     if "MET is the acceptor" not in excel_filename:
#         distinct_resnames = list()
#         for i in df15.itertuples():
#             info = df15.at[i[0], "Donor info"]
#             match1 = re.search(r"(\w+)(\.)(\w+)(\.)(\d+)(\.)(\w+)", info)
#             if match1:
#                 distinct_resnames.append(match1.group(3))
#         distinct_resnames = set(distinct_resnames)
#         distinct_resnames = list(distinct_resnames)

#         count = dict()
#         for i in distinct_resnames:
#             df_filtered = df15[df15["Donor info"].str.contains(i)]
#             for colname in df_filtered.columns:
#                 if colname[0:3] in amino_acids.keys():
#                     count.setdefault("Residue",  []).append(i + " donating to " + colname[0:3] + " sidechain group")
#                     count.setdefault("Residue2", []).append(colname[0:3])
#                     count.setdefault("Count",    []).append(df_filtered[colname].value_counts().get(True, 0)) #if no instances of True is found, it will append 0
#         df16 = pd.DataFrame(count)

# ######################################################################################################################################################################################################
# # Create a dataframe which stores, out of all the UNIQUE acceptor residues, how many are accepting from SER sidechain, how many from THR sidechain and so on. I will not be counting for mainchain:
# ######################################################################################################################################################################################################
#     #The code takes a lot of time to execute (for CYS/SER/THR and MET globular dataset) if I don't use "if" statments for the below codes
#     if "MET is the acceptor" not in excel_filename:
#         df19 =pd.DataFrame()
#     if "MET is the acceptor" in excel_filename:
#         distinct_resnames = list()
#         for i in df18.itertuples():
#             info = df18.at[i[0], "Acceptor info"]
#             match1 = re.search(r"(\w+)(\.)(\w+)(\.)(\d+)(\.)(\w+)", info)
#             if match1:
#                 distinct_resnames.append(match1.group(3))
#         distinct_resnames = set(distinct_resnames)
#         distinct_resnames = list(distinct_resnames)

#         count = dict()
#         for i in distinct_resnames:
#             df_filtered = df18[df18["Acceptor info"].str.contains(i)]
#             for colname in df_filtered.columns:
#                 if colname[0:3] in amino_acids.keys():
#                     count.setdefault("Residue",  []).append(i + " accepting from " + colname[0:3] + " sidechain group")
#                     count.setdefault("Residue2", []).append(colname[0:3])
#                     count.setdefault("Count",    []).append(df_filtered[colname].value_counts().get(True, 0)) #if no instances of True is found, it will append 0
#         df19 = pd.DataFrame(count)

#######################################################################################################################################################
#Create a datafame which stores DA distance in primary sequence for cases where DA occupancy is 1:
#######################################################################################################################################################

    def my_func(i):
        if i < 0:
            return "i" + str(i) # the "-" symbol is already part of "i"
        if i > 0:
            return "i+" + str(i)
        if i == 0:
            return "i"

    dict_counts = dict()
    for i in range(-9,9+1):

        if i == -9:
            dict_counts.setdefault("DA distance in primary sequence",              []).append("i<-8")
            dict_counts.setdefault("Total instances",                              []).append(df4[(df4["DA distance in primary sequence"] < -8)].shape[0]) # shape[0] counts the total number of rows
            if "sidechain is the acceptor" in excel_filename:
                dict_counts.setdefault("Mainchain donor instances",                []).append(df4[(df4["DA distance in primary sequence"] < -8) & (df4["Donor atom"].isin(["N"]))].shape[0])
                dict_counts.setdefault("Sidechain donor instances",                []).append(df4[(df4["DA distance in primary sequence"] < -8) & (~df4["Donor atom"].isin(["N"]))].shape[0])
                for x in ["O", "N", "S"]:
                    dict_counts.setdefault("Instances when donor atom is " + x,    []).append(df4[(df4["DA distance in primary sequence"] < -8) & (df4["Donor atom"].str[0] == x)].shape[0])
            if "sidechain is the donor" in excel_filename:
                dict_counts.setdefault("Mainchain acceptor instances",             []).append(df4[(df4["DA distance in primary sequence"] < -8) & (df4["Acceptor atom"].isin(["N", "O"]))].shape[0])
                dict_counts.setdefault("Sidechain acceptor instances",             []).append(df4[(df4["DA distance in primary sequence"] < -8) & (~df4["Acceptor atom"].isin(["N", "O"]))].shape[0])
                for x in ["O", "N", "S"]:
                    dict_counts.setdefault("Instances when acceptor atom is " + x, []).append(df4[(df4["DA distance in primary sequence"] < -8) & (df4["Acceptor atom"].str[0] == x)].shape[0])

        if i == 9:
            dict_counts.setdefault("DA distance in primary sequence",              []).append("i>+8")
            dict_counts.setdefault("Total instances",                              []).append(df4[(df4["DA distance in primary sequence"] > 8)].shape[0]) # shape[0] counts the total number of rows
            if "sidechain is the acceptor" in excel_filename:
                dict_counts.setdefault("Mainchain donor instances",                []).append(df4[(df4["DA distance in primary sequence"] > 8) & (df4["Donor atom"].isin(["N"]))].shape[0])
                dict_counts.setdefault("Sidechain donor instances",                []).append(df4[(df4["DA distance in primary sequence"] > 8) & (~df4["Donor atom"].isin(["N"]))].shape[0])
                for x in ["O", "N", "S"]:
                    dict_counts.setdefault("Instances when donor atom is " + x,    []).append(df4[(df4["DA distance in primary sequence"] > 8) & (df4["Donor atom"].str[0] == x)].shape[0])
            if "sidechain is the donor" in excel_filename:
                dict_counts.setdefault("Mainchain acceptor instances",             []).append(df4[(df4["DA distance in primary sequence"] > 8) & (df4["Acceptor atom"].isin(["N", "O"]))].shape[0])
                dict_counts.setdefault("Sidechain acceptor instances",             []).append(df4[(df4["DA distance in primary sequence"] > 8) & (~df4["Acceptor atom"].isin(["N", "O"]))].shape[0])
                for x in ["O", "N", "S"]:
                    dict_counts.setdefault("Instances when acceptor atom is " + x, []).append(df4[(df4["DA distance in primary sequence"] > 8) & (df4["Acceptor atom"].str[0] == x)].shape[0])

        if i not in [-9,9]:
            dict_counts.setdefault("DA distance in primary sequence",              []).append(my_func(i))
            dict_counts.setdefault("Total instances",                              []).append(df4[(df4["DA distance in primary sequence"] == i)].shape[0])
            if "sidechain is the acceptor" in excel_filename:
                dict_counts.setdefault("Mainchain donor instances",                []).append(df4[(df4["DA distance in primary sequence"] == i) & (df4["Donor atom"].isin(["N"]))].shape[0])
                dict_counts.setdefault("Sidechain donor instances",                []).append(df4[(df4["DA distance in primary sequence"] == i) & (~df4["Donor atom"].isin(["N"]))].shape[0])
                for x in ["O", "N", "S"]:
                    dict_counts.setdefault("Instances when donor atom is " + x,    []).append(df4[(df4["DA distance in primary sequence"] == i) & (df4["Donor atom"].str[0] == x)].shape[0])
            if "sidechain is the donor" in excel_filename:
                dict_counts.setdefault("Mainchain acceptor instances",             []).append(df4[(df4["DA distance in primary sequence"] == i) & (df4["Acceptor atom"].isin(["N", "O"]))].shape[0])
                dict_counts.setdefault("Sidechain acceptor instances",             []).append(df4[(df4["DA distance in primary sequence"] == i) & (~df4["Acceptor atom"].isin(["N", "O"]))].shape[0])
                for x in ["O", "N", "S"]:
                    dict_counts.setdefault("Instances when acceptor atom is " + x, []).append(df4[(df4["DA distance in primary sequence"] == i) & (df4["Acceptor atom"].str[0] == x)].shape[0])

    df20 = pd.DataFrame(dict_counts)

# ##################################################################################################################################################
# #Create a dataframes that contains all the information that I need to put in powerpoint for case when CYS SER THR sidechain groups act as donor:
# ##################################################################################################################################################

#     dict_ppt_sidechain_is_donor = dict()
#     if "sidechain is the donor" in excel_filename:

#         if "TM H-bonds " in excel_filename:
#             df_rescounts = pd.read_excel("D:\\My Data\\Study Stuff\\PhD Work\\Part 3\\files_excel\\Residuewise counts in TM segments of unique chains in tm dataset.xlsx", sheet_name = reskey)
#         if "Globular H-bonds " in excel_filename:
#             df_rescounts = pd.read_excel("D:\\My Data\\Study Stuff\\PhD Work\\Project 02 Selenium H-bonds & chalcogen bonds\\20210302\\files_excel\\Latest\\Globular residuewise info of all residues in unique chains.xlsx", sheet_name = reskey)

#         if reskey == "CYS":
#             sidechain_group = "SH"
#         if reskey in ["SER", "THR"]:
#             sidechain_group = "OH"

#         dict_ppt_sidechain_is_donor.setdefault("Col1", [
#             "Total count of " + reskey, 
#             "Total count of " + reskey + " engaged as donor in " + sidechain_group + "…S/O/N H-bonding", 
#             "Total count of " + reskey + " engaged as donor in " + sidechain_group + "…O H-bonding", 
#             "Total count of " + reskey + " engaged as donor in " + sidechain_group + "…N H-bonding", 
#             "Total count of " + reskey + " engaged as donor in " + sidechain_group + "…O mainchain H-bonding", 
#             "Total count of " + reskey + " engaged as donor in " + sidechain_group + "…N mainchain H-bonding", 
#             "Total count of " + reskey + " engaged as donor in " + sidechain_group + "…O sidechain H-bonding", 
#             "Total count of " + reskey + " engaged as donor in " + sidechain_group + "…N sidechain H-bonding", 
#             "Total count of " + reskey + " engaged as donor in " + sidechain_group + "…S sidechain H-bonding",
#             "Total count of " + reskey + " engaged as donor in " + sidechain_group + "…O self-contact H-bonding", 
#             "Total count of " + reskey + " engaged as donor in " + sidechain_group + "…N self-contact H-bonding"
#             ])

#         dict_ppt_sidechain_is_donor.setdefault("Count", [
#             len(df_rescounts), 
#             len(df15), 
#             len(df15[df15["Acceptor count mainchain O plus sidechain O"] > 0 ]), 
#             len(df15[df15["Acceptor count mainchain N plus sidechain N"] > 0 ]), 
#             len(df15[df15["Acceptor count mainchain O"] > 0 ]), 
#             len(df15[df15["Acceptor count mainchain N"] > 0 ]), 
#             len(df15[df15["Acceptor count sidechain O"] > 0 ]), 
#             len(df15[df15["Acceptor count sidechain N"] > 0 ]), 
#             len(df15[df15["Acceptor count sidechain S"] > 0 ]), 
#             len(df5[df5["Acceptor atom"].isin(["O"])]), 
#             len(df5[df5["Acceptor atom"].isin(["N"])])
#             ])

#         dict_ppt_sidechain_is_donor.setdefault("Percent of total " + reskey, [
#             round((len(df_rescounts) / len(df_rescounts)) * 100, 2), 
#             round((len(df15) / len(df_rescounts)) * 100, 2), 
#             round((len(df15[df15["Acceptor count mainchain O plus sidechain O"] > 0 ]) / len(df_rescounts)) * 100, 2), 
#             round((len(df15[df15["Acceptor count mainchain N plus sidechain N"] > 0 ]) / len(df_rescounts)) * 100, 2), 
#             round((len(df15[df15["Acceptor count mainchain O"] > 0 ]) / len(df_rescounts)) * 100, 2), 
#             round((len(df15[df15["Acceptor count mainchain N"] > 0 ]) / len(df_rescounts)) * 100, 2), 
#             round((len(df15[df15["Acceptor count sidechain O"] > 0 ]) / len(df_rescounts)) * 100, 2), 
#             round((len(df15[df15["Acceptor count sidechain N"] > 0 ]) / len(df_rescounts)) * 100, 2), 
#             round((len(df15[df15["Acceptor count sidechain S"] > 0 ]) / len(df_rescounts)) * 100, 2), 
#             round((len(df5[df5["Acceptor atom"].isin(["O"])]) / len(df_rescounts)) * 100, 2), 
#             round((len(df5[df5["Acceptor atom"].isin(["N"])]) / len(df_rescounts)) * 100, 2)
#             ])

#         dict_ppt_sidechain_is_donor.setdefault("Percent of total " + reskey + " engaged as donor", [
#             round((len(df_rescounts) / len(df15)) * 100, 2), 
#             round((len(df15) / len(df15)) * 100, 2), 
#             round((len(df15[df15["Acceptor count mainchain O plus sidechain O"] > 0 ]) / len(df15)) * 100, 2), 
#             round((len(df15[df15["Acceptor count mainchain N plus sidechain N"] > 0 ]) / len(df15)) * 100, 2), 
#             round((len(df15[df15["Acceptor count mainchain O"] > 0 ]) / len(df15)) * 100, 2), 
#             round((len(df15[df15["Acceptor count mainchain N"] > 0 ]) / len(df15)) * 100, 2), 
#             round((len(df15[df15["Acceptor count sidechain O"] > 0 ]) / len(df15)) * 100, 2), 
#             round((len(df15[df15["Acceptor count sidechain N"] > 0 ]) / len(df15)) * 100, 2), 
#             round((len(df15[df15["Acceptor count sidechain S"] > 0 ]) / len(df15)) * 100, 2), 
#             round((len(df5[df5["Acceptor atom"].isin(["O"])]) / len(df15)) * 100, 2), 
#             round((len(df5[df5["Acceptor atom"].isin(["N"])]) / len(df15)) * 100, 2)
#             ])

#     df21 = pd.DataFrame(dict_ppt_sidechain_is_donor)

# #########################################################################################################################################################
# #Create a dataframes that contains all the information that I need to put in powerpoint for case when CYS SER THR MET sidechain groups act as acceptor:
# #########################################################################################################################################################

#     dict_ppt_sidechain_is_acceptor = dict()
#     if "sidechain is the acceptor" in excel_filename:

#         if "TM H-bonds " in excel_filename:
#             df_rescounts = pd.read_excel("D:\\My Data\\Study Stuff\\PhD Work\\Part 3\\files_excel\\Residuewise counts in TM segments of unique chains in tm dataset.xlsx", sheet_name = reskey)
#         if "Globular H-bonds " in excel_filename:
#             df_rescounts = pd.read_excel("D:\\My Data\\Study Stuff\\PhD Work\\Project 02 Selenium H-bonds & chalcogen bonds\\20210302\\files_excel\\20230826\\Globular residuewise info of all residues in unique chains.xlsx", sheet_name = reskey)

#         if reskey in ["CYS", "MET"]:
#             sidechain_group = "S"
#         if reskey in ["SER", "THR"]:
#             sidechain_group = "O"

#         dict_ppt_sidechain_is_acceptor.setdefault("Col1", [
#             "Total count of " + reskey, 
#             "Total count of " + reskey + " engaged as acceptor in SH/OH/NH..." + sidechain_group + " H-bonding", 
#             "Total count of " + reskey + " engaged as acceptor in OH..."       + sidechain_group + " H-bonding", 
#             "Total count of " + reskey + " engaged as acceptor in NH..."       + sidechain_group + " H-bonding", 
#             "Total count of " + reskey + " engaged as acceptor in NH..."       + sidechain_group + " mainchain H-bonding", 
#             "Total count of " + reskey + " engaged as acceptor in SH..."       + sidechain_group + " sidechain H-bonding", 
#             "Total count of " + reskey + " engaged as acceptor in OH..."       + sidechain_group + " sidechain H-bonding", 
#             "Total count of " + reskey + " engaged as acceptor in NH..."       + sidechain_group + " sidechain H-bonding",
#             "Total count of " + reskey + " engaged as acceptor in NH..."       + sidechain_group + " self-contact H-bonding"
#             ])

#         dict_ppt_sidechain_is_acceptor.setdefault("Count", [
#             len(df_rescounts), 
#             len(df18), 
#             len(df18[df18["Donor count sidechain O"] > 0 ]), 
#             len(df18[df18["Donor count mainchain N plus sidechain N"] > 0 ]), 
#             len(df18[df18["Donor count mainchain N"] > 0 ]), 
#             len(df18[df18["Donor count sidechain S"] > 0 ]), #
#             len(df18[df18["Donor count sidechain O"] > 0 ]), #
#             len(df18[df18["Donor count sidechain N"] > 0 ]), 
#             len(df5[df5["Donor atom"].isin(["N"])])
#             ])

#         dict_ppt_sidechain_is_acceptor.setdefault("Percent of total " + reskey, [
#             round((len(df_rescounts) / len(df_rescounts)) * 100, 2), 
#             round((len(df18) / len(df_rescounts)) * 100, 2), 
#             round((len(df18[df18["Donor count sidechain O"] > 0 ]) / len(df_rescounts)) * 100, 2), 
#             round((len(df18[df18["Donor count mainchain N plus sidechain N"] > 0 ]) / len(df_rescounts)) * 100, 2), 
#             round((len(df18[df18["Donor count mainchain N"] > 0 ]) / len(df_rescounts)) * 100, 2), 
#             round((len(df18[df18["Donor count sidechain S"] > 0 ]) / len(df_rescounts)) * 100, 2), 
#             round((len(df18[df18["Donor count sidechain O"] > 0 ]) / len(df_rescounts)) * 100, 2), 
#             round((len(df18[df18["Donor count sidechain N"] > 0 ]) / len(df_rescounts)) * 100, 2),  
#             round((len(df5[df5["Acceptor atom"].isin(["N"])]) / len(df_rescounts)) * 100, 2)
#             ])

#         dict_ppt_sidechain_is_acceptor.setdefault("Percent of total " + reskey + " engaged as acceptor", [
#             round((len(df_rescounts) / len(df18)) * 100, 2), 
#             round((len(df18) / len(df18)) * 100, 2), 
#             round((len(df18[df18["Donor count sidechain O"] > 0 ]) / len(df18)) * 100, 2), 
#             round((len(df18[df18["Donor count mainchain N plus sidechain N"] > 0 ]) / len(df18)) * 100, 2), 
#             round((len(df18[df18["Donor count mainchain N"] > 0 ]) / len(df18)) * 100, 2), 
#             round((len(df18[df18["Donor count sidechain S"] > 0 ]) / len(df18)) * 100, 2), 
#             round((len(df18[df18["Donor count sidechain O"] > 0 ]) / len(df18)) * 100, 2), 
#             round((len(df18[df18["Donor count sidechain N"] > 0 ]) / len(df18)) * 100, 2),  
#             round((len(df5[df5["Acceptor atom"].isin(["N"])]) / len(df18)) * 100, 2)
#             ])

#     df22 = pd.DataFrame(dict_ppt_sidechain_is_acceptor)

#######################################################################################################################################################################################################################################################
# Create a dataframe that contains a summary of secondary structure info where acceptor atom (when dealing with Hbonds sidechain is donor file) or donor atom (when dealing with Hbonds sidechain is acceptor file) starts either with O or N for cases where DA occupancy is 1:
#######################################################################################################################################################################################################################################################

    if "sidechain is the donor" in excel_filename:
        filter_word1 = "Acceptor atom"
    if "sidechain is the acceptor" in excel_filename:
        filter_word1 = "Donor atom"

    listz = list()
    for atomz in ["O", "N"]:
        dict_ssec_summary = dict()
        df_ssec1          = df4[df4["Donor sec. str."]              != "Stride file not found"]   # this will also cover the condition where the column "Acceptor sec. str." equals "Stride file not found"
        df_ssec2          = df_ssec1[df_ssec1["Donor sec. str."]    != "Residue entry not found"]
        df_ssec3          = df_ssec2[df_ssec2["Acceptor sec. str."] != "Residue entry not found"]
        df_ssec4          = df_ssec3[df_ssec3[filter_word1].str.startswith(atomz)][["Donor sec. str.", "Acceptor sec. str."]] # the condition when donor atom is "O" will never be true when dealing with Hbonds sidechain is acceptor file because mainchain O cannot act as donor but still I am using it cuz it will be used for donor files

        ssec_don_unique   = sorted(df_ssec4["Donor sec. str."].unique().tolist())    # grab all the unique sec. str. values from donor and acceptor cols
        ssec_accep_unique = sorted(df_ssec4["Acceptor sec. str."].unique().tolist())

        for val_don in ssec_don_unique:
            dict_ssec_summary.setdefault("Col1",                         []).append("Donor in " + val_don)
            for val_accep in ssec_accep_unique:
                dict_ssec_summary.setdefault("Acceptor in " + val_accep, []).append(df_ssec3[(df_ssec3["Donor sec. str."] == val_don) & (df_ssec3["Acceptor sec. str."] == val_accep)].shape[0]) # ".shape[0]" will count the rows where te given condition is satisfied

        listz.append(pd.DataFrame(dict_ssec_summary))

    df_ssec_all_O = listz[0]
    df_ssec_all_N = listz[1]


##################################################################################################################################################
# Create a dataframe that contains a summary of secondary structure info of donor and acceptor residues for cases where DA occupancy is 1:
##################################################################################################################################################

    dict_ssec_summary = dict()
    df_ssec1          = df4[df4["Donor sec. str."]              != "Stride file not found"]   # this will also cover the condition where the column "Acceptor sec. str." equals "Stride file not found"
    df_ssec2          = df_ssec1[df_ssec1["Donor sec. str."]    != "Residue entry not found"]
    df_ssec3          = df_ssec2[df_ssec2["Acceptor sec. str."] != "Residue entry not found"]

    ssec_don_unique   = sorted(df_ssec3["Donor sec. str."].unique().tolist())    # grab all the unique sec. str. values from donor and acceptor cols
    ssec_accep_unique = sorted(df_ssec3["Acceptor sec. str."].unique().tolist())

    for val_don in ssec_don_unique:
        dict_ssec_summary.setdefault("Col1",                         []).append("Donor in " + val_don)
        for val_accep in ssec_accep_unique:
            dict_ssec_summary.setdefault("Acceptor in " + val_accep, []).append(df_ssec3[(df_ssec3["Donor sec. str."] == val_don) & (df_ssec3["Acceptor sec. str."] == val_accep)].shape[0]) # ".shape[0]" will count the rows where te given condition is satisfied

    df_ssec_all_ON = pd.DataFrame(dict_ssec_summary)


##################################################################################################################
# Create dataframes for charting:
##################################################################################################################

    if "sidechain is the donor" in excel_filename:
        filter_word1 = "Acceptor atom"
        filter_word2 = "donor"
    if "sidechain is the acceptor" in excel_filename:
        filter_word1 = "Donor atom"
        filter_word2 = "acceptor"
    ################################################################################################################
    # For DA occup 1 charting:
    ################################################################################################################
    df_da1 = df4[["DA distance (Å)"]]
    df_da2 = df4[df4[filter_word1].str.startswith("O")][["DA distance (Å)"]]
    df_da3 = df4[df4[filter_word1].str.startswith("N")][["DA distance (Å)"]]

    df_ha1 = df4[["HA distance (Å)"]]
    df_ha2 = df4[df4[filter_word1].str.startswith("O")][["HA distance (Å)"]]
    df_ha3 = df4[df4[filter_word1].str.startswith("N")][["HA distance (Å)"]]

    df_dha1 = df4[["DHA angle (degrees)"]]
    df_dha2 = df4[df4[filter_word1].str.startswith("O")][["DHA angle (degrees)"]]
    df_dha3 = df4[df4[filter_word1].str.startswith("N")][["DHA angle (degrees)"]]

    df_haan1 = df4[["HAAn angle/HABisector angle (degrees)"]]
    df_haan2 = df4[df4[filter_word1].str.startswith("O")][["HAAn angle/HABisector angle (degrees)"]]
    df_haan3 = df4[df4[filter_word1].str.startswith("N")][["HAAn angle/HABisector angle (degrees)"]]

    df_chi11 = df14_ON[["Chi1 angle " + filter_word2 + " (degrees)"]]
    df_chi12 = df14_O[["Chi1 angle " + filter_word2 + " (degrees)"]]
    df_chi13 = df14_N[["Chi1 angle " + filter_word2 + " (degrees)"]]

    df_prep1    = df14_ON[df14_ON["Phi angle "   + filter_word2 + " (degrees)"] != "None"] # first remove those rows where the phi angle is None
    df_prep2    = df_prep1[df_prep1["Psi angle " + filter_word2 + " (degrees)"] != "None"] # then remove those rows where the psi angle is None
    df_phi_psi1 = df_prep2[["Phi angle " + filter_word2 + " (degrees)", "Psi angle " + filter_word2 + " (degrees)"]]

    df_prep1    = df14_O[df14_O["Phi angle "     + filter_word2 + " (degrees)"] != "None"]
    df_prep2    = df_prep1[df_prep1["Psi angle " + filter_word2 + " (degrees)"] != "None"]
    df_phi_psi2 = df_prep2[["Phi angle " + filter_word2 + " (degrees)", "Psi angle " + filter_word2 + " (degrees)"]]

    df_prep1    = df14_N[df14_N["Phi angle "     + filter_word2 + " (degrees)"] != "None"]
    df_prep2    = df_prep1[df_prep1["Psi angle " + filter_word2 + " (degrees)"] != "None"]
    df_phi_psi3 = df_prep2[["Phi angle " + filter_word2 + " (degrees)", "Psi angle " + filter_word2 + " (degrees)"]]

    df_sasa1 = df14_ON[["SASA " + filter_word2 + " (Å)^2"]]
    df_sasa2 = df14_O[["SASA " + filter_word2 + " (Å)^2"]]
    df_sasa3 = df14_N[["SASA " + filter_word2 + " (Å)^2"]]

    ################################################################################################################
    # For self contacts charting:
    ################################################################################################################
    df_da4 = df5[["DA distance (Å)"]]
    df_da5 = df5[df5[filter_word1].str.startswith("O")][["DA distance (Å)"]]
    df_da6 = df5[df5[filter_word1].str.startswith("N")][["DA distance (Å)"]]

    df_ha4 = df5[["HA distance (Å)"]]
    df_ha5 = df5[df5[filter_word1].str.startswith("O")][["HA distance (Å)"]]
    df_ha6 = df5[df5[filter_word1].str.startswith("N")][["HA distance (Å)"]]

    df_dha4 = df5[["DHA angle (degrees)"]]
    df_dha5 = df5[df5[filter_word1].str.startswith("O")][["DHA angle (degrees)"]]
    df_dha6 = df5[df5[filter_word1].str.startswith("N")][["DHA angle (degrees)"]]

    df_haan4 = df5[["HAAn angle/HABisector angle (degrees)"]]
    df_haan5 = df5[df5[filter_word1].str.startswith("O")][["HAAn angle/HABisector angle (degrees)"]]
    df_haan6 = df5[df5[filter_word1].str.startswith("N")][["HAAn angle/HABisector angle (degrees)"]]

    df_chi14 = df5[["Chi1 angle " + filter_word2 + " (degrees)"]]
    df_chi15 = df5[df5[filter_word1].str.startswith("O")][["Chi1 angle " + filter_word2 + " (degrees)"]]
    df_chi16 = df5[df5[filter_word1].str.startswith("N")][["Chi1 angle " + filter_word2 + " (degrees)"]]

    df_prep4 = df5[df5["Phi angle " + filter_word2 + " (degrees)"] != "None"]
    df_prep5 = df_prep4[df_prep4["Psi angle " + filter_word2 + " (degrees)"] != "None"]

    df_phi_psi4 = df_prep5[["Phi angle " + filter_word2 + " (degrees)", "Psi angle " + filter_word2 + " (degrees)"]]
    df_phi_psi5 = df_prep5[df_prep5[filter_word1].str.startswith("O")][["Phi angle " + filter_word2 + " (degrees)", "Psi angle " + filter_word2 + " (degrees)"]]
    df_phi_psi6 = df_prep5[df_prep5[filter_word1].str.startswith("N")][["Phi angle " + filter_word2 + " (degrees)", "Psi angle " + filter_word2 + " (degrees)"]]

    df_sasa4 = df5[["SASA " + filter_word2 + " (Å)^2"]]
    df_sasa5 = df5[df5[filter_word1].str.startswith("O")][["SASA " + filter_word2 + " (Å)^2"]]
    df_sasa6 = df5[df5[filter_word1].str.startswith("N")][["SASA " + filter_word2 + " (Å)^2"]]

    df_ssec1 = df5[df5["Donor sec. str."]              != "Stride file not found"]   # this will also cover the condition where the column "Acceptor sec. str." equals "Stride file not found"
    df_ssec2 = df_ssec1[df_ssec1["Donor sec. str."]    != "Residue entry not found"]
    df_ssec3 = df_ssec2[df_ssec2["Acceptor sec. str."] != "Residue entry not found"]

    df_ssec_self_contacts_ON         = df_ssec3[filter_word2.title() + " sec. str."].value_counts().reset_index() # dont use "drop = True" else the dataframe's column headers will go away; I am using two sets of brackets around the colname was resulting in an incorrect dataframe
    df_ssec_self_contacts_ON.columns = [filter_word2.title() + " sec. str.", "Count"] # this dataframe will have 2 columns named "donor/acceptor sec. str." and "Count"

    df_ssec_self_contacts_O          = df_ssec3[df_ssec3[filter_word1].str.startswith("O")][filter_word2.title() + " sec. str."].value_counts().reset_index()
    df_ssec_self_contacts_O.columns  = [filter_word2.title() + " sec. str.", "Count"]

    df_ssec_self_contacts_N          = df_ssec3[df_ssec3[filter_word1].str.startswith("N")][filter_word2.title() + " sec. str."].value_counts().reset_index()
    df_ssec_self_contacts_N.columns  = [filter_word2.title() + " sec. str.", "Count"]

##################################################################################################################
#Write dataframes to excel file:
##################################################################################################################

    if "sidechain is the donor" in excel_filename:
        filter_word2 = " don. "
    if "sidechain is the acceptor" in excel_filename:
        filter_word2 = " accep. "

    list_of_worksheets = [
    [df1, "All inst. undesired bif excl."], 
    [df2,  "All inst. undesired bif incl."], 
    [df3,  "Unwanted bif rows"], 
    [df4,  "Instances DA occup one"], 
    [df5,  "Self contact instances"], 
    [df6,  "Mainchain donor instances"], 
    [df7,  "Mainchain acceptor instances"], 
    [df8,  "Sidechain donor instances"], 
    [df9,  "Sidechain acceptor instances"], 
    [df10, "Mainchain don inst cnts"], 
    [df11, "Mainchain accep inst cnts"],
    [df12, "Sidechain don inst cnts"], 
    [df13, "Sidechain accep inst cnts"], 
    [df14_SON, "Unique donor residues"], 
    [df14_ON, "ON_donor residues"], 
    [df14_O, "O_donor residues"], 
    [df14_N, "N_donor residues"], 
    [df15, "Accep per unique donor"], 
    [df_B_factors_of_unique_don_res_accep_is_ON, "ON_B factor donor res"], 
    [df_B_factors_of_unique_don_res_accep_is_O,  "O_B factor donor res"], 
    [df_B_factors_of_unique_don_res_accep_is_N,  "N_B factor donor res"], 
    [df_unique_self_cont_don_res_accep_is_mainchain_ON, "ON_B factor self cont. res"], 
    [df_unique_self_cont_don_res_accep_is_mainchain_O, "O_B factor self cont. res"], 
    [df_unique_self_cont_don_res_accep_is_mainchain_N, "N_B factor self cont. res"], 
    # [df17, "Unique acceptor resiudes"], 
    # [df18, "Donors per unique accep"], 
    # [df16, "AAAA"], 
    # [df19, "BBBB"], 
    # [df21, "PPT data when " + reskey + " is don"], 
    # [df22, "PPT data when " + reskey + " is accep"], 
    [df20, "Positional pref."], 
    [df_ssec_all_ON, "ON_ssec. all"], 
    [df_ssec_all_O, "O_ssec. all"], 
    [df_ssec_all_N, "N_ssec. all"], 
    [df_ssec_self_contacts_ON, "ON_ssec. self cont. res"], 
    [df_ssec_self_contacts_O, "O_ssec. self cont. res"], 
    [df_ssec_self_contacts_N, "N_ssec. self cont. res"], 
    [df_da1, "ON_DA distance all"], 
    [df_da2, "O_DA distance all"], 
    [df_da3, "N_DA distance all"], 
    [df_da4, "ON_DA distance self cont."], 
    [df_da5, "O_DA distance self cont."], 
    [df_da6, "N_DA distance self cont."], 
    [df_ha1, "ON_HA distance all"], 
    [df_ha2, "O_HA distance all"], 
    [df_ha3, "N_HA distance all"], 
    [df_ha4, "ON_HA distance self cont."], 
    [df_ha5, "O_HA distance self cont."], 
    [df_ha6, "N_HA distance self cont."], 
    [df_dha1, "ON_DHA angle all"], 
    [df_dha2, "O_DHA angle all"], 
    [df_dha3, "N_DHA angle all"], 
    [df_dha4, "ON_DHA angle self cont."], 
    [df_dha5, "O_DHA angle self cont."], 
    [df_dha6, "N_DHA angle self cont."], 
    [df_haan1, "ON_HAAn angle all"], 
    [df_haan2, "O_HAAn angle all"], 
    [df_haan3, "N_HAAn angle all"], 
    [df_haan4, "ON_HAAn angle self cont."], 
    [df_haan5, "O_HAAn angle self cont."], 
    [df_haan6, "N_HAAn angle self cont."], 
    [df_chi11, "ON_Chi1 angle" + filter_word2 + "res"], 
    [df_chi12, "O_Chi1 angle" + filter_word2 + "res"], 
    [df_chi13, "N_Chi1 angle" + filter_word2 + "res"], 
    [df_chi14, "ON_Chi1 angle self cont. res"], 
    [df_chi15, "O_Chi1 angle self cont. res"], 
    [df_chi16, "N_Chi1 angle self cont. res"], 
    [df_phi_psi1, "ON_Phi psi" + filter_word2 + "res"], 
    [df_phi_psi2, "O_Phi psi" + filter_word2 + "res"], 
    [df_phi_psi3, "N_Phi psi" + filter_word2 + "res"], 
    [df_phi_psi4, "ON_Phi psi self cont. res"], 
    [df_phi_psi5, "O_Phi psi self cont. res"], 
    [df_phi_psi6, "N_Phi psi self cont. res"], 
    [df_sasa1, "ON_SASA" + filter_word2 + "res"], 
    [df_sasa2, "O_SASA" + filter_word2 + "res"], 
    [df_sasa3, "N_SASA" + filter_word2 + "res"], 
    [df_sasa4, "ON_SASA self cont. res"], 
    [df_sasa5, "O_SASA self cont. res"], 
    [df_sasa6, "N_SASA self cont. res"],
    ]

    workbook_path = "D:\\" + excel_filename
    with pd.ExcelWriter(workbook_path, engine = "xlsxwriter") as writer:
        for df_name, worksheet_name in list_of_worksheets:
            df_name.to_excel(writer, sheet_name = worksheet_name, index = False, na_rep = "NA")
    range_to_table.range_to_table(workbook_path)
