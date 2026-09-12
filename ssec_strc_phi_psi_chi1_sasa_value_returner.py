import os
import re
import subprocess

def file_generator(path_molprobity_files, path_molprobity_files2, path_stride_files2, path_stride_files1):
	for pdb_file in os.listdir(path_molprobity_files):
		print("processing", pdb_file)
		output_file_name = pdb_file[0:len(pdb_file) - 4] + "_" + "stride.txt"
		p1 = subprocess.Popen(["C:\\Windows\\System32\\wsl.exe", "-d", "Ubuntu-18.04", "stride", "-f" + path_stride_files2 + output_file_name, path_molprobity_files2 + pdb_file], stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = False) #whitespaces are automatically added after each comma in the list
		out, err = p1.communicate()
		# print(out.decode("utf-8"))
		# print(err.decode("utf-8"))

	list_pdb    = set([pdb_file[0:len(pdb_file) - 4] for pdb_file in os.listdir(path_molprobity_files)])
	list_stride = set([stride_file[0:len(stride_file) - 11] for stride_file in os.listdir(path_stride_files1)])

	if len(list_pdb) == len(list_stride):
		print("All pdb files processed successfully.")
	if len(list_pdb) != len(list_stride):
		failed_files = list_pdb.difference(list_stride)
		print(len(failed_files), "PDB files weren't processed:", failed_files)

def sec_str_value_returner(desired_donor_atom, desired_acceptor_atom, pdb_id_molprobity, molprobity_files): # the variable "excel_filename" was not available at this stage so I had to use "molprobity_files" as the discriminating factor between tm and globular datasets
	if "Project 03" in molprobity_files:
		path_stride_files = "D:\\PhD Work\\Project 03\\files_stride_tm\\"
	if "Project 03" not in molprobity_files:
		path_stride_files = "D:\\PhD Work\\Project 02 Selenium H-bonds & chalcogen bonds\\20210302\\files_stride_globular\\"

	stride_filename   = next((stride_file for stride_file in os.listdir(path_stride_files) if stride_file[0:4].upper() == pdb_id_molprobity[0:4].upper()), None)
	if stride_filename == None:
		return False
	if stride_filename != None:
		with open(path_stride_files + stride_filename) as file1:
			list_stridelines = file1.readlines()
			file1.close()
			dict_ssec = {"Donor ssec" : "Residue entry not found", "Acceptor ssec" : "Residue entry not found"} #if a donor/acceptor residue's entry is not found in its parent stride file, its ssec value will become equal to "Residue entry not found"
			for stride_line in list_stridelines:
				match1 = re.search(r"(ASG)(\s{2})(\w{3})(\s|\w{2}\s)(\w)(\s+)(-?\d+)(\w?\s+)(\d+)(\s+)(\w)(\s+)(\w+)(.+)", stride_line)
				if match1:
					if (match1.group(3) == desired_donor_atom.get_parent().get_resname() and 
						match1.group(5) == str(desired_donor_atom.get_parent().get_full_id()[2]) and 
						int(match1.group(7)) == desired_donor_atom.get_parent().get_full_id()[3][1]):
						dict_ssec["Donor ssec"] = match1.group(13)

					if (match1.group(3) == desired_acceptor_atom.get_parent().get_resname() and 
						match1.group(5) == str(desired_acceptor_atom.get_parent().get_full_id()[2]) and 
						int(match1.group(7)) == desired_acceptor_atom.get_parent().get_full_id()[3][1]):
						dict_ssec["Acceptor ssec"] = match1.group(13)
			return dict_ssec

def chimera_phi_psi_chi1_value_returner(desired_donor_atom, desired_acceptor_atom, pdb_id_molprobity, molprobity_files): # the variable "excel_filename" was not available at this stage so I had to use "molprobity_files" as the discriminating factor
	#^################################################################################
	#^ Define the location of the directory in which phi psi chi1 files are located 
	#^################################################################################
	if "Project 03" in molprobity_files:
		path_phi_psi_chi1 = "D:\\PhD Work\\Project 03\\files_phi_psi_chi1_tm\\"
	if "Project 02" in molprobity_files:
		path_phi_psi_chi1 = "D:\\PhD Work\\Project 02 Selenium H-bonds & chalcogen bonds\\20210302\\files_phi_psi_chi1_globular\\"

	dict_phi_psi_chi1 = dict()
	for phi_psi_chi1_file_name in os.listdir(path_phi_psi_chi1):
		match_name = re.search(r"(\w{4})(_)(phi|psi|chi1)(\.defattr)", phi_psi_chi1_file_name) # this regex will find 3 files per PDB file: "1ATZ_phi", "1ATZ_psi", "1ATZ_chi1"
		if match_name:
			if pdb_id_molprobity[0:4].upper() == match_name.group(1).upper(): # this condition will become true three time, once for "1ATZ_phi", "1ATZ_psi", "1ATZ_chi1"
				colname = match_name.group(3).title() # "match_name.group(3)" will be phi/psi/chi1
				dict_phi_psi_chi1.setdefault(colname, {"donor" : "None", "acceptor" : "None"}) 
				# when "1ATZ_phi" file will be found, the dictionary will become {"phi" : {"donor" : "None", "acceptor" : "None"}}
				# when "1ATZ_psi" file will be found, the dictionary will become {"phi" : {"donor" : "None", "acceptor" : "None"}, {"psi" : {"donor" : "None", "acceptor" : "None"}}
				# when "1ATZ_chi1" file will be found, the dictionary will become {"phi" : {"donor" : "None", "acceptor" : "None"}, {"psi" : {"donor" : "None", "acceptor" : "None"}, {"chi1" : {"donor" : "None", "acceptor" : "None"}}

				with open(path_phi_psi_chi1 + phi_psi_chi1_file_name) as file_phi_psi_chi1:
					list_phi_psi_chi1 = file_phi_psi_chi1.readlines()
					file_phi_psi_chi1.close()
					for line in list_phi_psi_chi1:
						match_phi_psi_chi1 = re.search(r"(\t/)(\w+)(:)(-?\d+\w?)(\s+)(None|-?\d+\.\d+)", line)
						if match_phi_psi_chi1:
							if (match_phi_psi_chi1.group(4) == str(desired_donor_atom.get_parent().get_full_id()[3][1]) and #resno
								match_phi_psi_chi1.group(2) == str(desired_donor_atom.get_parent().get_full_id()[2])):      #reschain
								if match_phi_psi_chi1.group(6) != "None": # if the angle value is "None" (e.g., chi1 of GLY), its "chi1" value will be "None" in the dictionary (it is pre-filled)
									dict_phi_psi_chi1[colname]["donor"] = round(float(match_phi_psi_chi1.group(6)), 2)
							if (match_phi_psi_chi1.group(4) == str(desired_acceptor_atom.get_parent().get_full_id()[3][1]) and
								match_phi_psi_chi1.group(2) == str(desired_acceptor_atom.get_parent().get_full_id()[2])):
								if match_phi_psi_chi1.group(6) != "None":
									dict_phi_psi_chi1[colname]["acceptor"] = round(float(match_phi_psi_chi1.group(6)), 2)

	return dict_phi_psi_chi1

def sasa_value_returner(desired_donor_atom, desired_acceptor_atom, pdb_id_molprobity, molprobity_files):
	if "Project 03" in molprobity_files:
		path_naccess = "D:\\PhD Work\\Project 03\\files_naccess_tm\\"
	if "Project 02" in molprobity_files:
		path_naccess = "D:\\PhD Work\\Project 02 Selenium H-bonds & chalcogen bonds\\20210302\\files_naccess_globular\\"

		dict_sasa = {"SASA donor (Å)^2" : "None", "SASA acceptor (Å)^2" : "None"}
		for naccess_file_name in os.listdir(path_naccess):
			if naccess_file_name[-4:] == ".rsa" and pdb_id_molprobity[0:4].upper() == naccess_file_name[0:4].upper():
				with open(path_naccess + naccess_file_name) as file_naccess:
					list_naccess = file_naccess.readlines()
					file_naccess.close()
					if len(list_naccess) == 0: # some Naccess rsa files are empty
						return dict_sasa
					if len(list_naccess) != 0:
						for line in list_naccess:
							match_naccess = re.search(r"(RES)(\s)(\w{3})(\s)(\w)(\s*)(-?\d+)(\w?\s+)(\d+\.\d+)(\s+)(.+)", line) #this regex takes into consideration negative resnos (2OFZ), resno greater than 1000 (1C0P) and insertion codes (3BEI)
							if match_naccess:
								if str(match_naccess.group(7)) == str(desired_donor_atom.get_parent().get_full_id()[3][1]) and str(match_naccess.group(5)) == str(desired_donor_atom.get_parent().get_full_id()[2]): # match resno & reschain
									dict_sasa["SASA donor (Å)^2"] = float(match_naccess.group(9))
								if str(match_naccess.group(7)) == str(desired_acceptor_atom.get_parent().get_full_id()[3][1]) and str(match_naccess.group(5)) == str(desired_acceptor_atom.get_parent().get_full_id()[2]): # match resno & reschain
									dict_sasa["SASA acceptor (Å)^2"] = float(match_naccess.group(9))
	return dict_sasa
