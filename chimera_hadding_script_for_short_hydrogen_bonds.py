#NOTE: This script is called by the short hydrogen bond systems creating script of VS code
import os
from chimera import runCommand

path_unhadded_files = "D:\\unhadded_files\\"
for file in os.listdir(path_unhadded_files):
    runCommand("open " + path_unhadded_files + file)
    runCommand("addh")
    hadded_name = file[0 : len(file) - 4] + "_hadded.pdb"

    hadded_folder = "D:\\hadded_files\\"
    if not os.path.isdir(hadded_folder):
        os.makedirs(hadded_folder)

    runCommand("write 0 " + hadded_folder + hadded_name)
    runCommand("close 0")

runCommand("stop") # exits chimera