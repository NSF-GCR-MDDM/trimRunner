import os
import numpy as np
import pickle
import uuid
import shutil
import gzip
import tarfile
import trimUtils
import sys

#Incident ion
if len(sys.argv) >= 2:
  ion_name = sys.argv[1]
else:
  ion_name = "19F"
energy = 10 #keV
if len(sys.argv) >= 3:
  nps = int(sys.argv[2])
else:
  nps = 2000

#Target description
target_name = "Diamond"

if __name__ == "__main__":
  # Create a unique working directory for this ion run, cp srim to that folder
  unique_id = uuid.uuid4().hex[:8]
  temp_workdir = os.path.join("C:/Users/Sam/Desktop/SRIM_tmp", f"{ion_name}_{unique_id}")
  shutil.copytree("C:/Users/Sam/Desktop/SRIM_exe", temp_workdir)
  srimFolder = temp_workdir + "/"

  # Make output folder, output file names
  outputFolder = f"C:/Users/Sam/Documents/code/trimRunner/data/txt/{target_name}/"
  if not os.path.exists(outputFolder):
    os.mkdir(outputFolder)

  txt_name = "{0}_{1}.txt".format(target_name, ion_name)
  txt_path = os.path.join(outputFolder, txt_name)

  # Overwrite if exists
  if os.path.exists(txt_path):
    os.remove(txt_path)
  tar_path =  "{0}_{1}.tar.gz".format(target_name, ion_name)
  if os.path.exists(tar_path):
    os.remove(tar_path)

  # Create new SRIM input
  os.chdir(srimFolder)
  lines = trimUtils.makeTrimInputString(energy,target_name,ion_name,nps)
  with open(f"{srimFolder}TRIM.in", "w") as trimFile:
    for line in lines:
      if not line.endswith("\n"):
        line += "\n"
      trimFile.write(line)
    trimFile.close()

  # Run SRIM  
  exit_code = os.system("TRIM.exe")

  # Upon successful completion, mv collison.txt to file name, compress, remove temp folder
  if exit_code == 0:
    print("TRIM completed successfully.")

    # Move and compress output
    src = os.path.join(srimFolder, "SRIM Outputs", "COLLISON.txt")
    os.replace(src, txt_path)

    # Create .tar.gz
    with tarfile.open(tar_path, "w:gz") as tar:
      tar.add(txt_path, arcname=os.path.basename(txt_path))
      tar.close()

