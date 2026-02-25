import os
import numpy as np
import pickle
import uuid
import shutil
import gzip
import tarfile
import trimUtils
import sys

if len(sys.argv) < 2:
  print("Error! Pass in config!\n")
  sys.exit()

inpFilename = sys.argv[1]
mass_dict = trimUtils.createMassDict("mass_1.mas20.txt")
materials_dict = trimUtils.createMaterialsDict()
input_data = trimUtils.parseConfig(inpFilename)
trimUtils.checkTrimArgs(input_data,mass_dict,materials_dict)

SRIM_EXE_PATH= "C:/Users/Sam/Desktop/SRIM_exe"
SRIM_TMP_PATH = "C:/Users/Sam/Desktop/SRIM_tmp"

target_name = input_data["material"]
ion_name = input_data["ionSymbol"]
outputFolder = input_data["outputPath"]
runMode = input_data["runMode"]
energy = input_data["energy_keV"]
nps = input_data["nps"]

output_base_name = f"{target_name}_{ion_name}_{energy}keV"

# Create new SRIM run folder, go there
runFolder = trimUtils.makeTempSRIMFolder(f"{target_name}_{ion_name}",SRIM_TMP_PATH,SRIM_EXE_PATH)
try:
  os.chdir(runFolder)

  if runMode=="damage":
    txt_path = os.path.join(outputFolder, f"{output_base_name}.txt")
    tar_path = os.path.join(outputFolder, f"{output_base_name}.tar.gz")
    if os.path.exists(txt_path):
      os.remove(txt_path)
    if os.path.exists(tar_path):
      os.remove(tar_path)

    # Run SRIM return location of collison.txt [sic] 
    collison_path = trimUtils.runSRIM(runFolder,input_data,mass_dict,materials_dict)

    # Move and compress output
    os.replace(collison_path, txt_path)

    # Create .tar.gz
    with tarfile.open(tar_path, "w:gz") as tar:
      tar.add(txt_path, arcname=os.path.basename(txt_path))
    os.remove(txt_path)

  elif runMode=="efficiency":
    csv_path = os.path.join(outputFolder, f"{output_base_name}.csv")
    energies = np.arange(0.001,energy,0.001) #keV
    output_efficiencies = []

    for energy in energies:
      input_data["energy_keV"] = energy
      # Run SRIM, return location of collison.txt [sic] 
      collison_path = trimUtils.runSRIM(runFolder,input_data,mass_dict,materials_dict)

      #Parse output, essentially just count # of "Summary of Above Cascade"
      nDamagelessPrimaries = 0
      with open(collison_path, "r") as f:
        for line in f:
          # Every completed cascade ends with "Summary of Above Cascade"
          if "Vacancies         = 000000.0" in line:
            nDamagelessPrimaries += 1

      #Calculate fraction of # damage producing recoils / nps
      eff = float(nps-nDamagelessPrimaries)/float(nps)

      #Append to CSV
      output_efficiencies.append(eff)
      os.remove(collison_path)
      print(f"TRIM run completed successfully for E={energy} keV, eff={eff}")
      if eff >= 1.:    
        print(f"Efficiency saturated (>0.999) at {energy:.3f} keV — stopping sweep.")
        break
    
    with open(csv_path,"w") as outFile:
      line="Ion energy (keV),efficiency\n"
      outFile.write(line)
      for i in range(0,len(output_efficiencies)):
        line=f"{energies[i]:.4f},{output_efficiencies[i]:.5f}\n"
        outFile.write(line)
finally:
  os.chdir(SRIM_TMP_PATH)
  shutil.rmtree(runFolder, ignore_errors=True)
