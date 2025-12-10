import os
import numpy as np
import pickle
import uuid
import shutil
import gzip
import tarfile
import sys
import trimUtils

#Incident ion
if len(sys.argv)>1:
  ion_name=sys.argv[1]
else:
  ion_name = "19F"
energies = np.arange(0.001,2,0.001) #keV
nps = 20000

#Target description
target_name = "LiF"

if __name__ == "__main__":
    # Create a unique working directory for this ion run, cp srim to that folder
    unique_id = uuid.uuid4().hex[:8]
    temp_workdir = os.path.join("C:/Users/Sam/Desktop/SRIM_tmp", f"{ion_name}_{unique_id}")
    shutil.copytree("C:/Users/Sam/Desktop/SRIM_exe", temp_workdir)
    srimFolder = temp_workdir + "/"

    #Handle output:
    output_energies = []
    output_efficiencies = []

    # Create new SRIM input
    os.chdir(srimFolder)
    for energy in energies:
        output_energies.append(energy)

        lines = trimUtils.makeTrimInputString(energy,target_name,ion_name,nps)
        with open(f"{srimFolder}TRIM.in", "w") as trimFile:
            for line in lines:
                if not line.endswith("\n"):
                    line += "\n"
                trimFile.write(line)

        # Run SRIM  
        exit_code = os.system("TRIM.exe")

        # Upon successful completion, mv collison.txt to file name, compress, remove temp folder
        if exit_code == 0:
            #Output path
            src = os.path.join(srimFolder, "SRIM Outputs", "COLLISON.txt")

            #Parse output, essentially just count # of "Summary of Above Cascade"
            nDamagelessPrimaries = 0
            with open(src, "r") as f:
                for line in f:
                    # Every completed cascade ends with "Summary of Above Cascade"
                    if "Vacancies         = 000000.0" in line:
                        nDamagelessPrimaries += 1

            #Calculate fraction of # damage producing recoils / nps
            eff = float(nps-nDamagelessPrimaries)/float(nps)

            #Append to CSV
            output_efficiencies.append(eff)
            os.remove(src)

            print(f"TRIM run completed successfully for E={energy} keV, eff={eff}")

            if eff >= 1.:    
              print(f"Efficiency saturated (>0.999) at {energy:.3f} keV — stopping sweep.")
              break

    outFilename = f"C:\\Users\\Sam\\Documents\\code\\trimRunner\\data\\eff\\{target_name}_{ion_name}.csv"
    with open(outFilename,"w") as outFile:
        line="Ion energy (keV),efficiency\n"
        outFile.write(line)
        for i in range(0,len(output_energies)):
            line=f"{output_energies[i]:.4f},{output_efficiencies[i]:.5f}\n"
            outFile.write(line)
    #shutil.rmtree(temp_workdir, ignore_errors=True)
